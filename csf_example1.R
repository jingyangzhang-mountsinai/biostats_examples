# Save this as causal_survival_forest_simulation.R and run in R (>= 4.0)
# Required packages:
# install.packages(c("randomForestSRC","survival","pec","dplyr","data.table"))
library(randomForestSRC)
library(survival)
library(pec)
library(dplyr)
library(data.table)

set.seed(2025)

n <- 254
p_total <- 51   # final column count target

# 1) Treatment assignment (randomized ~50/50)
randomization_assignment <- rbinom(n, 1, 0.5)

# 2) Generate 11 numeric covariates
num_cov <- replicate(11, rnorm(n))
colnames(num_cov) <- paste0("num_", 1:11)
num_cov <- as.data.frame(num_cov)

# 3) Generate 37 categorical covariates (some binary, some 3-5 levels)
cat_cov <- list()
for (i in 1:37) {
  if (i <= 10) k <- 2
  else if (i <= 20) k <- 3
  else if (i <= 30) k <- 4
  else k <- 5
  levs <- paste0("cat", i, "_lvl", 1:k)
  cat_cov[[paste0("cat_", i)]] <- sample(levs, size = n, replace = TRUE)
}
cat_cov <- as.data.frame(cat_cov, stringsAsFactors = FALSE)

# 4) Assemble df and impose missingness on many covariates (but not on a few)
df <- cbind(randomization_assignment = randomization_assignment, num_cov, cat_cov)

# Always-observed covariates (no missingness)
always_observed <- c("num_1", "num_2", "cat_1", "cat_2", "cat_3")

# introduce missingness at random in others (5% - 35% missing)
for (col in setdiff(names(df), always_observed)) {
  if (col == "randomization_assignment") next
  miss_rate <- runif(1, 0.05, 0.35)
  is_na <- runif(n) < miss_rate
  df[is_na, col] <- NA
}

# 5) Simulate event times influenced by numerics and treatment
# linear predictor from numeric covariates
beta_nums <- c(0.3, -0.2, 0.1, 0.0, 0.05, -0.1, 0.15, 0, 0, 0.02, -0.05)
lp <- as.matrix(sapply(paste0("num_", 1:11), function(col) {
  x <- as.numeric(df[[col]])
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
})) %*% beta_nums

# treatment protective effect
treatment_coef <- -0.5
lp <- lp + treatment_coef * df$randomization_assignment

# baseline exponential event times (days)
scale <- 365.0
event_time <- rexp(n, rate = 1 / (scale * exp(-lp)))  # mean depends on lp

# independent censoring times adjusted to target ~70% missing in first_mace_day
# We'll find a censor scale that yields ~70% censored (so first_mace_day missing ~70%)
get_miss_prop <- function(censor_scale) {
  censor_time <- rexp(n, rate = 1 / censor_scale)
  obs_time <- pmin(event_time, censor_time)
  event_ind <- as.integer(event_time <= censor_time)
  first_mace_day <- ifelse(event_ind == 1, round(event_time), NA_real_)
  mean(is.na(first_mace_day))
}
# search for a censor_scale producing ~0.70 missing
censor_scale <- scale * 1.2
for (i in 1:10) {
  mp <- get_miss_prop(censor_scale)
  if (mp < 0.67) censor_scale <- censor_scale * 0.9 else if (mp > 0.73) censor_scale <- censor_scale * 1.1 else break
}

censor_time <- rexp(n, rate = 1 / censor_scale)
time_observed <- pmin(event_time, censor_time)
event_observed <- as.integer(event_time <= censor_time)
mace <- event_observed
first_mace_day <- ifelse(mace == 1, round(event_time), NA_real_)

cat(sprintf("Proportion missing in first_mace_day: %.3f (target ~0.70)\n", mean(is.na(first_mace_day))))

# 6) Final dataset: 51 columns (treatment, first_mace_day, mace, 48 covariates)
final_df <- df
final_df$first_mace_day <- first_mace_day
final_df$mace <- mace

stopifnot(ncol(final_df) == 51)

# Save dataset
write.csv(final_df, "synthetic_causal_survival_dataset_254x51.csv", row.names = FALSE)
cat("Saved dataset to synthetic_causal_survival_dataset_254x51.csv\n")

# 7) Preprocess for modelling:
#    - For randomForestSRC, convert factors and impute simple missingness
#    - We'll impute numeric NAs with median and convert categorical NA to "MISSING" level.
proc_df <- final_df
num_cols <- paste0("num_", 1:11)
for (c in num_cols) {
  med <- median(as.numeric(proc_df[[c]]), na.rm = TRUE)
  proc_df[[c]][is.na(proc_df[[c]])] <- med
}
cat_cols <- setdiff(names(final_df), c("randomization_assignment","first_mace_day","mace", num_cols))
for (c in cat_cols) {
  proc_df[[c]] <- as.factor(as.character(proc_df[[c]]))
  # make "NA" into a factor level "MISSING"
  levels(proc_df[[c]]) <- c(levels(proc_df[[c]]), "MISSING")
  proc_df[[c]][is.na(proc_df[[c]])] <- "MISSING"
}

# 8) Create a model matrix of covariates (exclude first_mace_day & mace)
X_vars <- setdiff(names(proc_df), c("first_mace_day", "mace"))
# We'll use the variables directly (randomForestSRC handles factors)
model_df <- proc_df[, c(X_vars, "mace")]
# But model_df needs the outcome as a Surv object; we use the internal time_observed and event_observed
model_df$time_observed <- time_observed
model_df$event_observed <- event_observed

# 9) CV + grid search
# We'll do 5-fold CV. For each fold, fit two RSFs (treated & control) on training fold,
# predict survival probabilities at a horizon on validation fold, and compute concordance index.
k <- 5
set.seed(2025)
folds <- sample(rep(1:k, length.out = n))

param_grid <- expand.grid(ntree = c(100, 200),
                          nodesize = c(3, 10),
                          mtry = c( floor(sqrt(length(X_vars)-1)), floor((length(X_vars)-1)/3) ),
                          stringsAsFactors = FALSE)

cv_results <- list()

horizon <- 365  # days for survival prob comparison

for (i in 1:nrow(param_grid)) {
  params <- param_grid[i, ]
  cis <- c()
  for (fold in 1:k) {
    train_idx <- which(folds != fold)
    val_idx   <- which(folds == fold)
    
    df_train <- model_df[train_idx, ]
    df_val   <- model_df[val_idx, ]
    
    # Split training by treatment
    train_treated <- df_train[df_train$randomization_assignment == 1, ]
    train_control <- df_train[df_train$randomization_assignment == 0, ]
    
    # require a minimum number of events in group; else skip fold
    if (sum(train_treated$event_observed) < 8 || sum(train_control$event_observed) < 8) {
      cis <- c(cis, NA)
      next
    }
    
    # Fit RSF for treated and control
    rf_t <- rfsrc(Surv(time_observed, event_observed) ~ . - randomization_assignment - time_observed - event_observed - mace,
                  data = train_treated,
                  ntree = params$ntree, nodesize = params$nodesize, mtry = params$mtry,
                  nsplit = 10, importance = "none", splitrule = "logrank")
    rf_c <- rfsrc(Surv(time_observed, event_observed) ~ . - randomization_assignment - time_observed - event_observed - mace,
                  data = train_control,
                  ntree = params$ntree, nodesize = params$nodesize, mtry = params$mtry,
                  nsplit = 10, importance = "none", splitrule = "logrank")
    
    # Predict survival probability at horizon for validation subjects using each model
    # randomForestSRC predict returns a matrix of survival probabilities (time × n)
    pred_t <- predict(rf_t, newdata = df_val)$survival  # matrix: n_timepoints x n_val
    pred_c <- predict(rf_c, newdata = df_val)$survival
    
    times_t <- predict(rf_t, newdata = df_val)$time.interest
    times_c <- predict(rf_c, newdata = df_val)$time.interest
    
    # function to evaluate survival prob at horizon from survival matrix
    eval_at_horizon <- function(surv_mat, times, horizon) {
      idx <- max(which(times <= horizon))
      if (is.infinite(idx) || length(idx) == 0) {
        # if horizon smaller than smallest time, return 1s
        return(rep(1, ncol(surv_mat)))
      } else {
        return(as.numeric(surv_mat[idx, ]))
      }
    }
    s_t_val <- eval_at_horizon(pred_t, times_t, horizon)
    s_c_val <- eval_at_horizon(pred_c, times_c, horizon)
    
    # combine predicted risk according to actual treatment in validation set
    pred_surv_by_obs <- ifelse(df_val$randomization_assignment == 1, s_t_val, s_c_val)
    # risk = 1 - survival
    risk_pred <- 1 - pred_surv_by_obs
    
    # compute c-index using pec::cindex or survival::concordance
    # pec::cindex requires a matrix of predictions; simpler: use survival::concordance
    conc <- concordance(Surv(df_val$time_observed, df_val$event_observed) ~ risk_pred)$concordance
    cis <- c(cis, conc)
  } # folds
  cv_results[[i]] <- list(params = params, mean_cindex = mean(cis, na.rm = TRUE), fold_cindices = cis)
  cat(sprintf("params %d/%d: ntree=%d nodesize=%d mtry=%d => mean C-index=%.4f\n",
              i, nrow(param_grid), params$ntree, params$nodesize, params$mtry, cv_results[[i]]$mean_cindex))
}

# Collect CV results
cv_summary <- do.call(rbind, lapply(cv_results, function(x) {
  data.frame(ntree = x$params$ntree, nodesize = x$params$nodesize, mtry = x$params$mtry, mean_cindex = x$mean_cindex)
}))
cv_summary <- cv_summary[order(-cv_summary$mean_cindex), ]
print(cv_summary)

best <- cv_results[[which.max(sapply(cv_results, function(x) x$mean_cindex))]]
cat("Best params (by mean CV C-index):\n"); print(best$params)

# 10) Fit final RSFs on full data using best params
rf_t_full <- rfsrc(Surv(time_observed, event_observed) ~ . - randomization_assignment - time_observed - event_observed - mace,
                   data = model_df[model_df$randomization_assignment == 1, ],
                   ntree = best$params$ntree, nodesize = best$params$nodesize, mtry = best$params$mtry)
rf_c_full <- rfsrc(Surv(time_observed, event_observed) ~ . - randomization_assignment - time_observed - event_observed - mace,
                   data = model_df[model_df$randomization_assignment == 0, ],
                   ntree = best$params$ntree, nodesize = best$params$nodesize, mtry = best$params$mtry)

# Predict survival at horizon for every subject under both models
pred_t_all <- predict(rf_t_full, newdata = model_df)$survival
pred_c_all <- predict(rf_c_full, newdata = model_df)$survival
times_t_all <- predict(rf_t_full, newdata = model_df)$time.interest
times_c_all <- predict(rf_c_full, newdata = model_df)$time.interest

eval_at_horizon_colwise <- function(surv_mat, times, horizon) {
  idxs <- which(times <= horizon)
  if (length(idxs) == 0) {
    return(rep(1, ncol(surv_mat)))
  }
  idx <- max(idxs)
  as.numeric(surv_mat[idx, ])
}

s_t_all <- eval_at_horizon_colwise(pred_t_all, times_t_all, horizon)
s_c_all <- eval_at_horizon_colwise(pred_c_all, times_c_all, horizon)
ite_surv <- s_t_all - s_c_all

# Attach to final_df and save
final_out <- final_df
final_out$pred_surv_if_treated_365d <- round(s_t_all, 4)
final_out$pred_surv_if_control_365d <- round(s_c_all, 4)
final_out$ite_surv_365d <- round(ite_surv, 4)

write.csv(final_out, "synthetic_causal_survival_dataset_with_ite_254x51.csv", row.names = FALSE)
cat("Saved dataset with ITEs to synthetic_causal_survival_dataset_with_ite_254x51.csv\n")

# Summaries
cat("Summary of ITE (S_treated - S_control at 365 days):\n")
print(summary(final_out$ite_surv_365d))

# Optional: simple diagnostic plot of ITE distribution
# Uncomment if you want to visualize
# hist(final_out$ite_surv_365d, breaks=30, main="Distribution of estimated ITE on 365d survival", xlab="ITE (S_treated - S_control)")

# Compute overall C-index on full data using predictions matched to observed treatment
pred_surv_obs <- ifelse(final_out$randomization_assignment == 1, final_out$pred_surv_if_treated_365d, final_out$pred_surv_if_control_365d)
risk_obs <- 1 - pred_surv_obs
conc_full <- concordance(Surv(time_observed, event_observed) ~ risk_obs)$concordance
cat(sprintf("C-index on full data (using predicted risk at horizon): %.4f\n", conc_full))
