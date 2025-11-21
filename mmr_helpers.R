# Author: Jingyang (Judy) Zhang
# Date: Nov 12th, 2025
# A set of helper functions.


## Function: find variables in a dataframe that has no variation. 

find_novar <- function(df){
  no_variation <- sapply(df, function(x) length(unique(x[!is.na(x)])) == 1)
  novar_vars <- names(df)[no_variation]
  return(novar_vars) 
}



# Function: find variables in a dataframe that has little variation. 
find_low_variation <- function(df, threshold = 0.95){
  low_var <- sapply(df, function(x){
    x_nonNA <- x[!is.na(x)]
    if(length(x_nonNA) <= 1){return(TRUE)}
    
    
    if (is.numeric(x_nonNA)){
      
      
      # For numeric: use relative variance (sd/range).
      sd_x <- sd(x_nonNA)
      range_x <- max(x_nonNA) - min(x_nonNA)
      
      if(range_x == 0){return(TRUE)} # x is constant.
      (sd_x / range_x) < (1-threshold) # x has low variation.
    }else{
      # For categorical/binary: use proportion of most frequently level.
      freq <- table(x_nonNA)
      max(freq)/sum(freq) >= threshold
    }
  })
  return(low_var)
}



## Function: check if a variable is factor (when it is likely coded as 0, 1 and being treated as numeric by R).

check_factor <- function(var, max_unique_level = 5){
    if(is.numeric(var)){
      ### Count number of unique non-NA values. 
      n_unique <- length(unique(na.omit(var)))
      
      ### Return TRUE if the number of unique levels in the variable is less than the max_unique_level.
      if(n_unique <= max_unique_level){
        return(TRUE)
      }else
        {return(FALSE)}
    }
}
  
  

# Function: convert numeric-coded categorical/binary variables in a dataframe to factor.
convert_to_factor <- function(df, max_unique = 5){
  for(col_name in names(df)){
    if(is.numeric(df[[col_name]])){
 
      if(check_factor(df[[col_name]])){
        print(col_name)
        df[[col_name]] <- factor(df[[col_name]])
      }
    }
  }
  
  return(df)
}



# Function: calculate the hazard ratio of treatment. 
cal_hr_trt <- function(df, Y, D, W, X){
  # Build formula as a string
  fmla <- as.formula(paste0("Surv(", Y, ", ", D, ") ~ ", W, " + ", paste(X, collapse = " + ")))
  
  # Fit Cox model
  cox_fit <- coxph(fmla, data = df)
  
  # Extract hazard ratio for treatment variable
  hr_trt <- round(exp(cox_fit$coefficients[names(cox_fit$coefficients) == W]), 2)
  
  return(unname(hr_trt))
}



# Function: calculates the proportion of events.
cal_event_prop <- function(df, D, W, trt = 0){
  num_events <- table(df[[D]][df[[W]] == trt])["1"]
  num_events <- ifelse(is.na(num_events), 0, unname(num_events))
  
  num_group <- nrow(df[df[[W]] == trt, ])
  
  prop_text <- paste(num_events, "/", num_group)
  return(prop_text)
}



# Function: plot survival curves by treatment and annotate with observed event rate by group.

plot_survival <- function(df, Y, D, W, ITEQ, cluster = FALSE){
  
  ctrl_prop <- cal_event_prop(df, D, W, trt = 0)
  trt_prop <- cal_event_prop(df, D, W, trt = 1)
  if (cluster){
    title = paste("Survival Probability for Estimated ITEs in Cluster", ITEQ)
  }else{
    
    title = paste("Survival Probability for Estimated ITEs in the Quantile", ITEQ)
    
    
  }
  # Build formula as a string
  fmla <- as.formula(paste0("Surv(", Y, ", ", D, ") ~ ", W))
  survfit2(fmla, data = df) %>% ggsurvfit() + 
    scale_color_manual(
      name = "Treatment Group",
      values = c("#4b5663", "#ff6666"), 
      labels = c("CABG Alone", "CABG + MV Repair")
    ) + 
    labs(title = title,
         y = "Survival Probability") + 
    annotate("text", x = 0, y = c(0, 0.1),
             label = c(paste("Alone:", ctrl_prop), paste("Repair:", trt_prop)),
             hjust = 0,  # left align
             size = 3) + 
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(size = 7),
      axis.title = element_text(size = 6),
      legend.title = element_text(size = 6),
      legend.text = element_text(size = 6),
      plot.margin = margin(5, 5, 5, 5)
    )
  
}



# Function: extract predictors names from a dataframe.

get_pred_names <- function(df, Y, D, W, others){

  X <- df %>% dplyr::select(-all_of(c(Y, D, W, others)))
  
  return(colnames(X))
}



# Function: annotation function for MOB tree. Splits label for inner nodes. Shows variable name and cutpoint. 

inner_fun <- function(node) {
  split <- node$split
  varname <- varid_to_name(mob_fit)[split$varid]
  cut <- split$breaks
  
  paste0(varname, " < ", round(cut, 3))
}


# Function: annotation function for MOB tree. Shows node sample size (n), intercept + treatment coefficient, R-sqred of local model. 
terminal_fun <- function(node) {
  n <- node$n
  m <- node$model
  
  cf <- coef(m)
  
  # Example: show first and second coefficients
  label <- paste0(
    "n = ", n,
    "\nIntercept = ", round(cf[1], 3),
    "\nCoef2 = ", round(cf[2], 3)
  )
  
  # You can compute anything: p-values, R2, treatment effects, etc.
  label
}