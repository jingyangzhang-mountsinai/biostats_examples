# Author: Jingyang (Judy) Zhang
# Date: Nov. 12th, 2025
# Purpose: an example of building causal survival forest targeting a Restricted Mean Survival Time (RMST) with maximum follow-up time set to `horizon`.
# Source: https://grf-labs.github.io/grf/articles/rate.html

# An Application of RATE Metric for Treatment Prioritization Rules to SPRINT and ACCORD
## Target Estimand: the difference in restricted mean survival time (RMST) conditional on covariates:
## tau(x) = E[T(1)^h - T(0)^h|X=x]
### T: the (censored) survival time.
### h: horizon, set to 3 years (because after 3 years the SPRINT trial data nearly stops exhibiting events).
### W.hat = E[W_i|X_i = x]: propensity scores.
#### Since the treatment was randomized, set the propensity score to the mean number of treated units 



# Read in semi-synthetic data from https://github.com/grf-labs/grf/tree/master/r-package/grf/vignettes.
library(tidyverse)
library(grf)
load("data/synthetic_SPRINT_ACCORD.RData")

df <- data.frame(Y = c(Y.sprint, Y.accord),
                 D = c(D.sprint, D.accord),
                 data = c(rep("synthetic-SPRINT", length(Y.sprint)),
                          rep("synthetic-ACCORD", length(Y.accord))))


df$Censored <- factor(df$D, labels = c("Yes", "No"))

ggplot(df, aes(x = Y, fill = Censored)) + 
  facet_wrap(data ~ .) + 
  geom_histogram(alapha = 0.5, bins = 30) + 
  xlab("Time Until Primary Outcome (days)") + 
  ylab("Frequency") + 
  theme_classic()

# 1. The SPRINT trial was halted early due to a lower number of events occurring -> very high censoring rate. 



# Use causal survival forest to estimated CATE. 
horizon <- 3 * 365

csf.sprint <- causal_survival_forest(X.sprint, Y.sprint, W.sprint, D.sprint,
                                     W.hat = mean(W.sprint),
                                     target = "survival.probability",
                                     horizon = horizon)


csf.accord <- causal_survival_forest(X.accord, Y.accord, W.accord, D.accord,
                                     W.hat = mean(W.accord),
                                     target = "survival.probability",
                                     horizon = horizon)


tau.hat.sprint <- predict(csf.accord, X.sprint)$predictions
tau.hat.accord <- predict(csf.sprint, X.accord)$predictions

# Estimate RATE on SPRINT and ACCORD, using estimated CATE functions from ACCORD and SPRINT. 

rate.sprint <- rank_average_treatment_effect(csf.sprint, tau.hat.sprint, target = "AUTOC")
rate.accord <- rank_average_treatment_effect(csf.accord, tau.hat.accord, target = "AUTOC")


rate.sprint

rate.accord



# Visualize the AUTOC.
par(mfrow = c(1, 2))
plot(rate.sprint, xlab = "Treated fraction", main = "TOC evaluated on SPRINT\n tau(X) estimated from ACCORD")
plot(rate.accord, xlab = "Treated fraction", main = "TOC evaluated on ACCORD\n tau(X) estimated from SPRINT")


# Interpretation:
## 1. In this semi-synthetic example, both AUTOCs are insignificant at conventional levels -> no evidence of significant HTEs in the two trials. 
### This can be attributed to 1) low power (perhaps the sample size is not large enough to detect HTEs); 2) HTE estimator does not detect HTE; 3) HTE along observable predictors are truly negligible. 