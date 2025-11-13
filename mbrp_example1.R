# Author: Jingyang (Judy) Zhang
# Date: Nov. 7th, 2025
# Purpose: an example of applying model-based recursive partitioning (mbrp).
library(partykit)
if(require("mlbench")){
  
  ## Pima Indians diabetes data. 
  data("PimaIndiansDiabetes", package = "mlbench")
  
  ## A simple basic fitting function (of type 1) for a logitic regression. 
  logit <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...){
    glm(y ~ 0 + x, family = binomial, start = start, ...)
  }


# Set up a logistic regression tree.

pid_tree <- mob(diabetes ~ glucose | pregnant + pressure + triceps + insulin + mass + pedigree + age, 
                data = PimaIndiansDiabetes, fit = logit)

## Print tree
print(pid_tree)

## Print info about (some) nodes.
print(pid_tree, node = 3:4)

## Visualization 
plot(pid_tree)

## Coefficients and summary.
coef(pid_tree)
coef(pid_tree, node = 1)
summary(pid_tree, node = 1)


## Compare: normal logistic regression
logistic <- glm(diabetes ~ glucose + pregnant + pressure + triceps + insulin + mass + pedigree + age, data = PimaIndiansDiabetes, family = binomial)
summary(logistic)

## Average deviance computed in different ways:

### Method 1
mean(residuals(pid_tree)^2)

### Method 2
deviance(pid_tree)/sum(weights(pid_tree))


### Method 3
deviance(pid_tree)/nobs(pid_tree)


## log-likelihood and info criteria.
logLik(pid_tree)
AIC(pid_tree)
BIC(pid_tree)

## Predicted nodes.
predict(pid_tree, newdata = head(PimaIndiansDiabetes, 6), type = "node")
### Other types of predictions are possible using lmtree() or glmtree().

}