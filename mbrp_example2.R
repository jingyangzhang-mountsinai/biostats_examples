# Author: Jingyang (Judy) Zhang
# Date: Nov. 7th, 2025
# Purpose: an example of applying model-based recursive partitioning (mbrp).
library(partykit)
library(TH.data)
library(survival)
# German women with positive node breast cancer. 
data("GBSG2", package = "TH.data")

# For regression, a parametric Weibull regression based on the `survreg()` function in `survival` package is used.
## Set up fitting function for `mob()`:
wbreg <- function(y, x, start = NULL, weights = NULL, offset = NULL, ...){
  survreg(y ~ 0 + x, weights = weights, dist = "weibull")
}

# As the `survival` package does not provide a `logLik()` method for `survreg` objects, this needs to be defined. 
logLik.survreg <- function(object,...){
  structure(object$loglik[2], df = sum(object$df), class = "logLik")

}


gbsg2_tree <- mob(Surv(time, cens) ~ horTh + pnodes | age + tsize + tgrade + progrec + estrec + menostat, data = GBSG2, 
                  fit = wbreg, control = mob_control(minsize = 80))




gbsg2_tree


plot(gbsg2_tree)
