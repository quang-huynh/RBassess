## ----setup, include = FALSE, echo = FALSE--------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(RBassess)
data(BC_lakes)
data_object <- BC_lakes[[69]]


## ---- eval = FALSE-------------------------------------------------------
#  data_object <- new("RBdata")

## ---- eval = FALSE-------------------------------------------------------
#  data(BC_lakes)
#  data_object <- BC_lakes[[69]]
#  plot(data_object)

## ------------------------------------------------------------------------
model <- fit_model(data_object)

## ------------------------------------------------------------------------
summary(model)

## ---- eval = FALSE-------------------------------------------------------
#  plot(model)

## ----eval = FALSE--------------------------------------------------------
#  samps <- run_mcmc(model)
#  summary(samps)

## ----echo = FALSE--------------------------------------------------------
res <- read.csv('mcmc_summary.csv')

res_show <- res[, -1]
rownames(res_show) <- res[, 1]
round(res_show, 2)

## ----eval = FALSE--------------------------------------------------------
#  rstan::stan_trace(samps)

## ----eval = FALSE--------------------------------------------------------
#  samps_prior <- run_mcmc(model, priors_only = TRUE)
#  plot(samps, samps_prior)

## ----eval = FALSE--------------------------------------------------------
#  report(model)
#  report(samps)

