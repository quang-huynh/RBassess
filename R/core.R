#' @name BC_lakes
#' @docType data
#' @title Dataset from 256 surveys of small lakes in British Columbia, Canada
#'
#' @examples
#' data(BC_lakes)
#'
#' plot(BC_lakes[[69]]) # Plot data
#'
#' mod <- fit_model(BC_lakes[[69]]) # Run model
#' plot(mod) # View model results
NULL

#' Fit the model
#'
#' Setup the TMB model and runs the optimization to return the posterior mode.
#'
#' @param RBdata An object of class \linkS4class{RBdata} containing data and priors for the model.
#' @param start An optional list of starting values. Otherwise, the prior means are used as starting values.
#' @param nit_harvest The number of iterations for solving F and the probability of harvesting in the model.
#' @param use_priors Logical, whether to use priors for the parameters in the model.
#' @param lower A vector of the lower bounds of the parameters. This argument overrides the default, and thus generally not recommended to use.
#' @param upper A vector of the upper bounds of the parameters. This argument overrides the default, and thus generally not recommended to use.
#' @param control A list of control arguments to be passed to \link[stats]{nlminb}.
#' @param ... More arguments to pass in.
#' @seealso \link{run_mcmc} \link{report} \link{simulation} \link{summary.RBfit} \link{plot.RBfit}
#' @useDynLib RBassess
#' @import TMB
#' @importFrom stats nlminb
#' @examples
#' data(BC_lakes)
#' mod <- fit_model(BC_lakes[[42]]) # Run model
#' plot(mod) # View model results
#' report(mod) # Generate HTML report
#' @export
fit_model <- function(RBdata, start = NULL, nit_harvest = 5L, use_priors = TRUE, lower = NULL, upper = NULL,
                      control = list(iter.max = 1e3, eval.max = 1e3), ...) {
  dots <- list(...)
  RBdata <- validate_RBdata(RBdata)

  model_par <- list()
  model_par$Linf <- if(!is.null(start$Linf)) start$Linf else RBdata@prior_Linf[1]
  model_par$K <- if(!is.null(start$K)) start$K else RBdata@prior_K[1]
  model_par$CV_Len <- if(!is.null(start$CV_Len)) start$CV_Len else RBdata@prior_CV_Len[1]
  model_par$M <- if(!is.null(start$M)) start$M else RBdata@prior_M[1]

  model_par$q <- if(!is.null(start$q)) start$q else RBdata@prior_q[1]
  model_par$Effort <- if(!is.null(start$Effort)) start$Effort else RBdata@prior_Effort[1]
  model_par$GN_SL50 <- if(!is.null(start$GN_SL50)) start$GN_SL50 else RBdata@prior_GN_SL50[1]
  model_par$GN_gamma <- if(!is.null(start$GN_gamma)) start$GN_gamma else RBdata@prior_GN_gamma[1]
  model_par$angler_SL50 <- if(!is.null(start$angler_SL50)) start$angler_SL50 else RBdata@prior_angler_SL50[1]
  model_par$angler_gamma <- if(!is.null(start$angler_gamma)) start$angler_gamma else RBdata@prior_angler_gamma[1]

  model_data <- list(Len_data = RBdata@Length, Len_age_data = RBdata@Age_length, age_adjust = RBdata@Age_adjust,
                     Lmid = RBdata@Length_bin, Lbin_width = unique(diff(RBdata@Length_bin)), stocking_density = RBdata@stock_density,
                     L_stock = RBdata@L_stock, bag_limit = RBdata@bag_limit, release_mortality = RBdata@release_mortality,
                     p_vrel = RBdata@p_vrel, use_likelihood = as.integer(TRUE), use_priors = as.integer(use_priors),
                     init_p_harvest = if(!is.null(start$p_harvest)) start$p_harvest else 0.5, nit_harvest = nit_harvest)

  model_priors <- generate_priors(RBdata)

  obj <- MakeADFun(data = c(model_data, model_priors), parameters = model_par, DLL = "RBassess", hessian = TRUE, silent = TRUE)

  if(is.null(lower)) lower <- rep(0, 10)
  if(is.null(upper)) upper <- rep(Inf, 10)

  opt <- try(nlminb(obj$par, obj$fn, obj$gr, control = control, lower = lower, upper = upper), silent = TRUE)
  SD <- if(!is.character(opt)) sdreport(obj, opt$par, obj$he(opt$par), getReportCovariance = FALSE) else sdreport(obj)

  output <- new("RBfit", obj = obj, opt = opt, SD = SD, report = obj$report(obj$env$last.par.best), RBdata = RBdata)
  return(output)
}

#' Run MCMC
#'
#' Set up a stan model and runs the MCMC from a fitted TMB model (see example).
#'
#' @param RBfit An object of class \linkS4class{RBfit}.
#' @param priors_only If \code{TRUE}, turns off the likelihood and samples the priors.
#' @param chains The number of chains for running the MCMC.
#' @param iter The number of total iterations of the MCMC.
#' @param warmup The number of burn-in iterations of the MCMC.
#' @param thin How often the iterations are saved.
#' @param seed An integer for indexing random number generation.
#' @param cores The number of CPU cores for running MCMC in parallel.
#' @param lower A vector of the lower bounds of the parameters. This argument overrides the default, and thus generally not recommended to use.
#' @param upper A vector of the upper bounds of the parameters. This argument overrides the default, and thus generally not recommended to use.
#' @param ... More arguments to pass to \code{rstan::sampling} via \link[tmbstan]{tmbstan}, for example, \code{init} for starting values for the MCMC.
#' @seealso \link{fit_model} \link{summary.stanfit} \link{plot.stanfit} \link{report}
#' @examples
#' \donttest{
#' data(BC_lakes)
#' mod <- fit_model(BC_lakes[[42]]) # Run model
#'
#' pr_only <- run_mcmc(mod, priors_only = TRUE) # Run MCMC on priors
#' plot(pr_only) # Shows density plots of priors
#'
#' samps <- run_mcmc(mod) # Run MCMC with likelihood for data
#' plot(samps) # Shows density plots of posteriors
#' plot(samps, pr_only) # Shows density plots of both priors and posteriors
#'
#' report(samps) # Generate report
#' report(samps, pr_only)
#'
#' ##### Any plotting diagnostic function from the rstan package can be used.
#' summary(samps) # Summary of parameter values
#' stan_ac(samps) # Autocorrelation plot of MCMC chains
#' stan_trace(samps) # Trace plot of posteriors
#' stan_dens(samps) # Density plot of posteriors
#' }
#' @importFrom tmbstan tmbstan
#' @export
run_mcmc <- function(RBfit, priors_only = FALSE, chains = 2L, iter = 2e4, warmup = 0.5 * iter, thin = 5,
                     seed = 1, cores = chains, lower = NULL, upper = NULL, ...) {
  attachNamespace("tmbstan")
  attachNamespace("rstan")

  obj <- RBfit@obj
  if(priors_only) {
    obj_data <- obj$env$data
    obj_data$use_likelihood <- as.integer(FALSE)
    attr(obj_data, "check.passed") <- NULL

    obj_pars <- obj$env$parameters
    attr(obj_pars, "check.passed") <- NULL

    obj2 <- MakeADFun(obj_data, obj_pars, DLL = obj$env$DLL, hessian = obj$env$hessian, silent = obj$env$silent)
  } else obj2 <- obj

  tmbstan_args <- list(...)
  tmbstan_args <- c(tmbstan_args, list(obj = obj2, chains = chains, iter = iter, warmup = warmup,
                                       thin = thin, seed = seed, cores = cores))
  tmbstan_args$lower <- if(is.null(lower)) rep(0, 10) else lower
  tmbstan_args$upper <- if(is.null(upper)) rep(Inf, 10) else upper

  stan_obj <- do.call(tmbstan, tmbstan_args)
  stan_obj <- posterior_derived(stan_obj, obj2)

  assign("obj2", obj2, envir = stan_obj@.MISC)
  assign("RBfit", RBfit, envir = stan_obj@.MISC)
  assign("call", match.call(), envir = stan_obj@.MISC)

  return(stan_obj)
}

#' Generates stochastic length and age-length samples to run a simulation
#'
#' Simulate data for the assessment model.
#'
#' @param RBdata \linkS4class{RBdata} object used for setting up the simulation.
#' @param pars Matrix of parameters. See details.
#' @param N Numeric vector of length two for the number of age-length and length samples.
#' @param CV_L_stock Numeric, the CV of length at stocking among simulations.
#' @param CV_stock_d Numeric, the CV of stocking density among simulations.
#' @param seed Integer, for replicating random number generation.
#' @details The \code{pars} matrix should be nsim rows long and have columns.
#' @import stats
#' @import TMB
#' @export
simulation <- function(RBdata, pars, N = c(50, 100), CV_L_stock = 0, CV_stock_d = 0, seed = 1L) {
  RBdata <- validate_RBdata(RBdata)
  if(!is.matrix(pars)) {
    if(is.data.frame(pars) && all(apply(pars, 2, is.numeric))) {
      pars <- as.matrix(pars)
    } else stop("pars is not a data frame that can be converted to a matrix.", call. = FALSE)
  }
  if(ncol(pars) != 10) stop("Ten (10) columns needed for `pars` matrix. See help file.", call. = FALSE)
  nsim <- nrow(pars)

  Monte_Carlo_fn <- function(x, RBdata) { # x indexes simulation number
    set.seed(seed + x)
    n_age <- length(RBdata@Age)
    if(CV_L_stock > 0) {
      Lstock_mu <- RBdata@L_stock
      par_samp <- lognormal_par(list(Lstock_mu, Lstock_mu * CV_L_stock))
      RBdata@L_stock <- rlnorm(n_age, par_samp[1:n_age], par_samp[(n_age+1):(2*n_age)])
    }
    if(CV_stock_d > 0) {
      stock_mu <- RBdata@stock_density
      par_samp <- lognormal_par(list(stock_mu, stock_mu * CV_stock_d))
      RBdata@stock_density <- rlnorm(n_age, par_samp[1:n_age], par_samp[(n_age+1):(2*n_age)])
    }

    model_data <- list(Len_data = RBdata@Length, Len_age_data = RBdata@Age_length, age_adjust = RBdata@Age_adjust,
                       Lmid = RBdata@Length_bin, Lbin_width = unique(diff(RBdata@Length_bin)), stocking_density = RBdata@stock_density,
                       L_stock = RBdata@L_stock, bag_limit = RBdata@bag_limit, release_mortality = RBdata@release_mortality,
                       p_vrel = RBdata@p_vrel, use_likelihood = as.integer(FALSE), use_priors = as.integer(TRUE),
                       init_p_harvest = 1, nit_harvest = 5L)

    model_priors <- generate_priors(RBdata)

    model_par <- list(Linf = RBdata@prior_Linf[1], K = RBdata@prior_K[1], CV_Len = RBdata@prior_CV_Len[1], M = RBdata@prior_M[1],
                      q = RBdata@prior_q[1], Effort = RBdata@prior_Effort[1], GN_SL50 = RBdata@prior_GN_SL50[1], GN_gamma = RBdata@prior_GN_gamma[1],
                      angler_SL50 = RBdata@prior_angler_SL50[1], angler_gamma = RBdata@prior_angler_gamma[1])

    obj <- MakeADFun(data = c(model_data, model_priors), parameters = model_par, DLL = "RBassess", silent = TRUE)

    report <- obj$report(pars[x, ])
    pred_N <- report$survey_N
    pred_L <- report$survey_NL
    pred_F <- report$F

    samp_N <- rmultinom(1, N[1], pred_N)[, 1]
    samp_N <- matrix(samp_N, nrow = nrow(pred_N), ncol = ncol(pred_N))
    samp_L <- rmultinom(1, N[2], pred_L)[, 1]

    RBdata@Length <- samp_L
    RBdata@Age_length <- samp_N
    RBdata@Lake <- paste("Simulation No.", x)
    output <- list(RBdata = RBdata, sim = report)
    return(output)
  }

  res <- lapply(1:nsim, Monte_Carlo_fn, RBdata = RBdata)
  return(res)
}


