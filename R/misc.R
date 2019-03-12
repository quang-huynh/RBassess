
validate_RBdata <- function(RBdata) {
  n_lengthbin <- length(RBdata@Length_bin)
  n_age <- length(RBdata@Age)

  if(length(RBdata@Age_adjust) != n_age) stop("Lengths of RBdata@Age and RBdata@Age_adjust are not equal.", call. = FALSE)
  if(length(RBdata@Length) != n_lengthbin) {
    stop("Lengths of RBdata@Length is not equal to the number of length bins (RBdata@Length_bin).", call. = FALSE)
  }
  if(nrow(RBdata@Age_length) != n_age) {
    stop("The number of rows of RBdata@Age_length is not equal to the number of ages (RBdata@Age).", call. = FALSE)
  }
  if(ncol(RBdata@Age_length) != n_lengthbin) {
    stop("The number of columns of RBdata@Age_length is not equal to the number of length bins (RBdata@Length_bin).", call. = FALSE)
  }
  if(length(RBdata@L_stock) != n_age) stop("Length of RBdata@L_stock is not equal to the number of ages (RBdata@Age).", call. = FALSE)
  if(length(RBdata@stock_density) != n_age) stop("Length of RBdata@stock_density is not equal to the number of ages (RBdata@Age).", call. = FALSE)

  return(RBdata)
}

generate_priors <- function(RBdata) {
  res <- list()
  res$prior_Linf <- lognormal_par(RBdata@prior_Linf)
  res$prior_K <- lognormal_par(RBdata@prior_K)
  res$prior_CV_Len <- lognormal_par(RBdata@prior_CV_Len)
  res$prior_M <- lognormal_par(RBdata@prior_M)

  res$prior_q <- lognormal_par(RBdata@prior_q)
  res$prior_Effort <- lognormal_par(RBdata@prior_Effort)

  res$prior_GN_SL50 <- lognormal_par(RBdata@prior_GN_SL50)
  res$prior_GN_gamma <- lognormal_par(RBdata@prior_GN_gamma)
  res$prior_angler_SL50 <- lognormal_par(RBdata@prior_angler_SL50)
  res$prior_angler_gamma <- lognormal_par(RBdata@prior_angler_gamma)

  return(res)
}

lognormal_par <- function(x) {
  x_mean <- x[1]
  x_sd <- x[2]

  lnx_mean <- log(x_mean) - 0.5 * log(1 + (x_sd/x_mean)^2)
  lnx_sd <- sqrt(log(1 + (x_sd/x_mean)^2))

  return(c(lnx_mean, lnx_sd))
}

#beta_par <- function(x) {
#  x_mean <- x[1]
#  x_sd <- x[2]
#
#  beta_a <- x_mean * (x_mean * (1 - x_mean)/x_sd/x_sd - 1)
#  beta_b <- (1 - x_mean) * (x_mean * (1 - x_mean)/x_sd/x_sd - 1)
#
#  return(c(beta_a, beta_b))
#}





posterior_derived <- function(stan_obj, obj) {
  sim <- stan_obj@sim
  sim$n_flatnames <- 14
  sim$pars_oi <- c(sim$pars_oi, "F", "p_harvest", "CPUE")
  sim$fnames_oi <- c(sim$fnames_oi, "F", "p_harvest", "CPUE")
  sim$dims_oi <- c(sim$dims_oi, list(F = numeric(0), p_harvest = numeric(0), CPUE = numeric(0)))

  generate_derived_samples <- function(x) {

    res <- lapply(1:length(x$Linf),
                  function(i) {
                    par <- c(x$Linf[i], x$K[i], x$CV_Len[i], x$M[i], x$Effort[i], x$q[i], x$GN_SL50[i], x$GN_gamma[i],
                             x$angler_SL50[i], x$angler_gamma[i])
                    rep <- obj$report(par)
                    return(c(rep$F, rep$p_harvest, rep$CPUE))
                  })
    res <- do.call(rbind, res)
    x$F <- res[, 1]
    x$p_harvest <- res[, 2]
    x$CPUE <- res[, 3]

    inits <- attr(x, "inits")
    rep_init <- obj$report(inits)
    new_inits <- c(inits, rep_init$F, rep_init$p_harvest, rep_init$CPUE)
    attr(x, "inits") <- new_inits

    mean_pars <- attr(x, "mean_pars")
    new_mean_pars <- c(mean_pars, x$F, x$p_harvest, x$CPUE)
    attr(x, "mean_pars") <- new_mean_pars

    return(x)
  }

  new_samps <- lapply(sim$samples, generate_derived_samples)

  sim$samples <- new_samps
  stan_obj@sim <- sim

  stan_obj@model_pars <- c(stan_obj@model_pars, "F", "p_harvest", "CPUE")
  stan_obj@par_dims <- c(stan_obj@par_dims, list(F = numeric(0), p_harvest = numeric(0), CPUE = numeric(0)))

  add_inits <- function(x) {
    pars <- do.call(c, x)
    rep <- obj$report(pars)

    res <- c(x, list(F = rep$F, p_harvest = rep$p_harvest, CPUE = rep$CPUE))
    return(res)
  }
  stan_obj@inits <- lapply(stan_obj@inits, add_inits)
  stan_obj@.MISC <- new.env(parent = emptyenv())

  return(stan_obj)
}


