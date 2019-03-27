#' @name plot.RBdata
#' @aliases plot,RBdata,missing-method
#' @title \code{plot} methods for S4 classes in \code{RBassess} package
#'
#' @description A suite of plot functions for various objects in the RBassess packages. For data, plots length frequency and age-length data, as well as priors.
#' For posterior mode and MCMC model results, plots model fits and posteriors.
#'
#' @param x An object of class \linkS4class{RBdata}, \linkS4class{RBfit}, or \linkS4class{stanfit}.
#' @param y Optional, another object of class \linkS4class{stanfit} which only sampled the priors.
#' @param bubble A number for adjusting the bubble size in the plot for age-length data.
#' @param data Logical, whether a figure showing the data will be plotted.
#' @param prior Logical, whether a figure showing priors will be plotted.
#' @param diagnostic Logical, whether a figure of model diagnostics (residuals, model fits, etc.) will be plotted.
#' @param posterior Logical, whether a figure of posterior distributions will be plotted.
#' @param ... Other miscellaneous variables (not currently used).
#' @seealso \link{fit_model} \link{run_mcmc}
#' @examples
#' ###### View data
#' data(BC_lakes)
#' plot(BC_lakes[[1]])
#'
#' ###### Run model
#' mod <- fit_model(BC_lakes[[1]])
#' plot(mod)
#'
#' \donttest{
#' ###### Run MCMC
#' pr_only <- run_mcmc(mod, priors_only = TRUE) # Run MCMC on priors
#' plot(pr_only) # Shows density plots of priors
#'
#' samps <- run_mcmc(mod) # Run MCMC with likelihood for data
#'
#' plot(samps) # Shows density plots of posteriors
#' plot(samps, pr_only) # Shows density plots of both priors and posteriors
#'
#' ##### Any plotting diagnostic function from the rstan package are also available.
#' summary(samps) # Summary of parameter values
#' rstan::stan_ac(samps) # Autocorrelation plot of MCMC chains
#' rstan::stan_trace(samps) # Trace plot of MCMC chains
#' rstan::stan_dens(samps) # Density plot of posteriors
#' }
#' @export
setMethod("plot", signature(x = "RBdata", y = "missing"),
          function(x, bubble = 7, data = TRUE, prior = TRUE, ...) {

            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))

            if(!data && !prior) dat <- TRUE

            Lake_name <- if(nchar(x@Lake) > 0) x@Lake else substitute(x)

            ############ Plot data
            if(data) {
              par(mfrow = c(2, 2), mar = c(5, 4, 1, 1), oma = c(0, 0, 3, 0))
              plot_age_length(x, bubble = bubble)
              plot_length(x)
              plot_stocking_density(x)
              plot_Lstart(x)
              mtext(paste0("Data for ", Lake_name), outer = TRUE, side = 3, font = 2)
            }

            ############ Plot priors
            if(prior) plot_pars(x, plot.title = paste("Model priors for", Lake_name), plot_type = "analytical_prior")

            return(invisible())

          })


#' @rdname plot.RBdata
#' @aliases plot.RBfit
#' @export
setMethod("plot", signature(x = "RBfit", y = "missing"),
          function(x, bubble = 7, diagnostic = TRUE, posterior = TRUE, ...) {
            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par))

            if(!diagnostic && !posterior) diagnostic <- TRUE

            Lake_name <- if(nchar(x@RBdata@Lake) > 0) x@RBdata@Lake else substitute(x)

            ############ Plot diagnostics
            if(diagnostic) {
              par(mfrow = c(3, 2), mar = c(5, 5, 1, 1), oma = c(0, 0, 3, 0))
              plot_age_length(x@RBdata, x, bubble = bubble)
              plot_length(x@RBdata, x)
              plot_age_length_residual(x, bubble = bubble)
              plot_length_residual(x)
              plot_selectivity(x)
              plot_F_diagnostic(x)
              mtext(paste0("Model fit for ", Lake_name), outer = TRUE, side = 3, font = 2)
            }

            ############ Plot priors and posteriors
            if(posterior) plot_pars(RBdata = x@RBdata, RBfit = x, plot.title = paste("Parameter estimates for", Lake_name),
                                    plot_type = "posterior_mode")

            return(invisible())
          })


#' @rdname plot.RBdata
#' @aliases plot.stanfit
#' @importClassesFrom rstan stanfit
setMethod("plot", signature(x = "stanfit", y = "missing"),
          function(x, ...) {
            Lake_name <- if(nchar(x@.MISC$RBfit@RBdata@Lake) > 0) x@.MISC$RBfit@RBdata@Lake else substitute(x)
            plot_pars(RBdata = x@.MISC$RBfit@RBdata, stan_obj = x, plot.title = paste("Posteriors for", Lake_name), plot_type = "MCMC")
            invisible()
          })



#' @rdname plot.RBdata
#' @importClassesFrom rstan stanfit
setMethod("plot", signature(x = "stanfit", y = "stanfit"),
          function(x, y, ...) {
            Lake_name <- if(nchar(x@.MISC$RBfit@RBdata@Lake) > 0) x@.MISC$RBfit@RBdata@Lake else substitute(x)
            plot_pars(stan_obj = x, stan_prior = y, plot.title = paste("Priors and posteriors for", Lake_name), plot_type = "MCMC_both")
            invisible()
          })


#' @name summary
#' @title \code{summary} method for S4 class \code{RBfit}
#' @aliases summary,RBfit-method summary.RBfit
#'
#' @description Returns parameter estimates and standard deviations.
#'
#' @param object An object of class \linkS4class{RBfit} or \linkS4class{stanfit}.
#' @param description Logical, indicating whether an additional column with parameter descriptions.
#' @param digits The number of decimal places for rounding. Use \code{NA} for no rounding.
#' @param full Logical, returns more information for a stanfit object: posterior quantiles, Rhat, and effective sample size.
#' @param ... Additional arguments
#' @return A data frame.
#' @aliases summary.RBfit
#' @examples
#' data(BC_lakes)
#' res <- fit_model(BC_lakes[[69]])
#' summary(res)
#'
#' \donttest{
#' samps <- run_mcmc(res)
#' summary(samps)
#' }
#' @export
setMethod("summary", signature(object = "RBfit"),
          function(object, description = FALSE, digits = 2, ...) {
            output <- as.data.frame(summary(object@SD))
            output$`CV Estimate` <- output[, 2]/output[, 1]
            if(!is.na(digits)) output <- round(output, digits)
            if(description) output$Description <- add_description_to_summary(output)
            return(output)
          })

# #' @rdname summary
# #' @aliases summary.stanfit
# #' @importClassesFrom rstan stanfit
# #' @importMethodsFrom rstan summary
# #' @export
# setMethod("summary", signature(object = "stanfit"),
#           function(object, description = FALSE, digits = 2, full = FALSE, ...) {
#             output <- do.call(rstan_summary, list(object))$summary
#             output <- as.data.frame(output)
#
#             if(!full) output <- output[, c(1, 3)]
#             output$cv <- output$sd/output$mean
#             if(!is.na(digits)) output <- round(output, digits)
#             if(description) output$Description <- add_description_to_summary(output)
#             colnames(output) <- c("Mean", "Std. Dev.", "CV")
#             return(output)
#           })

summary_internal <- function(object, description = FALSE, digits = 2, full = FALSE, ...) {
  output <- do.call(rstan_summary, list(object))$summary
  output <- as.data.frame(output)
  if(!full) output <- output[, c(1, 3)]
  output$cv <- output$sd/output$mean
  if(!is.na(digits)) output <- round(output, digits)
  if(description) output$Description <- add_description_to_summary(output)
  colnames(output) <- c("Mean", "Std. Dev.", "CV")
  return(output)
}

rstan_summary <- getMethod("summary", "stanfit", where = asNamespace("rstan"))@.Data

par_to_match <- c("Linf", "K", "CV_Len", "M", "Effort", "q", "GN_SL50", "GN_gamma", "angler_SL50", "angler_gamma", "F", "p_harvest", "CPUE", "lp__",
                  "u", "F_retain", "F_release")
desc <- c("Von Bertalanffy asymptotic length", "Von Bertalanffy growth coefficient", "Coefficient of variation in length at age",
          "Natural mortality (per year)", "Fishing effort (angler days per hectacre)", "Catchability coefficient", "Length of 50% selectivity of gillnet",
          "Slope of gillnet selectivity ogive", "Length of 50% selectivity of angler", "Slope of angler selectivity ogive", "Fishing mortality (per year)",
          "Probability of harvest", "Catch-per-unit-effort (fish per angler day)", "Log posterior", "Harvest rate",
          "Fishing mortality associated with retention", "Fishing mortality from release")

add_description_to_summary <- function(dframe) {
  ind <- match(rownames(dframe), par_to_match)
  desc_out <- desc[ind]
  return(desc_out)
}
