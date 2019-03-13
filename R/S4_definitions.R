
#' Class-\code{RBdata}
#'
#' An S4 class for setting up inputs for the rainbow trout assessment.
#'
#' @name RBdata-class
#' @docType class
#'
#' @slot Lake Name of lake.
#' @slot Year Year of survey.
#' @slot Length_bin A vector of midpoints of the length bins. Width of all length bins assumed to be equal.
#' @slot Age A vector of integer ages. Should be consecutive from 1 to the maximum age.
#' @slot Age_adjust A vector of accumulated growing degree days (converted to years) corresponding to the integer ages.
#' @slot Length A vector (single year sample) for the length composition data. The vector should be of
#' length equal to \code{length(Length_bins)}.
#' @slot Age_length A matrix (single year sample) for the age-length sammples. The matrix should have
#' \code{length(Age)} rows and \code{length(Length_bins)} columns.
#' @slot L_stock A vector of lengths at stocking, in order of the corresponding ages.
#' @slot stock_density A vector of stocking density,  in order of the corresponding ages.
#' @slot bag_limit The daily bag limit (per day) used to calculated the probability of harvesting.
#' @slot M_release The assumed mortality from catch and release.
#' @slot prior_Linf Von Bertalanffy asymptotic length. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_K Von Bertalanffy growth coefficient. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_CV_Len Coefficient in variation in length-at-age. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_M Natural mortality. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_Effort Fishing effort. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_q Catchability coefficient for scaling effort. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_GN_SL50 Gillnet length at 50\% selectivity for the survey. Vector of length two for mean and standard deviation, respectively.
#' Default priors are based on values from Askey et al. (2007).
#' @slot prior_GN_gamma Steepness of selectivity for gillnet (positive values). Vector of length two for mean and standard deviation, respectively.
#' Default priors are based on values from Askey et al. (2007).
#' @slot prior_angler_gamma Angler length at 50\% selectivity. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_angler_SL50 Steepness of selectivity for angler (positive values). Vector of length two for mean and standard deviation, respectively.
#' @export
#' @examples
#' new_data <- new("RBdata")
#' @references
#' Asky, P.J. Post, J.R., Parkinson, E.A., Rivot, E., Paul, A.J., and Biro, P.A. 2007. Estimation of gillnet efficiency and selectivity across multiple
#' sampling units: A hierarchical Bayesian analysis using mark-recapture data. Fisheries Research 83: 162-174.
#' @seealso \link{fit_model} \link{run_mcmc} \link{BC_lakes} \link{plot,RBdata,ANY-method}
#' @exportClass RBdata
RBdata <- setClass("RBdata",
                   slots = c(Lake = "character", Year = "numeric", Length_bin = "vector", Age = "vector", Age_adjust = "vector",
                             Length = "vector", Age_length = "matrix", L_stock = "vector",
                             stock_density = "vector", bag_limit = "numeric", M_release = "numeric",
                             prior_Linf = "numeric", prior_K = "numeric", prior_CV_Len = "numeric", prior_M = "numeric",
                             prior_Effort = "numeric", prior_q = "numeric", prior_GN_SL50 = "numeric", prior_GN_gamma = "numeric",
                             prior_angler_SL50 = "numeric", prior_angler_gamma = "numeric"))

default_prior <- list(Linf = c(50, 10), K = c(0.2, 0.15), CV_Len = c(0.1, 0.1), M = c(0.5, 0.2), Effort = c(1, 1), q = c(0.05, 0.1),
                      GN_SL50 = c(10, 10), # Based on Askey et al (2007) and delta method
                      GN_gamma = c(0.753, 0.192), # Based on Askey et al (2007)
                      angler_SL50 = c(20, 10), #angler_SL50 = c(2, 2)
                      angler_gamma = c(0.33, 0.5))
names(default_prior) <- paste0("prior_", names(default_prior))

setMethod("initialize", "RBdata",
          function(.Object, silent = FALSE, ...) {
            dots <- list(...)
            if(length(dots) > 0) {
              for(i in 1:length(dots)) {
                if(names(dots)[[i]] %in% slotNames("RBdata")) {
                  slot(.Object, names(dots)[[i]]) <- dots[[i]]
                } else {
                  if(!silent) warning(paste0("'", names(dots)[[i]], "' is not a slot in class 'RBdata'"), call. = FALSE)
                }
              }
            }

            ## Use default priors if missing
            names_default_prior <- NULL
            for(i in 1:length(default_prior)) {
              if(length(slot(.Object, names(default_prior)[i])) == 0) {
                slot(.Object, names(default_prior)[i]) <- default_prior[[i]]
                names_default_prior <- c(names_default_prior, names(default_prior)[i])
              }
            }
            if(!silent && !is.null(names_default_prior)) message("Default values used for:\n", paste(names_default_prior, collapse = "\n"))
            .Object
          })

#' \code{plot} method for S4 class \code{RBdata}
#'
#' Plots length frequency and age-length data, as well as priors.
#'
#' @param x An object of class \linkS4class{RBdata}.
#' @param bubble A number for adjusting the bubble size in the plot for age-length data.
#' @param data Logical, whether a figure showing the data will be plotted.
#' @param prior Logical, whether a figure showing priors will be plotted.
#' @param ... Other miscellaneous variables (not currently used).
#' @export
setMethod("plot", "RBdata",
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


setOldClass("sdreport")

#' Class-\code{RBfit}
#'
#' An S4 class for posterior mode output from rainbow trout assessment model.
#'
#' @name RBfit-class
#' @docType class
#'
#' @slot obj The TMB assessment object.
#' @slot opt Output from \link[stats]{nlminb}.
#' @slot SD Output from \link[TMB]{sdreport}.
#' @slot report Output from TMB model.
#' @slot RBdata Data objected used to generate the assessment.
#' @seealso \link{plot,RBfit,missing-method} \link{summary,RBfit-method}
#' @exportClass RBfit
RBfit <- setClass("RBfit", slots = c(obj = "list", opt = "list", SD = "sdreport", report = "list", RBdata = "RBdata"))


#' \code{plot} method for S4 class \code{RBfit}
#'
#' Plots length frequency and age-length data, as well as priors.
#'
#' @param x An object of class \linkS4class{RBfit}.
#' @param bubble A number for adjusting the bubble size in the plot for age-length data.
#' @param diagnostic Logical, whether a figure of model diagnostics (residuals, model fits, etc.) will be plotted.
#' @param posterior Logical, whether a figure of posterior distributions will be plotted.
#' @param ... Other miscellaneous variables (not currently used).
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

#' \code{summary} method for S4 class \code{RBfit}
#'
#' Returns parameter estimates and standard deviations.
#'
#' @param object An object of class \linkS4class{RBfit}.
#' @export
setMethod("summary", "RBfit", function(object, ...) round(summary(object@SD), 2))




#' @rdname plot-stanfit-stanfit-method
#' @importClassesFrom rstan stanfit
setMethod("plot", signature(x = "stanfit", y = "missing"),
          function(x, ...) {
            Lake_name <- if(nchar(x@.MISC$RBfit@RBdata@Lake) > 0) x@.MISC$RBfit@RBdata@Lake else substitute(x)
            plot_pars(RBdata = x@.MISC$RBfit@RBdata, stan_obj = x, plot.title = paste("Posteriors for", Lake_name), plot_type = "MCMC")
            invisible()
          })




#' \code{plot} method for S4 class \code{stanfit}
#'
#' Returns a figure parameter priors and posteriors. See example from \link{run_mcmc}.
#'
#' @param x An object of class \linkS4class{stanfit}.
#' @param y Another object of class \linkS4class{stanfit} which only sampled the priors.
#' @seealso \link{run_mcmc}
#' @examples
#' \donttest{
#' data(BC_lakes)
#' mod <- fit_model(BC_lakes[[1]]) # Run model
#'
#' pr_only <- run_mcmc(mod, priors_only = TRUE) # Run MCMC on priors
#' plot(pr_only) # Shows density plots of priors
#'
#' samps <- run_mcmc(mod) # Run MCMC with likelihood for data
#' plot(samps) # Shows density plots of posteriors
#' plot(samps, pr_only) # Shows density plots of both priors and posteriors
#'
#' ##### Any plotting diagnostic function from the rstan package can be used.
#' summary(samps) # Summary of parameter values
#' rstan::stan_ac(samps) # Autocorrelation plot of MCMC chains
#' rstan::stan_trace(samps) # Trace plot of MCMC chains
#' rstan::stan_dens(samps) # Density plot of posteriors
#' }
#' @importClassesFrom rstan stanfit
setMethod("plot", signature(x = "stanfit", y = "stanfit"),
          function(x, y, ...) {
            Lake_name <- if(nchar(x@.MISC$RBfit@RBdata@Lake) > 0) x@.MISC$RBfit@RBdata@Lake else substitute(x)
            plot_pars(stan_obj = x, stan_prior = y, plot.title = paste("Priors and posteriors for", Lake_name), plot_type = "MCMC_both")
            invisible()
          })


setGeneric("report", function(x, y, ...) standardGeneric("report"))

