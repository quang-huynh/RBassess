
#' Class-\code{RBdata}
#'
#' An S4 class for setting up inputs for the rainbow trout assessment.
#'
#' @name RBdata-class
#' @docType class
#'
#' @slot Lake Name of lake.
#' @slot Year Year of survey.
#' @slot Length_bin A vector of midpoints of the length bins. Width of all length bins assumed to be equal. Typically in centimeters.
#' @slot Age A vector of integer ages. Should be consecutive from 1 to the maximum age.
#' @slot Age_adjust A vector of accumulated growing degree days (converted to years) corresponding to the integer ages.
#' @slot Length A vector (single year sample) for the length composition data. The vector should be of
#' length equal to \code{length(Length_bins)}.
#' @slot Age_length A matrix (single year sample) for the age-length sammples. The matrix should have
#' \code{length(Age)} rows and \code{length(Length_bins)} columns.
#' @slot L_stock A vector of mean lengths at stocking, in order of the corresponding ages.
#' @slot stock_density A vector of stocking density, in order of the corresponding ages. In units of numbers per hectacre.
#' @slot bag_limit The daily bag limit (number of fish per angler day) used to calculated the probability of harvesting.
#' @slot release_mortality Release mortality, i.e., the proportion of fish that die due to catch and release. Default is 0.1.
#' @slot p_vrel The proportion of fish caught by anglers that are subsequently voluntarily released. Default is zero.
#' @slot prior_Linf Von Bertalanffy asymptotic length. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_K Von Bertalanffy growth coefficient. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_CV_Len Coefficient in variation in length-at-age. Vector of length two for mean and standard deviation, respectively.
#' @slot prior_M Instantaneous natural mortality (per year). Vector of length two for mean and standard deviation, respectively.
#' @slot prior_Effort Fishing effort. Vector of length two for mean and standard deviation, respectively. In units of angler days per hectacre.
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
#' @seealso \link{fit_model} \link{run_mcmc} \link{BC_lakes} \link{plot.RBdata}
#' @import methods stats graphics grDevices utils
#' @exportClass RBdata
RBdata <- setClass("RBdata",
                   slots = c(Lake = "character", Year = "numeric", Length_bin = "vector", Age = "vector", Age_adjust = "vector",
                             Length = "vector", Age_length = "matrix", L_stock = "vector",
                             stock_density = "vector", bag_limit = "numeric", release_mortality = "numeric", p_vrel = "numeric",
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

            if(!silent) {
              if(!is.null(names_default_prior)) message("Default values used for:\n", paste(names_default_prior, collapse = "\n"))
              if(length(.Object@p_vrel) == 0) message("p_vrel set to 0")
              if(length(.Object@release_mortality) == 0) message("Release mortality set to 0.1")
            }

            if(length(.Object@p_vrel) == 0) .Object@p_vrel <- 0
            if(length(.Object@release_mortality) == 0) .Object@release_mortality <- 0.1
            .Object
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
#' @seealso \link{plot.RBfit} \link{summary.RBfit}
#' @exportClass RBfit
RBfit <- setClass("RBfit", slots = c(obj = "list", opt = "list", SD = "sdreport", report = "list", RBdata = "RBdata"))






