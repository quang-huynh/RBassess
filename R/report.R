
#' @name report
#'
#' @title HTML report for assessment model
#'
#' @description Returns an HTML document with figures showing data, model fits, and diagnostics.
#'
#' @param x An object of class \code{RBfit} returned from \link{fit_model} or class \code{stanfit}
#' returned from \link{run_mcmc}.
#' @param y An optional object of class \code{stanfit} which only sampled the priors.
#' @param dir The directory for saving the HTML document.
#' @param author A character string for the author.
#' @param filename Optional, string for name of the HTML file.
#' @param bubble Numeric for scaling size of bubbles in bubble plot. See \link{plot.RBdata}.
#' @param ... Miscellaneous arguments to pass.
#' @export
setGeneric("report", function(x, y, ...) standardGeneric("report"))

#' @rdname report
#' @aliases report,RBfit-missing-method
#' @importFrom rmarkdown render
#' @export
setMethod("report", signature(x = "RBfit", y = "missing"),
          function(x, dir = tempdir(), author = NULL, filename = NULL, bubble = 7, ...) {
            report_internal(x = x, dir = dir, author = author, filename = filename, bubble = bubble, ...)
          })

#' @rdname report
#' @aliases report,stanfit-missing-method
#' @importFrom rmarkdown render
#' @export
setMethod("report", signature(x = "stanfit", y = "missing"),
          function(x, dir = tempdir(), author = NULL, filename = NULL, bubble = 7, ...) {
            report_internal(x = x, dir = dir, author = author, filename = filename, bubble = bubble, ...)
          })

#' @rdname report
#' @aliases report,stanfit-stanfit-method
#' @importFrom rmarkdown render
#' @export
setMethod("report", signature(x = "stanfit", y = "stanfit"),
          function(x, y, dir = tempdir(), author = NULL, filename = NULL, bubble = 7, ...) {
            report_internal(x = x, dir = dir, author = author, filename = filename, bubble = bubble, ...)
          })

report_internal <- function(x, y = NULL, dir, author, filename, bubble = 7, ...) {
  if(inherits(x, "RBfit")) {
    if(is.null(filename)) filename <- "report_RBfit"
    rep_file <- file.path(path.package("RBassess"), "report_RBfit.Rmd")
    Lake_name <- ifelse(nchar(x@RBdata@Lake) > 0, x@RBdata@Lake, substitute(x))
  }
  if(inherits(x, "stanfit")) {
    if(is.null(filename)) filename <- "report_stanfit"
    rep_file <- file.path(path.package("RBassess"), "report_stanfit.Rmd")
    Lake_name <- ifelse(nchar(x@.MISC$RBfit@RBdata@Lake) > 0, x@.MISC$RBfit@RBdata@Lake, substitute(x))
  }
  con <- file(rep_file, open = "r")
  rmd <- readLines(con)
  close(con)

  # Subtitle is either x or Lake name
  rmd[3] <- paste0("subtitle: \"", Lake_name, "\"")

  # Use author name if provided
  if(!is.null(author)) rmd[4] <- paste0("author: \"by ", author, "\"")

  # Generate data plots
  if(!dir.exists(dir)) {
    dir.create(dir)
    message(paste("Creating directory:", dir))
  }

  # Write Rmd to dir
  filename_html <- paste0(filename, ".html")
  filename_rmd <- paste0(filename, ".Rmd")
  write(rmd, file = file.path(dir, filename_rmd))

  # Render file
  message("Rendering HTML file...")
  rmarkdown::render(file.path(dir, filename_rmd), "html_document", filename_html, dir,
                    output_options = list(df_print = "paged"))

  # Open html file
  browseURL(file.path(dir, filename_html))
  invisible()
}


generate_summary_table <- function(x) {
  if(inherits(x, "RBfit")) res <- as.data.frame(summary(x))
  if(inherits(x, "stanfit")) res <- as.data.frame(summary_internal(x))

  pars <- rownames(res)
  ind <- match(pars, par_to_match)
  res$Description <- desc[ind]
  return(res[, c(4, 1:3)])
}


