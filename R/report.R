

#' @name report
#' @title HTML report for assessment model
#'
#' @description Returns an HTML document with figures showing data, model fits, and diagnostics.
#'
#' @param x An object of class `RBfit`, returned from \link{fit_model}.
#' @param dir The directory for saving the HTML document.
#' @param author A character string for the author.
#' @importFrom rmarkdown render
#' @export
setMethod("report", signature(x = "RBfit", y = "missing"),
          function(x, dir = tempdir(), author = NULL, bubble = 7, ...) {
            #old_dir <- getwd()
            #on.exit(setwd(old_dir))

            old_par <- par(no.readonly = TRUE)
            on.exit(par(old_par), add = TRUE)

            rep_file <- file.path(path.package("RBassess"), "report.Rmd")
            con <- file(rep_file, open = "r")
            rmd <- readLines(con)
            close(con)

            # Subtitle is either x or Lake name
            Lake_name <- ifelse(nchar(x@RBdata@Lake) > 0, x@RBdata@Lake, substitute(x))
            rmd[3] <- paste0("subtitle: \"", Lake_name, "\"")

            # Use author name if provided
            if(!is.null(author)) rmd[4] <- paste0("author: \"by ", author, "\"")

            # Generate data plots
            if(!dir.exists(dir)) {
              dir.create(dir)
              message(paste("Creating directory:", dir))
            }
            if(!dir.exists(file.path(dir, "figures"))) dir.create(file.path(dir, "figures"))

            message("Generating data plots...")

            png(file.path(dir, "figures", "data_age_length.png"), units = "in", res = 400,
                width = 5, height = 3.5)
            par(mar = c(5, 4, 1, 1))
            plot_age_length(x@RBdata, x, bubble = bubble)
            dev.off()

            png(file.path(dir, "figures", "data_age_length_residuals.png"), units = "in", res = 400,
                width = 5, height = 3.5)
            par(mar = c(5, 4, 1, 1))
            plot_age_length_residual(x, bubble = bubble)
            dev.off()

            png(file.path(dir, "figures", "data_length.png"), units = "in", res = 400,
                width = 5, height = 3.5)
            par(mar = c(5, 4, 1, 1))
            plot_length(x@RBdata, x)
            dev.off()

            png(file.path(dir, "figures", "data_length_residuals.png"), units = "in", res = 400,
                width = 5, height = 3.5)
            par(mar = c(5, 4, 1, 1))
            plot_length_residual(x)
            dev.off()

            png(file.path(dir, "figures", "data_stock_density.png"), units = "in", res = 400,
                width = 5, height = 3.5)
            par(mar = c(5, 4, 1, 1))
            plot_stocking_density(x@RBdata)
            dev.off()

            png(file.path(dir, "figures", "data_length_stock.png"), units = "in", res = 400,
                width = 5, height = 3.5)
            plot_Lstart(x@RBdata)
            dev.off()

            # Write Rmd to dir
            write(rmd, file = file.path(dir, "report.Rmd"))

            # Render file
            message("Rendering HTML file...")
            rmarkdown::render(file.path(dir, "report.Rmd"), "html_document", "report.html", dir)

            # Open html file
            browseURL(file.path(dir, "report.html"))
            invisible()
          })






