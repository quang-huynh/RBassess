
plot_age_length <- function(RBdata, RBfit = NULL, stan_obj = NULL, bubble = 7, sdlen = 0.95) {
  if(all(is.na(RBdata@Age_length))) {
    plot(NULL, NULL, typ = "n", xlab = "Age (accumulated degree years)", ylab = "Length", xlim = c(0, 1), ylim = c(0, 1))
    text(x = 0.5, y = 0.5, labels = "No age-length data.")
    return(invisible())
  }
  range_obs <- pretty(RBdata@Age_length)
  n1 <- range_obs[2]
  n2 <- pretty(quantile(RBdata@Age_length, na.rm = TRUE, probs = 0.9))[2]
  if(n2 < n1) n1 <- 0.5 * n2
  diameter_max <- bubble/n2

  plot(NULL, NULL, typ = "n", xlab = "Age (accumulated degree years)", ylab = "Length",
       xlim = range(RBdata@Age_adjust), ylim = range(RBdata@Length_bin))
  for(i in 1:length(RBdata@Age)) {
    for(j in 1:length(RBdata@Length_bin)) {
      if(!is.na(RBdata@Age_length[i, j]) && RBdata@Age_length[i, j] > 0) {
        points(RBdata@Age_adjust[i], RBdata@Length_bin[j], cex = 0.5 * diameter_max * pmin(RBdata@Age_length[i, j], n2), pch = 21,
               bg = "grey80")
      }
    }
  }
  if(!is.null(RBfit)) {
    order_vec <- order(RBdata@Age_adjust)
    lines(RBdata@Age_adjust[order_vec], RBfit@report$Len_age[order_vec], col = "red", lwd = 3)
    len_lims <- c(0.5 * (1 - sdlen), 1 - 0.5 * (1 - sdlen))
    lines(RBdata@Age_adjust[order_vec], qnorm(len_lims[1], RBfit@report$Len_age, RBfit@report$sd)[order_vec], col = "red", lty = 2)
    lines(RBdata@Age_adjust[order_vec], qnorm(len_lims[2], RBfit@report$Len_age, RBfit@report$sd)[order_vec], col = "red", lty = 2)
  }
  if(!is.null(stan_obj)) {
    order_vec <- order(RBdata@Age_adjust)

    sim <- stan_obj@sim
    ind <- 1:(sim$warmup %/% sim$thin)
    mcmc_samps <- lapply(sim$samples, function(x) lapply(x, `[`, -ind))

    get_length_at_age <- function(Linf, K, RBdata) {
      RBdata@L_stock * exp(-K * RBdata@Age_adjust) + Linf * (1 - exp(-K * RBdata@Age_adjust))
    }
    map_fn <- function(y, RBdata) {
      res <- Map(get_length_at_age, Linf = y$Linf, K = y$K, MoreArgs = list(RBdata = RBdata))
      do.call(rbind, res)
    }
    get_matrix <- function(x, RBdata) {
      res <- lapply(x, map_fn, RBdata = RBdata)
      do.call(rbind, res)
    }
    res <- get_matrix(mcmc_samps, RBdata)
    lines(RBdata@Age_adjust[order_vec], colMeans(res)[order_vec], col = "red", lwd = 3)
  }

  legend("topright", legend = paste("N =", sum(RBdata@Age_length, na.rm = TRUE)), bty = "n", text.font = 2)
  legend("bottomright", legend = c(n1, paste0(">", n2)), pt.cex = 0.5 * diameter_max * c(n1, n2), pt.bg = "grey80", pch = 21, horiz = TRUE)

  return(invisible())

}

plot_length <- function(RBdata, RBfit = NULL) {
  Length_bin <- RBdata@Length_bin
  Length <- RBdata@Length
  Length[is.na(Length) | Length <= 0] <- NA

  if(!is.null(RBfit)) {
    pred_N <- RBfit@report$survey_NL/sum(RBfit@report$survey_NL)
    if(!all(is.na(Length))) pred_N <- sum(Length, na.rm = TRUE) * pred_N
    ymax <- 1.1 * max(c(Length, pred_N), na.rm = TRUE)
    barplot(Length, names.arg = Length_bin, space = 0, xlab = "Length", ylab = "Numbers-at-length", ylim = c(0, ymax))
    lines(c(1:length(Length_bin)) - 0.5, pred_N, col = "red", lwd = 2, typ = "o", pch = 16)
    if(sum(Length, na.rm = TRUE) > 0) legend("topright", legend = paste("N =", sum(Length, na.rm = TRUE)), bty = "n", text.font = 2)
    #legend("topleft", legend = c("Observed", "Predicted"), col = c("black", "red"), pt.bg = c("grey80", "red"), pch = c(22, 16))

  } else {
    if(all(is.na(Length))) ymax <- 1 else ymax <- 1.1 * max(Length, na.rm = TRUE)
    barplot(Length, names.arg = Length_bin, space = 0, xlab = "Length", ylab = "Numbers-at-length", ylim = c(0, ymax))

    if(all(is.na(Length))) {
      text(x = 0.5 * length(Length), y = 0.5, labels = "No length data.")
    } else legend("topright", legend = paste("N =", sum(Length, na.rm = TRUE)), bty = "n", text.font = 2)
  }
  box()
  return(invisible())
}

plot_stocking_density <- function(RBdata) {
  plot(RBdata@Age, RBdata@stock_density, xlab = "Age", ylab = "Stocking density", typ = "o", pch = 16)
}

plot_Lstart <- function(RBdata) {
  plot(RBdata@Age, RBdata@L_stock, xlab = "Age", ylab = "Length at stocking", typ = "o", pch = 16)
}

plot_age_length_residual <- function(RBfit, bubble = 7, add_title = TRUE) {
  if(all(is.na(RBfit@RBdata@Age_length))) {
    plot(NULL, NULL, typ = "n", xlab = "Age (accumulated degree years)", ylab = "Length", xlim = c(0, 1), ylim = c(0, 1))
    text(x = 0.5, y = 0.5, labels = "No age-length data.")
    return(invisible())
  }

  obs <- RBfit@RBdata@Age_length
  obs[obs <= 0] <- NA
  pred_N <- sum(RBfit@RBdata@Age_length, na.rm = TRUE) * RBfit@report$survey_N/sum(RBfit@report$survey_N)
  resids <- (obs - pred_N)/sqrt(pred_N)

  range_resids <- quantile(pretty(resids))[c(2, 4)]
  diameter_max <- bubble/max(abs(range_resids))

  plot(NULL, NULL, typ = "n", xlab = "Age (accumulated degree years)", ylab = "Length",
       xlim = range(RBfit@RBdata@Age_adjust), ylim = range(RBfit@RBdata@Length_bin))
  for(i in 1:length(RBfit@RBdata@Age)) {
    for(j in 1:length(RBfit@RBdata@Length_bin)) {
      if(!is.na(obs[i, j])) {
        points(RBfit@RBdata@Age_adjust[i], RBfit@RBdata@Length_bin[j],
               cex = 0.5 * diameter_max * pmin(abs(resids[i, j]), max(abs(range_resids))), pch = 21,
               bg = ifelse(resids[i, j] > 0, "grey80", "white"))
      }
    }
  }
  if(add_title) title("Age-length residual")

  legend("bottomright", legend = paste(c("<", ">"), range_resids), pt.cex = 0.5 * diameter_max * abs(range_resids),
         pt.bg = ifelse(range_resids > 0, "grey80", "white"), pch = 21, horiz = TRUE)

  return(invisible())

}

plot_length_residual <- function(RBfit, ylab = "Standardized\nPearson Residual") {
  obs <- RBfit@RBdata@Length
  obs[obs <= 0] <- NA
  pred_N <- sum(obs, na.rm = TRUE) * RBfit@report$survey_NL/sum(RBfit@report$survey_NL)
  resids <- (obs - pred_N)/sqrt(pred_N)

  if(all(is.na(resids))) ylim <- c(0, 1) else {
    if(all(resids < 0 | is.na(resids))) {
      ylim <- c(1.1 * min(resids, na.rm = TRUE), 0.1)
    } else if(all(resids > 0 | is.na(resids))) {
      ylim <- c(-0.1, 1.1 * max(resids, na.rm = TRUE))
    } else ylim <- 1.1 * range(resids, na.rm = TRUE)
  }

  barplot(rep(0, length(obs)), names.arg = RBfit@RBdata@Length_bin, space = 0, border = NA,
          xlab = "Length", ylim = ylim)
  if(all(is.na(obs))) {
    text(x = 0.5 * length(obs), y = 0.5, labels = "No length data.")
  } else {
    lines(c(1:length(obs)) - 0.5, resids, lwd = 2, typ = "o", pch = 16)
    abline(h = 0, lty = 3)
  }
  box()

  return(invisible())
}

plot_selectivity <- function(RBfit, stan_obj = NULL) {
  Length_bin <- RBfit@RBdata@Length_bin

  barplot(rep(0, length(Length_bin)), names.arg = Length_bin, space = 0, border = NA,
          xlab = "Length", ylab = "Selectivity", ylim = c(-0.1, 1.1), yaxp = c(0, 1, 2))
  abline(h = 0, col = "grey")

  if(is.null(stan_obj)) {
    GN <- RBfit@report$sel_GN
    ang <- RBfit@report$sel_angler
  } else {

    sim <- stan_obj@sim
    ind <- 1:(sim$warmup %/% sim$thin)
    mcmc_samps <- lapply(sim$samples, function(x) lapply(x, `[`, -ind))

    get_sel <- function(L50, gamma, Len_bin) {
      sel <- 1/(1 + exp(-gamma * (Len_bin - L50))); return(sel/max(sel))
    }
    map_fn <- function(y, Len, L50c, gammac) {
      res <- Map(get_sel, L50 = getElement(y, L50c), gamma = getElement(y, gammac), MoreArgs = list(Len = Len))
      do.call(rbind, res)
    }
    get_matrix <- function(x, Len, L50c, gammac) {
      res <- lapply(x, map_fn, Len = Len, L50c = L50c, gammac = gammac)
      do.call(rbind, res)
    }
    GN_mat <- get_matrix(mcmc_samps, Len = Length_bin, L50c = "GN_SL50", gammac = "GN_gamma")
    GN <- colMeans(GN_mat)
    ang_mat <- get_matrix(mcmc_samps, Len = Length_bin, L50c = "angler_SL50", gammac = "angler_gamma")
    ang <- colMeans(ang_mat)
  }

  lines(1:length(Length_bin) - 0.5, GN, lwd = 3)
  lines(1:length(Length_bin) - 0.5, ang, lwd = 3, lty = 3)
  legend("bottomright", legend = c("Gillnet", "Angler"), lty = c(1, 3), lwd = 2)
  box()

  return(invisible())
}

plot_F_diagnostic <- function(RBfit) {
  plot(RBfit@report$F_vec, typ = "o", xlab = "Iteration", ylab = "F", pch = 16, lwd = 2)
}

plot_pars <- function(RBdata = NULL, RBfit = NULL, stan_obj = NULL, stan_prior = NULL, plot.title = "Parameters",
                      plot_type = c("analytical_prior", "posterior_mode", "MCMC", "MCMC_both")) {
  plot_type <- match.arg(plot_type)

  prior_names <- substring(names(default_prior), 7)

  if(plot_type == "posterior_mode") {
    SD <- summary(RBfit@SD)
    if(any(is.nan(SD))) return(invisible())
  }

  if(grepl("MCMC", plot_type)) {
    sim <- stan_obj@sim
    ind <- 1:(sim$warmup %/% sim$thin)
    mcmc_samps <- lapply(sim$samples, function(x) lapply(x, `[`, -ind))
  }

  if(plot_type == "MCMC_both") {
    sim2 <- stan_prior@sim
    ind <- 1:(sim2$warmup %/% sim2$thin)
    mcmc_prior_samps <- lapply(sim2$samples, function(x) lapply(x, `[`, -ind))
  }

  par(mfrow = c(4, 3), mar = c(4, 3, 1, 1), oma = c(ifelse(is.null(RBfit) & is.null(stan_obj), 0, 2), 3, 3, 0))
  for(i in 1:length(prior_names)) {

    # Priors
    x_min <- 0
    x_max <- 1e-3
    if(plot_type == "analytical_prior" || plot_type == "posterior_mode" || plot_type == "MCMC") {
      prior_par <- slot(RBdata, paste0("prior_", prior_names[i]))

      x_min <- max(0, prior_par[1] - 2 * prior_par[2])
      x_max <- prior_par[1] + 2 * prior_par[2]
    } else if(plot_type == "MCMC_both") {
      dens_prior <- lapply(mcmc_prior_samps, getElement, prior_names[i])
      dens_prior <- do.call(c, dens_prior)
      get_dens_prior <- density(dens_prior, from = min(dens_prior), to = max(dens_prior))

      prior_par <- c(mean(dens_prior), sd(dens_prior))

      x_min <- min(0, min(dens_prior))
      x_max <- max(dens_prior)
    }

    # Posteriors
    if(plot_type == "posterior_mode") {
      post_par <- SD[match(prior_names[i], rownames(SD)), ]
      x_min_post <- max(0, post_par[1] - 2 * post_par[2])
      x_max_post <- post_par[1] + 2 * post_par[2]

      x_min <- min(x_min, x_min_post)
      x_max <- max(x_max, x_max_post)
    }
    if(grepl("MCMC", plot_type)) {
      dens <- lapply(mcmc_samps, getElement, prior_names[i])
      dens <- do.call(c, dens)
      get_dens <- density(dens, from = min(dens), to = max(dens))
      post_par <- c(mean(dens), sd(dens))

      x_min <- min(x_min, min(dens))
      x_max <- max(x_max, max(dens))
    }

    # Find appropriate limits for x-axis
    if(x_min > x_max) x_min <- 0.2 * x_max
    if((plot_type == "posterior_mode" | plot_type == "MCMC_both") && x_max > 20 * max(prior_par[1], post_par[1])) {
      x_max <- 20 * max(prior_par[1], post_par[1])
    }

    # Get probability density for prior values of x-axis
    if(plot_type == "analytical_prior" | plot_type == "posterior_mode") {
      x_vec <- seq(x_min, x_max, length.out = 100)
      prior_par2 <- lognormal_par(prior_par)
      y_vec <- dlnorm(x_vec, prior_par2[1], prior_par2[2])
    }

    ####### Make plot
    plot(NULL, NULL, xlab = prior_names[i], ylab = "", typ = "n", yaxp = c(0, 1, 2), xlim = c(x_min, x_max), ylim = c(0, 1))
    abline(h = 0, col = "grey")

    # Plot prior
    if(!grepl("MCMC", plot_type)) {
      lines(x_vec, y_vec/max(y_vec), lwd = 2)
      abline(v = prior_par[1], lwd = 2, lty = 2)
    }
    if(plot_type == "MCMC_both") {
      lines(get_dens_prior$x, get_dens_prior$y/max(get_dens_prior$y), lwd = 2)
      abline(v = mean(dens_prior), lwd = 2, lty = 2)
    }

    # Plot posterior
    if(plot_type == "posterior_mode") {
      y_vec_post <- dnorm(x_vec, post_par[1], post_par[2])
      lines(x_vec, y_vec_post/max(y_vec_post), lwd = 2, col = "red")
      abline(v = post_par[1], lwd = 2, lty = 2, col = "red")
    } else if(grepl("MCMC", plot_type)) {
      lines(get_dens$x, get_dens$y/max(get_dens$y), lwd = 2, col = "red")
      abline(v = mean(dens), lwd = 2, lty = 2, col = "red")
    }

    # Plot legend
    if(plot_type == "analytical_prior") {
      legend("topright", legend = c(paste("Mean =", round(prior_par[1], 2)), paste("SD =", round(prior_par[2], 2))))
    }
  }

  ############ F and p_harvest plots
  if(plot_type != "analytical_prior") {

    par_vec <- c("F", "p_harvest")
    for(i in 1:length(par_vec)) {

      # Prior
      if(plot_type == "MCMC_both") {
        dens_prior <- lapply(mcmc_prior_samps, getElement, par_vec[i])
        dens_prior <- do.call(c, dens_prior)
        get_dens_prior <- density(dens_prior, from = min(dens_prior), to = max(dens_prior))

        prior_par <- c(mean(dens_prior), sd(dens_prior))

        x_min <- min(dens_prior)
        x_max <- max(dens_prior)
      }

      # posterior
      if(plot_type == "posterior_mode") {
        post_par <- SD[match(par_vec[i], rownames(SD)), ]
        x_min <- max(0.01, post_par[1] - 2 * post_par[2])
        x_max <- post_par[1] + 2 * post_par[2]

        if(x_min > x_max) x_min <- 0.2 * x_max
        x_vec <- seq(x_min, x_max, length.out = 100)
        y_vec_post <- dnorm(x_vec, post_par[1], post_par[2])
      } else {
        dens <- lapply(mcmc_samps, getElement, par_vec[i])
        dens <- do.call(c, dens)
        get_dens <- density(dens, from = min(dens), to = max(dens))
        post_par <- c(mean(dens), sd(dens))

        x_min <- min(dens)
        x_max <- max(dens)
        if(x_max > 20 * post_par[1]) x_max <- 20 * post_par[1]

        x_vec <- get_dens$x
        y_vec_post <- get_dens$y
      }

      plot(NULL, NULL, xlab = par_vec[i], ylab = "", typ = "n", yaxp = c(0, 1, 2),
           xlim = c(x_min, x_max), ylim = c(0, 1))
      abline(h = 0, col = "grey")

      # Plot prior
      if(plot_type == "MCMC_both") {
        lines(get_dens_prior$x, get_dens_prior$y/max(get_dens_prior$y), lwd = 2)
        abline(v = prior_par[1], lwd = 2, lty = 2)
      }

      # Plot posterior
      lines(x_vec, y_vec_post/max(y_vec_post), lwd = 2, col = "red")
      abline(v = post_par[1], lwd = 2, lty = 2, col = "red")

    }
  }

  mtext(plot.title, outer = TRUE, side = 3, font = 2)
  mtext("Relative density", outer = TRUE, side = 2)

  if(plot_type == "posterior_mode" || plot_type == "MCMC_both") {
    par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)

    plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
    legend("bottom", legend = c("Prior", "Posterior"), col = c("black", "red"), lwd = 2, bty = "n", horiz = TRUE)
  }

  invisible()
}
