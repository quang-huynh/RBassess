
library(RBassess)
#data(BC_lakes)
#ind <- 69
#res <- fit_model(BC_lakes[[ind]])
#res_pr <- run_mcmc(res, priors_only = TRUE)
#res_mcmc_pathological <- run_mcmc(res)
#res_mcmc <- run_mcmc(res, seed = 10)
#save(res, res_pr, res_mcmc, res_mcmc_pathological, file = "vignettes/vignette_example.Rdata")
load('vignettes/vignette_example.Rdata')

png("vignettes/figures/data.png", res = 300, width = 6, height = 5, units = "in")
plot(BC_lakes[[69]], prior = FALSE)
dev.off()

png("vignettes/figures/prior.png", res = 300, width = 6, height = 7, units = "in")
plot(BC_lakes[[69]], data = FALSE)
dev.off()

png("vignettes/figures/fits.png", res = 300, width = 6, height = 7, units = "in")
plot(res, posterior = FALSE)
dev.off()

png("vignettes/figures/posterior_mode.png", res = 300, width = 6, height = 7, units = "in")
plot(res, diagnostic = FALSE)
dev.off()

png("vignettes/figures/posterior_sampled.png", res = 300, width = 6, height = 7, units = "in")
plot(res_mcmc, res_pr)
dev.off()

# MCMC summary table
write.csv(summary(res_mcmc)[[1]][-11, c(1, 3)], file = "vignettes/figures/mcmc_summary.csv")

pars <- c("Linf", "K", "CV_Len", "M", "Effort", "q", "GN_SL50", "GN_gamma", "angler_SL50", "angler_gamma",
          "F", "p_harvest")
rstan::stan_trace(res_mcmc, pars)
ggsave("vignettes/figures/trace.png", width = 10, height = 7)
