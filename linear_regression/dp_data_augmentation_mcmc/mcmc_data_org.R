library(sns)
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(writexl)

sim <- "conj"
save_loc <- ""

if (sim == "conj") {
  setwd("linear_regression/mcmc/conjugate_prior/data")
  save_loc <- "linear_regression/mcmc/conjugate_prior/results/"
} else {
  setwd("linear_regression/mcmc/non_conjugate_prior/data")
  save_loc <- "linear_regression/mcmc/non_conjugate_prior/results/"
}

sum_stat <- function(name) {
  file <- paste(name, ".RData", sep = "")
  load(file)

  # Extract the method
  method <- sub("^(.*?)_n\\d+e.*", "\\1", name)

  # Extract sample_size (n)
  sample_size <- as.numeric(sub(".*_n(\\d+)e.*", "\\1", name))

  # Extract epsilon (e)
  epsilon <- as.numeric(sub(".*e(\\d*\\.?\\d+).*", "\\1", name))

  # Extract seed
  seed <- as.numeric(sub(".*seed(\\d+).*", "\\1", name))

  # Extract particles (p)
  iter <- as.numeric(sub(".*i(\\d+).*", "\\1", name))

  chain <- get("chain")

  ## Post-simulation: We have 10,000 particles but treat the first 5,000 as burn-in (Save all for future reference)
  mu_chain = mean(chain$mu_chain[5001:10000])
  tau_chain = mean(chain$tau_chain[5001:10000])
  phi_chain = mean(chain$phi_chain[5001:10000]) #can check other parameters' behavior

  betaMat = as.matrix(chain$beta_chain) #beta_0, beta_1 are our focus
  ess_DAMCMC_beta = ess(t(betaMat[,5001:10000])) #ESS of beta_0 and beta_1 for DP-DA-MCMC from library(sns): Record them as usual
  mean_estimate = rowMeans(betaMat[,5001:10000]) #posterior mean estimate: each particle has equal weights
  beta_0 = betaMat[1,5001:10000]
  beta_1 = betaMat[2,5001:10000]
  beta_0_CI = quantile(beta_0, c(0.025, 0.975))
  beta_0_CI = beta_0_CI[2] - beta_0_CI[1]
  beta_1_CI = quantile(beta_1, c(0.025, 0.975)) #95% credible interval
  beta_1_CI = beta_1_CI[2] - beta_1_CI[1]

  time <- chain$time

  # Summarize the results
  ss <- c(mu_chain, tau_chain, phi_chain, ess_DAMCMC_beta, mean_estimate, beta_0_CI, beta_1_CI, time)
  c(method, sample_size, epsilon, iter, ss)
}

rdata_files <- list.files(pattern = "\\.RData$")
rdata_names <- sub("\\.RData$", "", rdata_files)

df <- sapply(rdata_names, function(name){
  sum_stat(name)
})

df <- as.data.frame(t(df))
df[, 2:14] <- lapply(df[, 2:14], as.numeric)
colnames(df) <- c("method", "sample_size", "epsilon", "iterations", "mu_chain",
    "tau_chain", "phi_chain", "ess_DAMCMC_beta_0", "ess_DAMCMC_beta_1", "mean_estimate_0", "mean_estimate_1", "CI_0", "CI_1", "time")

# Write the results to a file
write_xlsx(df, paste(save_loc, sim, "_sim_summary.xlsx", sep = ""))
