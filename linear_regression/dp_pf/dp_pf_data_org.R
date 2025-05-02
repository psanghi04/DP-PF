
library(ggplot2)
library(reshape2)
library(gridExtra)
library(gtable)
library(writexl)

sim <- "conj" # "conj" or "nonconj"
sdp_type <- "fixed" # "rand" or "fixed"
save_loc <- ""

if (sim == "conj") {
  setwd("linear_regression/dp_pf/conjugate_prior/data")
  save_loc <- "linear_regression/dp_pf/conjugate_prior/summary/"
} else {
  setwd("linear_regression/dp_pf/non_conjugate_prior/data")
  save_loc <- "linear_regression/dp_pf/non_conjugate_prior/summary/"
}

switch(sim,
  "conj" = {
    # 5 Predictors simulations
    true_mu = 1
    true_Phi = 1
    true_tau = 1
    true_beta = c(0, 2)
  },
  "nonconj" = {
    # 5 Predictors simulations
    true_mu = 1
    true_Phi = 1
    true_tau = 1
    true_beta =  c(0, 2)
  }
)

if (sdp_type == "fixed") {
  true_beta = c(0, 0)
  true_tau = 1
  true_mu = 0
  true_Phi = 1
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
  particles <- as.numeric(sub(".*p(\\d+).*", "\\1", name))

  lr <- get("res")

  #effective sample size
  sim_weight = sapply(lr[[11]], function(x) x$weightOld)
  eff_sample_size <- 1 / sum(sim_weight^2)

  sim_mu = sapply(lr[[11]], function(x) x$mu)
  sim_Phi = sapply(lr[[11]], function(x) x$Phi)
  sim_tau = sapply(lr[[11]], function(x) x$tau)
  sim_beta = sapply(lr[[11]], function(x) x$beta)

  # Function to calculate weighted quantiles
  weighted_quantile <- function(x, w, probs) {
    sorted_indices <- order(x)
    x_sorted <- x[sorted_indices]
    w_sorted <- w[sorted_indices]
    cum_weights <- cumsum(w_sorted) / sum(w_sorted)
    return(approx(cum_weights, x_sorted, xout = probs)$y)
  }

  switch(sim,
    "conj" = {
      # Calculate weighted quantiles for mu
      mu_quantile_low = weighted_quantile(sim_mu, sim_weight, 0.025)
      mu_quantile_high = weighted_quantile(sim_mu, sim_weight, 0.975)

      # Calculate weighted quantiles for beta
      beta_quantile_low = apply(sim_beta, 1, function(x, w) weighted_quantile(x, w, 0.025), w = sim_weight)
      beta_quantile_high = apply(sim_beta, 1, function(x, w) weighted_quantile(x, w, 0.975), w = sim_weight)
      
      # Calculate weighted quantiles for tau
      tau_quantile_low = weighted_quantile(sim_tau, sim_weight, 0.025)
      tau_quantile_high = weighted_quantile(sim_tau, sim_weight, 0.975)
    },
    "nonconj" = {
      # Calculate weighted quantiles for mu
      mu_quantile_low = weighted_quantile(sim_mu, sim_weight, 0.025)
      mu_quantile_high = weighted_quantile(sim_mu, sim_weight, 0.975)

      # Calculate weighted quantiles for beta
      beta_quantile_low = apply(sim_beta, 1, function(x, w) weighted_quantile(x, w, 0.025), w = sim_weight)
      beta_quantile_high = apply(sim_beta, 1, function(x, w) weighted_quantile(x, w, 0.975), w = sim_weight)
      
      # Calculate weighted quantiles for tau
      tau_quantile_low = weighted_quantile(sim_tau, sim_weight, 0.025)
      tau_quantile_high = weighted_quantile(sim_tau, sim_weight, 0.975) 
    }
  )

  switch(sdp_type, 
    "rand" = {
      # Calculate coverage for mu
      coverage_mu = mean(ifelse(mu_quantile_low < true_mu & true_mu < mu_quantile_high, 1, 0))

      # Calculate CI length for mu
      width_mu = mean(mu_quantile_high - mu_quantile_low)

      # Calculate coverage for beta
      coverage_beta = ifelse(beta_quantile_low < true_beta & true_beta < beta_quantile_high, 1, 0)

      # Calculate CI length for beta
      width_beta = beta_quantile_high - beta_quantile_low

      # Calculate coverage for tau
      coverage_tau = mean(ifelse(tau_quantile_low < true_tau & true_tau < tau_quantile_high, 1, 0))

      # Calculate CI length for tau
      width_tau = mean(tau_quantile_high - tau_quantile_low)

      # Calculate time
      time <- lr$time

      # Calculate second per ESS
      seconds_per_ess <- time / eff_sample_size

      sdp <- lr$sdp

      posterior_mean_beta = rowSums(sim_weight * sim_beta)
      posterior_variance_beta = rowSums(sim_weight * (sim_beta - posterior_mean_beta)^2)

      posterior_mean_tau <- sum(sim_weight * sim_tau)
      posterior_variance_tau <- sum(sim_weight * (sim_tau - posterior_mean_tau)^2)
      
      ss = c(eff_sample_size, posterior_mean_beta, posterior_variance_beta, posterior_mean_tau, posterior_variance_tau,
        coverage_mu, width_mu, coverage_tau, width_tau, coverage_beta, width_beta, time, seconds_per_ess, sdp)
    },

    "fixed" = {
      mu = rowSums(sim_weight * sim_beta)
      v = rowSums(sim_weight * (sim_beta - mu)^2)
      se = sqrt(v/eff_sample_size)
      alpha = 0.05
      z = qnorm(1-alpha/2)
      quantile_low = mu - z*se
      quantile_high = mu + z*se
      coverage_beta = ifelse(quantile_low < true_beta & true_beta < quantile_high, 1, 0)
      width_beta = quantile_high - quantile_low

      # Calculate time
      time <- lr$time

      # Calculate ESS per second
      seconds_per_ess <- time / eff_sample_size

      sdp <- lr$sdp

      ss = c(eff_sample_size, mu, v, width_beta, coverage_beta, time, seconds_per_ess, sdp)
    }
  )
  c(method, sample_size, epsilon, particles, ss)
}

rdata_files <- list.files(pattern = "\\.RData$")
rdata_names <- sub("\\.RData$", "", rdata_files)

df <- sapply(rdata_names, function(name){
  sum_stat(name)
})

df <- as.data.frame(t(df))

switch(sdp_type,
  "rand" = {
    switch(sim,
      "conj" = {
        df[, 2:26] <- lapply(df[, 2:26], as.numeric)
        colnames(df) <- c("method", "sample_size", "epsilon", "particles", "ESS", "post_mu_beta_0", "post_mu_beta_1", "post_var_beta_0", "post_var_beta_1",
                   "post_mean_tau", "post_var_tau", "mu_coverage", "mu_CI", "tau_coverage", "tau_CI", "beta_coverage_0", "beta_coverage_1", "beta_CI_0",
                   "beta_CI_1", "time", "seconds_per_ess", "sdp_0", "sdp_1", "sdp_2", "sdp_3", "sdp_4")
      },
      "nonconj" = {
        df[, 2:26] <- lapply(df[, 2:26], as.numeric)
        colnames(df) <- c("method", "sample_size", "epsilon", "particles", "ESS", "post_mu_beta_0", "post_mu_beta_1", "post_var_beta_0", "post_var_beta_1",
                   "post_mean_tau", "post_var_tau", "mu_coverage", "mu_CI", "tau_coverage", "tau_CI", "beta_coverage_0", "beta_coverage_1", "beta_CI_0",
                   "beta_CI_1", "time", "seconds_per_ess", "sdp_0", "sdp_1", "sdp_2", "sdp_3", "sdp_4")
      }
    )
  },
  "fixed" = {
    # Regular
    df[, 2:16] <- lapply(df[, 2:16], as.numeric)
    colnames(df) <- c("method", "sample_size", "epsilon", "particles", "ESS", "post_mu_0", "post_mu_1", "post_var_0", "post_var_1",
                   "beta_CI_0", "beta_CI_1", "beta_coverage_0", "beta_coverage_1", "time", "seconds_per_ess", "sdp_0", "sdp_1", "sdp_2", "sdp_3", "sdp_4")
  }
)

# Save the data
write_xlsx(df, paste(save_loc, sim, "_", sdp_type, "_sdp_sim_summary.xlsx", sep = ""))