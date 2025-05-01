
setwd("location_scale/dp_pf")

library(writexl)

sum_stat <- function(name) {
  file <- paste(name, ".RData", sep = "")
  load(file)

  # Extract the method
  method <- sub("^(.*?)_n\\d+e.*", "\\1", name)

  # Extract sample_size (n)
  sample_size <- as.numeric(sub(".*_n(\\d+)e.*", "\\1", name))

  # Extract epsilon (e)
  epsilon <- as.numeric(sub(".*e(\\d*\\.?\\d+).*", "\\1", name))

  # Extract p
  particles <- as.numeric(sub(".*p(\\d+).*", "\\1", name))

  iter <- as.numeric(sub(".*i(\\d+).*", "\\1", name))
  
  abc_smc_ci <- get("smc_abc_ci")

  weight_old <- sapply(abc_smc_ci[[3]], function(x) x$weight_old)
  eff_sample_size <- 1 / sum(weight_old^2)

  theta_old <- sapply(abc_smc_ci[[3]], function(x) x$theta_old)
  theta_old_order <- order(theta_old[1, ], decreasing = FALSE)
  theta_old_sorted <- theta_old[, theta_old_order]

  cdf_mu <- cumsum(weight_old[theta_old_order])

  post_mu <- sum(weight_old * theta_old[1, ])
  variance_mu <- sum(weight_old * (theta_old[1, ] - post_mu)^2) / eff_sample_size

  # Interpolate the CDF values to find the quantiles corresponding to the probabilities
  quantile_low <- approx(cdf_mu, theta_old_sorted[1, ],
                          xout = 0.025)$y
  quantile_low <- ifelse(is.na(quantile_low),
                          0.025 * theta_old_sorted[1, 1] / theta_old_sorted[1, 3],
                          quantile_low)
  quantile_high <- approx(cdf_mu, theta_old_sorted[1, ], xout = 0.975)$y
  width_mu <- quantile_high - quantile_low
  coverage_mu <- ifelse(quantile_low < 1 && 1 < quantile_high, 1, 0)
  
  theta_old_order <- order(theta_old[2, ], decreasing = FALSE)
  theta_old_sorted <- theta_old[, theta_old_order]

  cdf_sigma <- cumsum(weight_old[theta_old_order])

  post_sigma <- sum(weight_old * theta_old[2, ])
  variance_sigma <- sum(weight_old * (theta_old[1, ] - post_sigma)^2) / eff_sample_size

  # Interpolate the CDF values to find the quantiles corresponding to the probabilities
  quantile_low <- approx(cdf_sigma, theta_old_sorted[2,], xout = 0.025)$y
  quantile_low <- ifelse(is.na(quantile_low), 0.025*theta_old_sorted[1,1]/theta_old_sorted[1,3], quantile_low)
  quantile_high <- approx(cdf_sigma, theta_old_sorted[2,], xout = 0.975)$y
  width_sigma <- quantile_high - quantile_low
  coverage_sigma <- ifelse(quantile_low < 1 && 1 < quantile_high, 1, 0)
  
  time <- abc_smc_ci$time
  sdp <- abc_smc_ci$sdp
  time_per_ess <- time / eff_sample_size
  c(method, iter, sample_size, epsilon, particles, eff_sample_size, width_mu, coverage_mu, width_sigma, coverage_sigma, time, sdp, post_mu, variance_mu, post_sigma, variance_sigma, time_per_ess)  
}

rdata_files <- list.files(pattern = "\\.RData$")
rdata_names <- sub("\\.RData$", "", rdata_files)

df <- sapply(rdata_names, function(name) {
  t(sum_stat(name))
}, simplify = FALSE)

df <- do.call(rbind, df)
df <- as.data.frame(df)

df[, 2:18] <- lapply(df[, 2:18], as.numeric)
str(df)
colnames(df) <- c("method", "iter", "sample_size", "epsilon", "particles", "ESS", "mu_CI", "mu_coverage", "sigma_CI", "sigma_coverage", "time", "sdp0", "sdp1",
                  "post_mean_mu", "post_var_mu", "post_mean_sigma", "post_var_sigma", "time_per_ess")

write_xlsx(df, "results/dp_pf_location_scale_sim_summary.xlsx")
