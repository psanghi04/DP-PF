
setwd("location_scale/rejection_abc/data")

library(writexl)

sum_stat <- function(name) {
  file <- paste(name, ".RData", sep = "")
  load(file)

  is_present = grepl("intermediate", name, ignore.case = TRUE)

  # Extract the method
  method <- sub("^(.*?)_n\\d+e.*", "\\1", name)

  # Extract sample_size (n)
  sample_size <- as.numeric(sub(".*_n(\\d+)e.*", "\\1", name))

  # Extract epsilon (e)
  epsilon <- as.numeric(sub(".*e(\\d*\\.?\\d+).*", "\\1", name))

  # Extract p
  particles <- as.numeric(sub(".*p(\\d+).*", "\\1", name))

  iter <- as.numeric(sub(".*i(\\d+).*", "\\1", name))
  
  rejection_ci <- get("rejection_abc_ci")

  eff_sample_size <- particles

  theta <- sapply(1:particles, function(j) {
      rejection_ci[[j]]$theta
  })

  post_mu <- mean(theta[1, ])
  variance_mu <- var(theta[1, ])
  quantile_high <- quantile(theta[1, ], 0.975)
  quantile_low <- quantile(theta[1, ], 0.025)
  width_mu <- quantile_high - quantile_low
  coverage_mu <- ifelse(quantile_low < 1 && 1 < quantile_high, 1, 0)

  post_sigma <- mean(theta[2, ])
  variance_sigma <- var(theta[2, ])
  quantile_high <- quantile(theta[2, ], 0.975)
  quantile_low <- quantile(theta[2, ], 0.025)
  width_sigma <- quantile_high - quantile_low
  coverage_sigma <- ifelse(quantile_low < 1 && 1 < quantile_high, 1, 0)

  time <- rejection_ci$time
  sdp <- rejection_ci$sdp
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

write_xlsx(df, "location_scale/rejection_abc/results/rejection_abc_final_sim_summary.xlsx")