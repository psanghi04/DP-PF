
setwd("linear_regression/dp_rejection_abc/data")
save_loc <- "linear_regression/dp_rejection_abc/summary/"

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

  rejection_ci <- get("rejection_lr_abc_ci")

  rejection_ci_result <- sapply(1:10, function(i) {
  
    eff_sample_size <- particles
    
    mu <- numeric(particles)
    tau <- numeric(particles)
    beta <- matrix(-99, ncol = 2, nrow = particles)
    
    for (j in 1:particles) {
      if (is_present && (length(rejection_ci[[i]]) < j || is.numeric(rejection_ci[[i]][[j]]))) {
        mu[j] <- -99
        tau[j] <- -99
        beta[j, ] <- c(-99, -99)
      } else {
        mu[j] <- rejection_ci[[i]][[j]][[1]]$mu
        tau[j] <- rejection_ci[[i]][[j]][[1]]$tau
        beta[j, ] <- rejection_ci[[i]][[j]][[1]]$beta
      }
    }

    mu <- mu[mu != -99]
    tau <- tau[tau != -99]
    beta <- beta[!(beta[, 1] == -99 & beta[, 2] == -99), ]
    
    # For identifying if the rejection abc did not finish
    if (length(mu) == 0) {
      return(c(method, i, sample_size, epsilon, particles, eff_sample_size, -1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    }
 
    post_mu <- mean(mu)
    variance_mu <- var(mu)
    quantile_high <- quantile(mu, 0.975)
    quantile_low <- quantile(mu, 0.025)
    width_mu <- quantile_high - quantile_low
    coverage_mu <- ifelse(quantile_low < 1 && 1 < quantile_high, 1, 0)

    post_tau <- mean(tau)
    variance_tau <- var(tau)
    quantile_high <- quantile(tau, 0.975)
    quantile_low <- quantile(tau, 0.025)
    width_tau <- quantile_high - quantile_low
    coverage_tau <- ifelse(quantile_low < 1 && 1 < quantile_high, 1, 0)
   
    post_beta_1 <- mean(beta[, 2])
    variance_beta_1 <- var(beta[, 2])
    quantile_high <- quantile(beta[, 2], 0.975)
    quantile_low <- quantile(beta[, 2], 0.025)
    width_beta_1 <- quantile_high - quantile_low
    coverage_beta_1 <- ifelse(quantile_low < 2 && 2 < quantile_high, 1, 0)
 
    time <- rejection_ci$time
    sdp <- rejection_ci$sdp

    if (time == -1 && length(mu) > 0) {
      time = 100 * 14400 / length(mu)
    }

    c(method, i, sample_size, epsilon, particles, eff_sample_size, time, post_mu, post_tau, post_beta_1, variance_mu, variance_tau, variance_beta_1, width_mu, width_tau, width_beta_1, coverage_mu, coverage_tau, coverage_beta_1, sdp) 
  })
}

rdata_files <- list.files(pattern = "\\.RData$")
rdata_names <- sub("\\.RData$", "", rdata_files)

df <- sapply(rdata_names, function(name) {
  t(sum_stat(name))
}, simplify = FALSE)

df <- do.call(rbind, df)
df <- as.data.frame(df)

df[, 2:24] <- lapply(df[, 2:24], as.numeric)
str(df)
colnames(df) <- c("method", "iter", "sample_size", "epsilon", "particles", "ESS", "time", "post_mu", "post_tau", "post_beta_1", "var_mu", "var_tau", "var_beta_1", "width_mu", "width_tau", "width_beta_1", "coverage_mu", "coverage_tau", "coverage_beta_1", "sdp0", "sdp1", "sdp2", "sdp3", "sdp4")
#df <- aggregate(. ~ epsilon + sample_size + particles, data = df[, -1], FUN = mean, na.rm = TRUE)

write_xlsx(df, paste(save_loc, "rejection_abc_sim_summary.xlsx", sep = ""))
