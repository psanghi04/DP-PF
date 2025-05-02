library(mcmc)

load(file = "data/vicData.RData")

## Choose variables
vicMat <- vicData %>% filter(Retir != 88888888, AGEGRP != 88)
vicMat <- vicMat %>% dplyr::select(Retir, AGEGRP)
vicMat$Retir <- ifelse(vicMat$Retir<99999999, 1, 0)
vicMat$AGEGRP <- ifelse(vicMat$AGEGRP<8, vicMat$AGEGRP*2.5+1, vicMat$AGEGRP*5-18)

vicMat = vicMat[2001:3000,] #### vicMat = vicMat[2901:3000,]

## Pre-process data for DP
zAGEGRP <- (vicMat$AGEGRP - min(vicMat$AGEGRP))/(max(vicMat$AGEGRP) - min(vicMat$AGEGRP)) 
xAGEGRP <- 2*zAGEGRP - 1# Your data: x and y, n = 1000
# Assuming x and y are already in your workspace

set.seed(49)
y = vicMat$Retir
x = xAGEGRP

# Define log-posterior function
log_posterior <- function(params, x, y) {
  beta0 <- params[1]
  beta1 <- params[2]
  a <- params[3]
  b <- params[4]
  
  #if (a <= 0 || b <= 0){return(-Inf)}
  
  eta <- beta0 + beta1 * x
  log_lik <- sum(y * eta - log(1 + exp(eta)))
  
  # Compute z = (x + 1)/2
  z <- (x + 1) / 2
  
  # Induced density from Beta on x
  epsilon <- 1e-6
  z_adj <- pmax(pmin(z, 1 - epsilon), epsilon)
  log_density_x <- sum(dbeta(z_adj, a, b, log = TRUE) - log(2))
  
  # Priors
  log_prior_beta0 <- dnorm(beta0, 0, 4, log = TRUE)
  log_prior_beta1 <- dnorm(beta1, 0, 4, log = TRUE)
  log_prior_a <- dgamma(a, 6, 4, log = TRUE)
  log_prior_b <- dgamma(b, 6, 4, log = TRUE)
  
  return(log_lik + log_density_x +
           log_prior_beta0 + log_prior_beta1 +
           log_prior_a + log_prior_b)
}

# Initial parameter vector
init_params <- c(beta0 = 0, beta1 = 0, a = 2, b = 2)

# Proposal standard deviations (you may tune these)
proposal_sd <- c(0.1, 0.1, 0.05, 0.05)

# Wrapper for use with metropolis() – returns *unnormalized* log posterior
posterior_wrapper <- function(params) {
  log_posterior(params, x, y)
}

# Run Metropolis sampler
set.seed(123)
samples <- metrop(
  obj = posterior_wrapper,
  initial = init_params,
  nbatch = 10000,
  blen = 1,
  scale = proposal_sd
)

samples$accept #0.2072

samples$batch[5000+(1:200)*25,1:2]

data <- data.frame(
  x = samples$batch[5000+(1:200)*25,1],
  y = samples$batch[5000+(1:200)*25,2],
  weight = rep(1/200,200)
)

# Extract beta0, beta1, and weights
beta <- as.matrix(data[, 1:2])
weights <- data$weight

mean_beta <- colSums(beta * weights) / sum(weights)

centered <- sweep(beta, 2, mean_beta)  # Center data
cov_beta <- t(centered) %*% (centered * weights) / sum(weights)

data$mahalanobis_dist <- mahalanobis(beta, center = mean_beta, cov = cov_beta)

data_sorted <- data %>% arrange(mahalanobis_dist)
data_sorted <- data_sorted %>%
  mutate(cum_weight = cumsum(weight) / sum(weight))
data_95 <- data_sorted %>% filter(cum_weight <= 0.95)

sum(data_95$weight)  # Should be close to 0.95

# Generate fitted logistic regression values and Plot the fitted logistic regression curves
x_seq <- seq(-1, 1, length.out = 100)
fitted_values <- data.frame(
  expand.grid(x = x_seq, id = 1:nrow(data_95)) # Grid for each beta pair
) %>%
  mutate(
    beta0 = data_95[id,1],
    beta1 = data_95[id,2],
    weight = data_95[id,3],
    logistic_pred = 1 / (1 + exp(-(beta0 + beta1 * x)))
  )

# Compute mean Retir probability/proportion for each AGEGRP
prob_data <- vicEmpPlot %>%
  group_by(AGEGRP) %>%
  summarise(mean_retir = mean(Retir))

# Create stepwise data for piecewise constant function
step_data <- prob_data %>%
  mutate(x_end = lead(AGEGRP, default = max(x_seq))) %>% # Extend to the next point
  rename(x_start = AGEGRP)

# Plot the logistic regression curves with the step function
# Create stepwise data for piecewise constant function
step_data <- prob_data %>%
  mutate(x_end = lead(AGEGRP, default = max(x_seq))) %>%
  rename(x_start = AGEGRP)

true_values <- data.frame(
  x = x_seq,
  logistic_pred = 1 / (1 + exp(-(-3.825 + 6.822 * x_seq))) # Fitted value again
)

# Plot with stepwise mean Retir (because of its interval data) overlay
step_plot <- ggplot() +
  geom_line(data = fitted_values, 
            aes(x = x, y = logistic_pred, group = id, alpha = weight),
            color = "blue") +
  geom_segment(data = step_data, 
               aes(x = x_start, xend = x_end, 
                   y = mean_retir, yend = mean_retir),
               color = "black", size = 1.2) +
  geom_line(data = true_values, 
            aes(x = x, y = logistic_pred), 
            color = "red", size = 1.2, linetype = "dashed") +
  scale_alpha_continuous(range = c(0.1, 0.8)) +
  labs(title = "95% Non-Private Posterior Curves",
       x = "x", 
       y = "Probability") +
  theme_minimal() + 
  geom_hline(yintercept=0.5, linetype = "dotted", linewidth = 2.5) +
  theme(
    axis.text.x = element_text(size = 17),  
    axis.text.y = element_text(size = 17),
    axis.title.x = element_text(size = 19),
    axis.title.y = element_text(size = 19),
    legend.text = element_text(size = 17),
    legend.title = element_text(size = 17),
    plot.title = element_text(size = 23),  
  )

ggsave("plots/0.95_non-private_posterior_curves", plot = step_plot, width = 8, height = 6)
