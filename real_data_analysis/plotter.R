library(ggplot2)
library(dplyr)
library(MASS)

# Loading individual chains
load("data/res_explore_n1000e0.5q0.5a6r4mo0so16mt0st16.RData")
beta0 <- res[[11]][1, ]
beta1 <- res[[11]][2, ]
weights <- res[[11]][5, ]

cat("beta0", sum(beta0 * weights), "\n")
cat("beta1", sum(beta1 * weights), "\n")
cat("ess", 1 / sum(weights^2), "\n")

# # Load all chains
# chain_dir <- "."
# chain_files <- list.files(chain_dir, pattern = "res_explore_n1000.*\\.RData", full.names = TRUE)

# # Initialize vectors to store combined data
# beta0 <- c()
# beta1 <- c()
# weights <- c()

# chain_beta0 <- c()
# chain_beta1 <- c()
# chain_ess <- c()

# # Iterate through all chain files
# for (chain_file in chain_files) {
#   # Load the chain
#   load(chain_file)
  
#   beta0 <- c(beta0, res[[11]][1,]/3)
#   beta1 <- c(beta1, res[[11]][2,]/3)
#   weights <- c(weights, res[[11]][5,])

#   print(chain_file)
#   print(sum(res[[11]][1,]*res[[11]][5,]))
#   print(sum(res[[11]][2,]*res[[11]][5,]))
#   print(1/sum(res[[11]][5,]^2))

#   chain_beta0 <- c(chain_beta0, sum(res[[11]][1,]*res[[11]][5,]))
#   chain_beta1 <- c(chain_beta1, sum(res[[11]][2,]*res[[11]][5,]))
#   chain_ess <- c(chain_ess, 1/sum(res[[11]][5,]^2))
# }

# # Print the final results
# cat("beta0:", mean(chain_beta0), "\n")
# cat("beta1:", mean(chain_beta1), "\n")
# cat("ess:", sum(chain_ess), "\n")

# weights <- weights / sum(weights)  # Normalize weights

# Fitted values
data <- data.frame(
  x = beta0,
  y = beta1,
  weight = weights
)

# Plot for beta1
beta_1_plot <- ggplot(data, aes(x = y)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 30, fill = "lightblue", color = "blue", alpha = 0.5) +  # Histogram
  geom_density(aes(weight = weight), fill = "red", alpha = 0.3) +  # Density plot (red)
  geom_vline(xintercept = 6.822, color = "red", linetype = "dashed", linewidth = 1) +  # True beta1 (red line)
  scale_y_continuous(name = "weight") +
  labs(title = expression(paste("Posterior Distribution of ", beta[1])) , x = expression(beta[1]))


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

cat("data_95 =", sum(data_95$weight), "\n")  # Should be close to 0.95


# Point Cloud
point_cloud <- ggplot(data, aes(x = x, y = y)) +
  geom_point(aes(size = weight, alpha = weight), color = "blue") +
  annotate("point", x = -3.825, y = 6.822, 
           color = "red", size = 4, shape = 8) +  # Non-private fitted value: a single point
  scale_size_continuous(range = c(1, 5)) +
  scale_alpha_continuous(range = c(0.1, 0.8)) +
  labs(title = expression(paste("Joint Posterior Distribution of ", beta[0], " and ", beta[1])),
       x = expression(beta[0]), y = expression(beta[1]))

# Generate fitted logistic regression values
# Plot the fitted logistic regression curves
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

curves <- ggplot(fitted_values, aes(x = x, y = logistic_pred, group = id, alpha = weight)) +
  geom_line(color = "blue") +
  scale_alpha_continuous(range = c(0.1, 0.8)) +
  labs(title = "Fitted Private and Non-private Curves",
       x = "x", 
       y = "Probability") +
  theme_minimal()

# Load the orginal data
load(file = "data/vicData.RData")


## Choosing features
vicMat <- vicData %>% filter(Retir != 88888888, AGEGRP != 88)
vicMat <- vicMat %>% dplyr::select(Retir, AGEGRP)
vicMat$Retir <- ifelse(vicMat$Retir<99999999, 1, 0)
vicMat$AGEGRP <- ifelse(vicMat$AGEGRP<8, vicMat$AGEGRP*2.5+1, vicMat$AGEGRP*5-18)

vicMat = vicMat[2001:3000,]

## Pre-process data for DP
zAGEGRP <- (vicMat$AGEGRP - min(vicMat$AGEGRP))/(max(vicMat$AGEGRP) - min(vicMat$AGEGRP)) 
xAGEGRP <- 2*zAGEGRP - 1
vicEmpPlot <- as.data.frame(cbind("Retir"=vicMat$Retir, "AGEGRP"=xAGEGRP))

glmfit <- glm(vicMat$Retir ~ xAGEGRP, family = binomial) 
print(summary(glmfit))


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
  logistic_pred = 1 / (1 + exp(-(-3.8 + 6.8 * x_seq))) # Fitted value again
)

# Plot with stepwise mean Retir (because of its interval data) overlay
step_plot <- ggplot() +
  geom_line(data = fitted_values, 
            aes(x = x, y = logistic_pred, group = id, alpha = weight),
            color = "blue") +
  geom_segment(data = step_data, 
               aes(x = x_start, xend = x_end, 
                   y = mean_retir, yend = mean_retir),
               color = "black", linewidth = 1.2) +
  geom_line(data = true_values, 
            aes(x = x, y = logistic_pred), 
            color = "red", linewidth = 1.2, linetype = "dashed") +
  scale_alpha_continuous(range = c(0.1, 0.8)) +
  labs(title = "95% Private Curves and Non-private Curve",
       x = "x",
       y = "Probability")


ggplot_show <- function() {
  print(beta_1_plot)
  print(point_cloud)
  print(step_plot)
}
ggplot_show()

# Save the plots
ggsave("plots/beta_1_plot.png", plot = beta_1_plot, width = 8, height = 6)
ggsave("plots/point_cloud.png", plot = point_cloud, width = 8, height = 6)
ggsave("plots/step_plot.png", plot = step_plot, width = 8, height = 6)