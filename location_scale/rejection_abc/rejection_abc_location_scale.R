library(doParallel)
library(emdbook)
library(foreach)
library(invgamma)
library(LaplacesDemon)
library(MASS)
library(tictoc)
library(truncnorm)
library(VGAM)
library(jsonlite)
library(optparse)

# Path to the JSON configuration file
# config_file_path <- system.file("extdata", "dp_config.json",
#                                 package = "smcDiffPrivacy")

config_file_path <- "~/smcDiffPrivacy/inst/extdata/rejection_abc_config.json"

# Read the configuration file
config <- fromJSON(config_file_path)

# Accessing configuration values
seed <- config$seed
mu <- config$mu
sigma_square <- config$sigma_square

dp_info <- config$dp_info
dp_info$sample_size <- config$sample_size

hyper_params <- config$hyper_params
hyper_params$hp <- config$hyper_params$hp
hyper_params_bound <- config$hyper_params_bound

#TEST:
#set.seed(seed)

# Adjusting for 'no_cores' based on the configuration or dynamically
if (config$no_cores == "auto") {
  no_cores <- parallel::detectCores() - 1
} else {
  no_cores <- config$no_cores
}

# Register parallel backend using the determined number of cores
parallel::detectCores()
doParallel::registerDoParallel(cores = no_cores)

####
option_list <- list(
  make_option(c("-n", "--sample_size"), type="integer", default=dp_info$sample_size, 
              help="Sample size", metavar="number"),
  make_option(c("-m", "--mu"), type="numeric", default=1, 
              help="Mean of the distribution", metavar="number"),
  make_option(c("-s", "--sigma_square"), type="numeric", default=1, 
              help="Variance of the distribution", metavar="number"),
  make_option(c("-l", "--lower_bound"), type="numeric", default=dp_info$lower_bound, 
              help="Lower bound for privatization", metavar="number"),
  make_option(c("-u", "--upper_bound"), type="numeric", default=dp_info$upper_bound, 
              help="Upper bound for privatization", metavar="number"),
  make_option(c("-e", "--epsilon"), type="numeric", default=dp_info$epsilon, 
              help="Privacy budget", metavar="number"),
  make_option(c("-d", "--seed"), type="numeric", default=-1, 
              help="Seed for reproducibility of data", metavar="number"),
  make_option(c("-p", "--particles"), type = "numeric", default = 100, help = "Number of particles"),
  make_option(c("-i", "--iter"), type = "numeric", default = 1, help = "Run number")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

# Use the arguments or defaults
n <- args$sample_size
mu <- args$mu
sigma_square <- args$sigma_square
lower_bound <- args$lower_bound
upper_bound <- args$upper_bound
epsilon <- args$epsilon
seed <- args$seed
particles <- args$particles
iter <- args$iter

#################################
#NOTE: For Fixed Setting
#mu = 0
#sigma_square = 1
#lower_bound = -3
#upper_bound = 3
#################################

# Set seed for reproducibility
#if (seed != -1) {
#  set.seed(seed)
#}

if (n != dp_info$sample_size) dp_info$sample_size <- n
if (lower_bound != dp_info$lower_bound) dp_info$lower_bound <- lower_bound
if (upper_bound != dp_info$upper_bound) dp_info$upper_bound <- upper_bound
if (epsilon != dp_info$epsilon) dp_info$epsilon <- epsilon
if (particles != hyper_params$no_particles) hyper_params$no_particles <- particles


clamp <- function(x_data, lower_bound, upper_bound) {
  x_data <- ifelse(x_data < lower_bound, lower_bound, x_data)
  x_data <- ifelse(x_data > upper_bound, upper_bound, x_data)
  return(x_data)
}

privatized_stat.lap <- function(x_data, lower_bound, upper_bound, epsilon, sample_size) {
  clamp_data <- clamp(x_data, lower_bound, upper_bound)

  # NEW:
  # Z-Score Normalization
  norm_clamp_data <- clamp_data / 5
  p <- 0
  deltaa <- p^2 + 4*p + 3
  clamp_sum <- sum(norm_clamp_data) + 
    LaplacesDemon::rlaplace(n = 1, location = 0, scale = deltaa / epsilon)
  clamp_sum_sq <- sum(norm_clamp_data^2) +
    LaplacesDemon::rlaplace(n = 1, location = 0, scale = deltaa / epsilon)

  return(c(clamp_sum, clamp_sum_sq))
}

set.seed(1)
## simulated missing data
missing_data <- rnorm(n = dp_info$sample_size,
                      mean = mu,
                      sd = sqrt(sigma_square))
## s_{dp} in the paper
sdp <- privatized_stat.lap(x_data = missing_data,
                       lower_bound = lower_bound,
                       upper_bound = upper_bound,
                       epsilon = epsilon,
                       sample_size = dp_info$sample_size) #=c(1.004, 0.609)
#TEST:
#set.seed(seed)

#################################
# NOTE: For Fixed Setting
# sdp <- c(0, 1) # Scale the second term
#################################

dp_info$sdp <- sdp

rmodel <- function(theta, dp_info) {
  return(rnorm(n = dp_info$sample_size,
               mean = theta[1],
               sd = sqrt(theta[2])))
}

rprior <- function(hyper_params, particle) {
  theta_prior <- matrix(data = 0,
                        nrow = hyper_params$dim,
                        ncol = 1)
  for (d in 1:hyper_params$dim) {
    if (d == 1){
      theta_prior[d, ] <- rnorm(n = 1,
                                mean = hyper_params$hp[d, 1],
                                sd = hyper_params$hp[d, 2])
    } else if (d == 2){
      theta_prior[d, ] <- LaplacesDemon::rinvgamma(n = 1,
                                                   shape = hyper_params$hp[d, 1],
                                                   scale = 1/hyper_params$hp[d, 2])
    } else {  }
  }

  particle[[1]]$theta <- theta_prior

  return(particle)
}

threshold.lap <- function(privacy_budget, x_data, dp_info) {
  clamp_data <- clamp(x_data = x_data,
                      lower_bound = dp_info$lower_bound,
                      upper_bound = dp_info$upper_bound)

  norm_clamp_data = clamp_data / 5
  clamp_sum = sum(norm_clamp_data)
  clamp_sum_sq = sum(norm_clamp_data^2)

  p <- 0
  deltaa <- p^2 + 4*p + 3
  
  log_ratio <- LaplacesDemon::dlaplace(x = dp_info$sdp[1],
                                       location = clamp_sum,
                                       scale = deltaa / privacy_budget,
                                       log = TRUE) +
    LaplacesDemon::dlaplace(x = dp_info$sdp[2],
                            location = clamp_sum_sq,
                            scale = deltaa / privacy_budget,
                            log = TRUE) -
    (
      LaplacesDemon::dlaplace(x = clamp_sum,
                              location = clamp_sum,
                              scale = deltaa / privacy_budget,
                              log = TRUE) +
        LaplacesDemon::dlaplace(x = clamp_sum_sq,
                                location = clamp_sum_sq,
                                scale = deltaa / privacy_budget,
                                log = TRUE)
    )
  return(log_ratio)
}

# NOTE:
# Generate Prior
# Draw samples from the prior (1 particle at a time)
# rejection

rejection_abc_exact.lap <- function(dp_info, hyper_params, hyper_params_bound) {

  #NOTE: Parallel
  #results <- foreach::foreach(i = 1:hyper_params$no_particles) %dopar% {
  results <- sapply(1:hyper_params$no_particles, function(j) {
    particle <- replicate(n = 1,
                         expr = list(attempt = 1,
                                     theta = 0),
                         simplify = FALSE)
    repeat {
      particle <- rprior(hyper_params, particle) # Prior generating 1 particle
      theta <- particle[[1]]$theta
      
      ### We use rejection sampler instead of a rho measurement
      model_data <- rmodel(theta = theta,
                           dp_info = dp_info)
      r <- threshold.lap(privacy_budget = dp_info$epsilon,
                         x_data = model_data,
                         dp_info = dp_info)
      u <- runif(n = 1,
                 min = 0,
                 max = 1)

      if (r > log(u)) {
        break
      }

      particle[[1]]$attempt <- particle[[1]]$attempt + 1
    }
    return(particle)
  })
  
  return(results)
}

#TEST:
#rejection_abc_ci <- replicate(n = 100,
#                                 expr = 1,
#                                 simplify = FALSE)

#rejection_abc_ci$time <- -1
#rejection_abc_ci$sdp <- dp_info$sdp

#tictoc::tic()
#for (l in 1:100) {
#  rejection_abc_ci[[l]] <- rejection_abc_exact.lap(dp_info = dp_info,
#                              hyper_params = hyper_params,
#                              hyper_params_bound = hyper_params_bound)
#  
#  save(rejection_abc_ci, file = paste("rejection_abc_lap_type_1_intermediate_n", dp_info$sample_size, "e", dp_info$epsilon, "seed", seed, "p", hyper_params$no_particles, "i", iter, ".RData", sep = ""))
#
#}
#time <- tictoc::toc()

tictoc::tic()
rejection_abc_ci <- rejection_abc_exact.lap(dp_info = dp_info,
                                     hyper_params = hyper_params,
                                     hyper_params_bound = hyper_params_bound)
time <- tictoc::toc()

rejection_abc_ci$time <- unname(time$toc - time$tic)
rejection_abc_ci$sdp <- dp_info$sdp

#save(rejection_abc_ci, file = paste("rejection_abc_lap_type_1_final_n", dp_info$sample_size, "e", dp_info$epsilon, "seed", seed, "p", hyper_params$no_particles, "i", iter, ".RData", sep = ""))
save(rejection_abc_ci, file = paste("rejection_abc_lap_type_1_final_n", dp_info$sample_size, "e", dp_info$epsilon, "p", hyper_params$no_particles, "i", iter, ".RData", sep = ""))
