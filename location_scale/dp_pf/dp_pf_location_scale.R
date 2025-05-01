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

config_file_path <- "~/smcDiffPrivacy/inst/extdata/dp_abc_config.json"

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
#lower_bound = -5
#upper_bound = 5

# Set seed for reproducibility
#if (seed != -1) {
#  set.seed(seed)
#}

lower_bound = -5
upper_bound = 5

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
  # Normalize
  # p = 0
  # deltaa = p^2+4*p+3
  # sum of clamp data and clamp data^2

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
# sdp <- c(0, sample_size / 10) # Scale the second term
#################################

dp_info$sdp <- sdp

rmodel <- function(theta_new, dp_info) {
  return(rnorm(n = dp_info$sample_size,
               mean = theta_new[1],
               sd = theta_new[2]))
}

rprior <- function(hyper_params, particles) {
  theta_prior <- matrix(data = 0,
                        nrow = hyper_params$dim,
                        ncol = hyper_params$no_particles)
  for (d in 1:hyper_params$dim) {
    if (d == 1){
      theta_prior[d, ] <- rnorm(n = hyper_params$no_particles,
                                mean = hyper_params$hp[d, 1],
                                sd = hyper_params$hp[d, 2])
    } else if (d == 2){
      theta_prior[d, ] <- LaplacesDemon::rinvgamma(n = hyper_params$no_particles,
                                                   shape = hyper_params$hp[d, 1],
                                                   scale = 1/hyper_params$hp[d, 2])
    } else {  }
  }

  particles <- lapply(seq_along(particles), function(index) {
    particle <- particles[[index]]
    particle$theta_old <- theta_prior[, index]
    return(particle)
  })
}

dprior <- function(d, theta_new, hyper_params) {
  if (d == 1) {
    dst <- dnorm(x = theta_new[d],
                 mean = hyper_params$hp[d, 1],
                 sd = hyper_params$hp[d, 2])
  } else if (d == 2) {
    dst <- LaplacesDemon::dinvgamma(x = theta_new[d],
                     shape = hyper_params$hp[d, 1],
                     scale = 1/hyper_params$hp[d, 2])
  } else {  }
  return(dst)
}

rperturb <- function(t, hyper_params) {
  if (t <= 0) stop("Parameter 't' must be positive")
  variation <- 1 / t
  perturbation <- rnorm(n = hyper_params$dim,
                        mean = 0,
                        sd = variation)
  return(perturbation)
}

dperturb <- function(t, theta_new, theta_old, hyper_params_bound) {
  if (t <= 0) stop("Parameter 't' must be positive")
  variation <- 1 / t
  dst <- dnorm(x = theta_new[1],
               mean = theta_old[1, ],
               sd = variation) *
    dnorm(x = theta_new[2],
          mean = theta_old[2, ],
          sd = variation) /
    pnorm(q = hyper_params_bound[2, 1],
          mean = theta_old[2, ],
          sd = variation,
          lower.tail = FALSE)
  return(dst)
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

weight_update <- function(t, theta_new, theta_old, weight_old,
                              hyper_params, hyper_params_bound, dprior) {
  dprior(d = 1,
         theta_new = theta_new,
         hyper_params = hyper_params) *
    dprior(d = 2,
           theta_new = theta_new,
           hyper_params = hyper_params) /
    sum(weight_old *
          dperturb(t = t,
                       theta_new = theta_new,
                       theta_old = theta_old,
                       hyper_params_bound = hyper_params_bound))
}

verification_fn <- function(picked_particle, hyper_params_bound) {
  verify <- (picked_particle[2] > hyper_params_bound[2, 1])
  return (verify)
}

theta_update <- function(t, particles, hyper_params, hyper_params_bound) {
  index <- sample(x = 1:hyper_params$no_particles,
                  size = 1,
                  replace = TRUE,
                  prob = sapply(particles, function(x) x$weight_old))
  repeat {
    picked_particle <- particles[[index]]$theta_old + rperturb(t, hyper_params)

    # In this case, we only need to examine sigma^2 > 0
    # verify <- (picked_particle[2] > hyper_params_bound[2, 1])
    # if (verify == TRUE) {
    #   break
    # }

    verify = verification_fn(picked_particle = picked_particle,
                     hyper_params_bound = hyper_params_bound)
    if (verify == TRUE) {
      break
    }
  }
  return(picked_particle)
}

schedule <- function(dp_info, hyper_params) {
  sche <- seq(by = dp_info$epsilon / hyper_params$stop_time,
              from = dp_info$epsilon / hyper_params$stop_time,
              to = dp_info$epsilon)
  return(sche)
}

abc_smc_exact.lap <- function(dp_info, hyper_params, hyper_params_bound) {
  ### S1: Schedule computation

  ep_schedule <- schedule(dp_info = dp_info, hyper_params = hyper_params)

  particles <- replicate(n = hyper_params$no_particles,
                         expr = list(weight_old = 1 / hyper_params$no_particles,
                                     attempt = 1,
                                     theta_old = 0,
                                     weight_new = 1),
                         simplify = FALSE)
  particles <- rprior(hyper_params, particles)

  ### S2: Sequential processing
  t <- 1
  step_result_store <- replicate(n = hyper_params$stop_time + 1,
                                 expr = particles,
                                 simplify = FALSE)

  repeat {
    theta_old <- sapply(particles, function(x) x$theta_old)
    weight_old <- sapply(particles, function(x) x$weight_old)
    # Parallelized computation for each particle
    results <- foreach::foreach(i = 1:hyper_params$no_particles) %dopar% {
      a <- 1
      repeat {
        theta_new <- theta_update(t = t,
                                      particles = particles,
                                      hyper_params = hyper_params,
                                      hyper_params_bound = hyper_params_bound)

        ### We use rejection sampler instead of a rho measurement
        model_data <- rmodel(theta_new = theta_new,
                                 dp_info = dp_info)
        r <- threshold.lap(privacy_budget = ep_schedule[t],
                           x_data = model_data,
                           dp_info = dp_info)
        u <- runif(n = 1,
                   min = 0,
                   max = 1)

        if (r > log(u)) {
          break
        }
        a <- a + 1
      }

      weight_new <- weight_update(t = t,
                                      theta_new = theta_new,
                                      theta_old = theta_old,
                                      weight_old = weight_old,
                                      hyper_params = hyper_params,
                                      hyper_params_bound = hyper_params_bound,
                                      dprior = dprior)

      list(theta_new = theta_new,
           weight_new = weight_new,
           attempt = a)
    }

    # Update particles with results
    for (i in 1:hyper_params$no_particles) {
      particles[[i]]$theta_old <- results[[i]]$theta_new
      particles[[i]]$weight_new <- results[[i]]$weight_new
      particles[[i]]$attempt <- results[[i]]$attempt
    }

    sum_of_weights <- sum(sapply(particles, function(y) y$weight_new))

    particles <- lapply(particles, function(x) {
      x$weight_old <- x$weight_new / sum_of_weights
      return(x)
    })

    step_result_store[[t + 1]] <- particles

    ### S3: Check stop condition
    if (t >= hyper_params$stop_time) {
      break
    }
    print(t)
    t <- t + 1
  }

  print(t)
  return(step_result_store)
}

#TEST:
#tictoc::tic()
#smc_abc_ci <- lapply(1:100,
#                     function(l) {
#                       abc_smc_exact.lap(dp_info = dp_info,
#                                     hyper_params = hyper_params,
#                                     hyper_params_bound = hyper_params_bound)
#                     })
#time <- tictoc::toc()

tictoc::tic()
smc_abc_ci <- abc_smc_exact.lap(dp_info = dp_info,
                                hyper_params = hyper_params,
                                hyper_params_bound = hyper_params_bound)
time <- tictoc::toc()

smc_abc_ci$time <- unname(time$toc - time$tic)
smc_abc_ci$sdp <- dp_info$sdp
#save(smc_abc_ci, file = paste("abc_smc_lap_type_1_final_n", dp_info$sample_size, "e", dp_info$epsilon, "seed", seed, "p", hyper_params$no_particles, "i", iter, ".RData", sep = ""))
#TEST:
save(smc_abc_ci, file = paste("data/abc_smc_lap_type_1_final_n", dp_info$sample_size, "e", dp_info$epsilon, "p", hyper_params$no_particles, "i", iter, ".RData", sep = ""))

