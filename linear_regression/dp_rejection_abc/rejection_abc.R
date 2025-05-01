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

no_cores <- 100
max_cores = parallel::detectCores()
if (no_cores <= max_cores) {
  doParallel::registerDoParallel(cores = no_cores)
} else {
  doParallel::registerDoParallel(cores = max_cores)
}

# Define options for command line arguments
option_list <- list(
  make_option(c("-n", "--sample_size"), type="integer", default=200, 
              help="Sample size", metavar="number"),
  make_option(c("-e", "--epsilon"), type="numeric", default=1, 
              help="Privacy budget", metavar="number"),
  make_option(c("-d", "--seed"), type="numeric", default=731, 
              help="Seed for reproducibility of data", metavar="number"),
  make_option(c("-p", "--particles"), type="numeric", default=100, 
              help="Particles", metavar="number"),
  make_option(c("-i", "--iter"), type="numeric", default=1, 
              help="Run Number", metavar="number")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

# Use the arguments or defaults
n <- args$sample_size
epsilon <- args$epsilon
seed <- args$seed
iter <- args$iter
N <- args$particles

# Random sdp setting
p = 1
mu = 1
Phi = 1
tau = 1
beta = c(0, 2) 

# Fixed sdp setting
# p = 1
# mu = 0
# Phi = 1
# tau = 1
# beta = c(0, 0)

deltaa = p^2 + 4*p + 3

df = 2
theta = 0
Sigma = 1
sc = 1.25
sh = 2
a = 0
B = 1
m = c(0,0)
V = diag(p+1)

L = -5
U = 5


clamp <- function(x_data, L, U) {
  x_data <- ifelse(x_data < L, L, x_data)
  x_data <- ifelse(x_data > U, U, x_data)
  return(x_data)
}

privatized_stat.lap <- function(x_data, L, U) {
  clamp_data <- clamp(x_data, L, U)  
  
  ydp <- clamp_data[, 1, drop = FALSE] 
  xdp <- cbind(1, clamp_data[,-1, drop = FALSE]) 
  s1 <- t(xdp) %*% ydp 
  s2 <- t(ydp) %*% ydp 
  s3 <- t(xdp) %*% xdp 
  
  ur_s1 <- c(s1) 
  ur_s2 <- c(s2) 
  ur_s3 <- s3[upper.tri(s3, diag = TRUE)][-1] 
  
  return(c(ur_s1,ur_s2,ur_s3) + rlaplace(((p+2)^2+p)/2, location = 0, scale = deltaa/epsilon))
}

# Random sdp generation
set.seed(1)
missingX <- mvrnorm(n, mu = mu, Sigma = 1/Phi) #model matrix generation~N(mu, 1) ##
missingY <- cbind(1, missingX) %*% beta + rnorm(n, sd = sqrt(1/tau)) #model error~N(0, 1)

sdp <- privatized_stat.lap(cbind(missingY, missingX), L=L, U=U) 

# Fixed sdp generation
# sdp <- (0, 0, n/10, 0, n/10)

set.seed(seed)

dp_info = list(sdp = sdp,
               L=L,
               U=U,
               epsilon=epsilon,
               n=n,
               p=p,
               deltaa=deltaa)

hp = list(df = df,
          theta = theta,
          Sigma = Sigma,
          sc = sc,
          sh = sh,
          a = a,
          B = B,
          m = m,
          V = V)


hyper_params = list(hp = hp,
                    stopTime = 10,
                    N=N,
                    dim=5)

params_bound = list(mu = c(-Inf, Inf),
                    Phi = c(0, Inf),
                    tau = c(0, Inf),
                    beta = rbind(c(-Inf, Inf), c(-Inf, Inf)))

rmodel <- function(mu, Phi, tau, beta, dp_info){
  daX = mvrnorm(dp_info$n, mu = mu, Sigma = 1/Phi)
  daY = cbind(1, daX) %*% beta[1:2] + rnorm(dp_info$n, sd = sqrt(1/tau))
  
  return(cbind(daY, daX))
}

rprior <- function(hyper_params, particles){ #OK
  
  muPrior = t(rmvt(n = 1, mu = hyper_params$hp$theta, S = hyper_params$hp$Sigma, df = hyper_params$hp$df))
  PhiPrior = t(abs(rmvt(n = 1, mu = hyper_params$hp$theta, S = hyper_params$hp$Sigma, df = hyper_params$hp$df)))
  tauPrior = t(rweibull(n = 1, scale = hyper_params$hp$sc, shape = hyper_params$hp$sh))
  betaPrior = rmvt(n = 1, mu = hyper_params$hp$m, S = hyper_params$hp$V, df = hyper_params$hp$df)
    

  particles[[1]]$mu <- muPrior
  particles[[1]]$Phi <- PhiPrior
  particles[[1]]$tau <- tauPrior
  particles[[1]]$beta <- betaPrior
    
  return(particles)
}

threshold <- function(privacyBudget, Data, dp_info){ #OK
  
  cData <- clamp(Data, dp_info$L, dp_info$U) 
  ydp <- cData[,1, drop = FALSE] 
  xdp <- cbind(1, cData[,-1, drop = FALSE]) 
  s1 <- t(xdp) %*% ydp 
  s2 <- t(ydp) %*% ydp 
  s3 <- t(xdp) %*% xdp 
  s3 <- s3[upper.tri(s3, diag = TRUE)][-1]
  clampStat = c(s1, s2, s3)
  
  #log ratio
  Delta = dp_info$deltaa
  r = sum(dlaplace(x = dp_info$sdp, location = clampStat, scale = Delta / privacyBudget, log = T)) - sum(dlaplace(x = dp_info$sdp, location = dp_info$sdp, scale = Delta / privacyBudget, log = T))
  return(r)
}

rejection_lr_abc_exact.lap <- function(dp_info, hyper_params, hyper_params_bound) {
  particle <- replicate(n = 1,
                         expr = list(attempt = 1,
                                     mu = 0,
                                     Phi = 0,
                                     tau = 0,
                                     beta = 0),
                         simplify = FALSE)

  step_result_store <- replicate(n = hyper_params$N,
                                 expr = particle,
                                 simplify = FALSE)

  results <- foreach(i=1:hyper_params$N) %dopar% {
    repeat {
      particle <- rprior(hyper_params, particle) # Prior generating 1 particle
      mu <- particle[[1]]$mu
      Phi <- particle[[1]]$Phi
      tau <- particle[[1]]$tau
      beta <- particle[[1]]$beta
      
      ### We use rejection sampler instead of a rho measurement
      model_data <- rmodel(mu = mu, Phi = Phi, tau = tau, beta = beta,
                           dp_info = dp_info)
      r <- threshold(privacyBudget = dp_info$epsilon,
                         Data = model_data,
                         dp_info = dp_info)
      u <- runif(n = 1,
                 min = 0,
                 max = 1)

      if (r > log(u)) {
        break
      }

      particle[[1]]$attempt <- particle[[1]]$attempt + 1
    }
    #step_result_store[[i]] <- particle
    return(particle)
  }

  return(results)
  #return(step_result_store)
}

# rejection_lr_abc_ci <- replicate(n = 100,
#                                  expr = 1,
#                                  simplify = FALSE)

#rejection_lr_abc_ci$time <- -1
#rejection_lr_abc_ci$sdp <- dp_info$sdp

#tictoc::tic()
#for (l in 1:100) {
#rejection_lr_abc_ci[[l]] <- rejection_lr_abc_exact.lap(dp_info = dp_info,
#                              hyper_params = hyper_params,
#                              hyper_params_bound = params_bound)
  
# save(rejection_lr_abc_ci, file = paste("rejection_lr_abc_lap_type_1_intermediate_n", dp_info$n, "e", dp_info$epsilon, "seed", seed, "p", hyper_params$N, "i", iter, ".RData", sep = ""))

#}
#time <- tictoc::toc()

tictoc::tic()
rejection_abc_ci <- rejection_abc_exact.lap(dp_info = dp_info,
                                     hyper_params = hyper_params,
                                     hyper_params_bound = hyper_params_bound)
time <- tictoc::toc()

rejection_lr_abc_ci$time <- unname(time$toc - time$tic)
rejection_lr_abc_ci$sdp <- dp_info$sdp

save(rejection_lr_abc_ci, file = paste("linear_regression/rejection_abc/data/rejection_lr_abc_lap_type_3_final_n", dp_info$n, "e", dp_info$epsilon, "seed", seed, "p", hyper_params$N, "i", iter, ".RData", sep = ""))

