library(MASS)
library(VGAM)
library(emdbook)
library(LaplacesDemon)
library(truncnorm)
library(tictoc)
library(foreach)
library(doParallel)
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

n <- args$sample_size
epsilon <- args$epsilon
seed <- args$seed
iter <- args$iter
N <- args$particles

#set.seed(seed)

# Random sdp setting
#p = 1
#mu = 1
#Phi = 1
#tau = 1
#beta =  c(0, 2) 

# Fixed sdp setting
p = 1
mu = 0
Phi = 1
tau = 1
beta = c(0, 0)

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

clamp.data <- function(Data, L, U) { pmin(pmax(Data, L), U) / (U/2 - L/2) }

privatized.stat <- function(Data, L, U){ 
  
  cData <- clamp.data(Data, L, U) 
  ydp <- cData[,1, drop = FALSE] 
  xdp <- cbind(1, cData[,-1, drop = FALSE]) 
  s1 <- t(xdp) %*% ydp 
  s2 <- t(ydp) %*% ydp 
  s3 <- t(xdp) %*% xdp 
  
  ur_s1 <- c(s1) 
  ur_s2 <- c(s2) 
  ur_s3 <- s3[upper.tri(s3, diag = TRUE)][-1] 
  
  return( c(ur_s1,ur_s2,ur_s3) + rlaplace(((p+2)^2+p)/2, location = 0, scale = deltaa/epsilon) )
}

# Random sdp generation
# set.seed(1)
# missingX <- mvrnorm(n, mu = mu, Sigma = 1/Phi) #model matrix generation~N(mu, 1) ##
# missingY <- cbind(1, missingX) %*% beta + rnorm(n, sd = sqrt(1/tau)) #model error~N(0, 1)

# sdp <- privatized.stat(cbind(missingY, missingX), L=-5, U=5) 

# Fixed sdp generation
sdp <- c(0, 0, n/10, 0, n/10)

# set.seed(seed)
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


#####################################################################################################################################

rmodel <- function(mu, Phi, tau, beta, dp_info){ #OK
  daX = mvrnorm(dp_info$n, mu = mu, Sigma = 1/Phi)
  daY = cbind(1, daX) %*% beta[1:2] + rnorm(dp_info$n, sd = sqrt(1/tau))
  
  return(cbind(daY, daX))
}

threshold <- function(privacyBudget, Data, dp_info){ #OK
  
  cData <- clamp.data(Data, dp_info$L, dp_info$U) 
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

dprior <- function(thetaNew, hyper_params){ #OK
  
  logDst = dmvt(x = thetaNew$mu, mu = hyper_params$hp$theta, S = matrix(hyper_params$hp$Sigma), df = hyper_params$hp$df, log = T) +
    dmvt(x = thetaNew$Phi, mu = hyper_params$hp$a, S = matrix(hyper_params$hp$B), df = hyper_params$hp$df, log = T) + log(2) +
    dweibull(x = thetaNew$tau, scale = hyper_params$hp$sc, shape = hyper_params$hp$sh, log = T) +
    dmvt(x = thetaNew$beta, mu = hyper_params$hp$m, S = hyper_params$hp$V, df = hyper_params$hp$df, log = T)
  return(logDst)
}

rprior <- function(hyper_params, particles){ #OK
  
  muPrior = t(rmvt(n = hyper_params$N, mu = hyper_params$hp$theta, S = hyper_params$hp$Sigma, df = hyper_params$hp$df))
  PhiPrior = t(abs(rmvt(n = hyper_params$N, mu = hyper_params$hp$theta, S = hyper_params$hp$Sigma, df = hyper_params$hp$df)))
  tauPrior = t(rweibull(n = hyper_params$N, scale = hyper_params$hp$sc, shape = hyper_params$hp$sh))
  betaPrior = sapply(1:hyper_params$N, function(k){
    rmvt(n = 1, mu = hyper_params$hp$m, S = hyper_params$hp$V, df = hyper_params$hp$df)
  })
    
  particles <- lapply(seq_along(particles), function(index) {
    particle <- particles[[index]]

    particle$mu <- muPrior[, index]
    particle$Phi <- PhiPrior[, index]
    particle$tau <- tauPrior[, index]
    particle$beta <- betaPrior[, index]
    
    return(particle)
  })
}

rperturb <- function(t, mu, Phi, tau, beta, hyper_params, params_bound){ #OK
  variation = 1/t
  
  muNew = mvrnorm(n = 1, mu = mu, Sigma = variation*diag(1,1))
  log_PhiNew = mvrnorm(n = 1, mu = log(Phi), Sigma = variation*diag(1,1))
  PhiNew = exp(log_PhiNew)
  log_tauNew = mvrnorm(n = 1, mu = log(tau), Sigma = variation*diag(1,1))
  tauNew = exp(log_tauNew)
  betaNew = mvrnorm(n = 1, mu = beta, Sigma = variation*diag(1,2))
  
  perturbed = list(mu = muNew, Phi = PhiNew, tau = tauNew, beta = betaNew)
  return(perturbed)
}

dperturb <- function(t, thetaNew, muOld, PhiOld, tauOld, betaOld, hyper_params, params_bound){
  variation = 1/t
  
  logDst = sapply(1:hyper_params$N, function(j){
    dmvnorm(x = thetaNew$mu, mu = muOld[j], Sigma = variation*diag(1,1), log = T) +
      dmvnorm(x = log(thetaNew$Phi), mu = log(PhiOld[j]), Sigma = variation*diag(1,1), log = T) +  
      dmvnorm(x = log(thetaNew$tau), mu = log(tauOld[j]), Sigma = variation*diag(1,1), log = T) +  
      dmvnorm(x = thetaNew$beta[1:2], mu = betaOld[1:2,j], Sigma = variation*diag(1,2), log = T)
  })

  return(logDst)
}

theta.update <- function(t, particles, hyper_params, params_bound){ #OK
  index = sample((1:hyper_params$N), size = 1, replace = T, prob = sapply(particles, function(x) x$weightOld))
  pickedParticle = rperturb(t=t, mu = particles[[index]]$mu, Phi = particles[[index]]$Phi, tau = particles[[index]]$tau, beta = particles[[index]]$beta, hyper_params = hyper_params, params_bound = params_bound)
  return(pickedParticle)
}

weight.update <- function(t, particles, thetaNew, hyper_params, params_bound, dprior){ #OK
  muOld = sapply(particles, function(x) x$mu)
  PhiOld = sapply(particles, function(x) x$Phi)
  tauOld = sapply(particles, function(x) x$tau)
  betaOld = sapply(particles, function(x) x$beta)
  exp(dprior(thetaNew = thetaNew, hyper_params = hyper_params) - log(sum(exp(log(sapply(particles, function(x) x$weightOld)) + dperturb(t=t, thetaNew = thetaNew, muOld = muOld, PhiOld = PhiOld, tauOld = tauOld, betaOld = betaOld, hyper_params = hyper_params, params_bound = params_bound)))))
}

schedule <- function(dp_info, hyper_params){ #OK
  sche = seq(by=dp_info$epsilon/hyper_params$stopTime, from=dp_info$epsilon/hyper_params$stopTime, to=dp_info$epsilon)
  return(sche)
}

ABC.SMC.Exact <- function(dp_info, hyper_params, params_bound){ 
  
  stepResultStore <- list()
  
  ### S1
  
  # outputs a vector
  epSchedule = schedule(dp_info = dp_info, hyper_params = hyper_params)
  
  particles <- replicate(hyper_params$N, list(mu = 0, Phi = 0, tau = 0, beta = 0, weightOld = 1/hyper_params$N, attempt = 1, weightNew = 1, thetaNew = 1), simplify = FALSE)
  particles <- rprior(hyper_params, particles)
  
  ### S2
  t = 1
  stepResultStore[[t]] = particles
  i=1
  repeat {
    results <- foreach(i=1:hyper_params$N) %dopar% {
      a = 1
      ### S2.1
      repeat {
        thetaNew = theta.update(t=t, particles = particles, hyper_params = hyper_params, params_bound = params_bound)
        modelData = rmodel(mu = thetaNew$mu, Phi = thetaNew$Phi, tau = thetaNew$tau, beta = thetaNew$beta, dp_info = dp_info)
        r = threshold(privacyBudget = epSchedule[t], Data = modelData, dp_info = dp_info)
        u = runif(n = 1, min = 0, max = 1)
        #c
        if (r > log(u)) {
          break
        }
        a = a + 1
      }

      weightNew = weight.update(t=t, thetaNew = thetaNew, particles = particles, hyper_params = hyper_params, params_bound = params_bound, dprior = dprior)

      list(weightNew = weightNew, attempt = a, mu = thetaNew$mu, Phi = thetaNew$Phi, tau = thetaNew$tau, beta = thetaNew$beta)
    }

    for (i in 1:hyper_params$N) {
      particles[[i]]$weightNew = results[[i]]$weightNew
      particles[[i]]$mu = results[[i]]$mu
      particles[[i]]$Phi = results[[i]]$Phi
      particles[[i]]$tau = results[[i]]$tau
      particles[[i]]$beta = results[[i]]$beta
    }

    sum_of_weights = sum(sapply(particles, function(y) y$weightNew))

    particles <- lapply(particles, function(x) {
        x$weightOld = x$weightNew / sum_of_weights
        return(x)
    })
    
    stepResultStore[[t+1]] = particles
    
    ### S3
    if (t >= hyper_params$stopTime) {
      break
    }
    print(t)
    t = t + 1
  }
  
  print(t)
  return(stepResultStore)
}

tic()
res <- ABC.SMC.Exact(dp_info = dp_info, hyper_params = hyper_params, params_bound = params_bound)
time <- toc()

res$time <- unname(time$toc - time$tic)
res$sdp <- dp_info$sdp

print("Saving results")
save(res, file = paste("linear_regression/dp_pf/non_conjugate_prior/results/lr_nonconj_n", n, "e", epsilon, "seed", seed, "p", hyper_params$N, "i", iter, ".RData", sep = ""))
print("Saved Results")