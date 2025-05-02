library(dplyr)
library(emdbook)
library(MASS)
library(progressr)
library(tictoc)
library(doParallel)
library(foreach)
library(optparse)
  

no_cores <- 100
max_cores = parallel::detectCores()
if (no_cores <= max_cores) {
  doParallel::registerDoParallel(cores = no_cores)
} else {
  doParallel::registerDoParallel(cores = max_cores)
}
no_chains = 1


option_list <- list(
  make_option(c("--shape"), type="integer", default=6,
              help="Shape for Gamma Dist", metavar="number"),
  make_option(c("--scale"), type="integer", default=4,
              help="Scale for Gamma Dist", metavar="number"),
  make_option(c("--sigma1"), type="integer", default=16,
              help="Variance for 1st", metavar="number"),
  make_option(c("--sigma2"), type="integer", default=16,
              help="Variance for 2nd", metavar="number"),
  make_option(c("--id"), type="integer", default=1,
              help="Chain ID", metavar="number"),
  make_option(c("--epsilon"), type="numeric", default=0.5,
              help="Privacy Budget", metavar="number")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

# Use the arguments or defaults
sigma_square_1 <- args$sigma1
sigma_square_2 <- args$sigma2
shape <- args$shape
scale <- args$scale
epsilon <- args$epsilon
chain_id <- args$id

## Select and store the data of Vancouver
# rawData <- read.csv("data_donnees_2021_ind_v2.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
#vcrData <- rawData[rawData$CMA==933,]
#save(vcrData, file = "vcrData.RData")
#load(file = "vcrData.RData")
# vicData <- rawData[rawData$CMA==935,]
#save(vicData, file = "vicData.RData")
load(file = "vicData.RData")


## Choose variables
vicMat <- vicData %>% filter(Retir != 88888888, AGEGRP != 88)
vicMat <- vicMat %>% dplyr::select(Retir, AGEGRP)
vicMat$Retir <- ifelse(vicMat$Retir<99999999, 1, 0)
vicMat$AGEGRP <- ifelse(vicMat$AGEGRP<8, vicMat$AGEGRP*2.5+1, vicMat$AGEGRP*5-18)

# hist(vicMat[vicMat$Retir==0,][,2])
# hist(vicMat[vicMat$Retir==1,][,2])

vicMat = vicMat[2001:3000,]

## Pre-process data for DP
zAGEGRP <- (vicMat$AGEGRP - min(vicMat$AGEGRP))/(max(vicMat$AGEGRP) - min(vicMat$AGEGRP)) 
xAGEGRP <- 2*zAGEGRP - 1

## EDA: true beta and true (sum z_i, sum z_i^2)
glm(vicMat$Retir ~ xAGEGRP, family = binomial) #beta = (-3.950,  7.304) #### to check the non-private beta
c(sum(zAGEGRP), sum(zAGEGRP^2)) #(5021.79, 3188.06) #### to check the non-private (sum zi, sum zi^2)

# lower and upper bound of tightest K-Norm for (sum z_i, sum z_i^2) sensitivity

# lower = function(x){
#   bottom = (x+1)^2-1
#   middle = x-1/4
#   top = x^2
#   return(bottom*(x<=-1/2)+middle*(x>-1/2 & x<1/2)+top*(x>=1/2))
# }
# 
# upper = function(x){
#   bottom = -x^2
#   middle = x+1/4
#   top = -(x-1)^2+1
#   return(bottom*(x<=-1/2)+middle*(x>-1/2 & x<1/2)+top*(x>=1/2))
# }

expit = function(x){
  return(exp(x)/(1+exp(x)))
}

expitit = function(x){
  return(exp(x)/(1+exp(x))^2)
}

# Infty-Norm extended Objective Perturbation

ObjPertLog = function(ep,norm,q,X,Y){
  m = dim(X)[2]
  n = dim(X)[1]
  lambda = (1/4)*(m)
  xi=1
  b=rep(0,m)
  if(ep==0){
    gamma=0
  } else {
    if(norm==1){
      xi = (m)*2
    }
    if(norm==2){
      xi = sqrt(m)*2
    }
    if(norm==3){# 3 represents infinity.
      xi = 2
    }
    
    ### Infty-Norm mechanism
    U1 = runif(2,min=-1,max=1)
    G1 = rgamma(1,shape=3,rate = ep/1) #rate=epsilon/sensitivity
    N1 = G1*U1
    ###
    
    b = (xi/(q*ep))*N1
    
    gamma = lambda/(exp((1-q)*ep)-1)
  }
  
  
  obj = function(beta){
    return( -(1/n)*sum(Y*(X%*%beta) - log(1+exp(X%*%beta)))+gamma/(2*n)*sum(beta^2) + sum(b*beta)/n)
  }
  grad = function(beta){
    logLike = (apply((X)*matrix(Y-exp(X%*%beta)/(1+exp(X%*%beta)),nrow=n,ncol=m),2,function(x) sum(x)))
    # 2 means cols
    return( -(1/n)*logLike    + (gamma/n)*beta + b/n)
  }
  #min = optim(par = rep(0,m),fn=obj,gr = grad,method = "L-BFGS-B")
  min = optim(par = rep(0,m),fn=obj,gr = grad,method = "BFGS")
  
  if(min$convergence!=0){
    print("Obj-Pert did not converge") 
  }
  return(min$par)
}###   END OBJ  PERT

# V and nabla V for Objective Perturbation and rejection threshold (V for numerator and nabla V for denominator's bound)
V <- function(ep, q, X, Y, beta) {
  m <- ncol(X)
  n <- nrow(X)
  lambda <- m / 4
  gamma <- lambda / (exp((1 - q) * ep) - 1)
  
  # Precompute the logistic term
  exp_term <- exp(X %*% beta)
  logistic_term <- Y - exp_term / (1 + exp_term)
  # Compute the gradient term directly
  gradient_term <- as.vector(t(X) %*% logistic_term)
  # Compute the final result
  -(gradient_term + gamma * beta)
}

nablaV <- function(ep, q, X, beta) {
  m <- ncol(X)
  n <- nrow(X)
  lambda <- m / 4
  gamma <- lambda / (exp((1 - q) * ep) - 1)

  # Compute the logistic scaling factor for all rows at once
  exp_term <- exp(X %*% beta)
  scaling_factors <- exp_term / (1 + exp_term)^2
  # Compute the weighted sum of outer products
  tempMat <- t(X) %*% sweep(X, 1, scaling_factors, "*")
  
  # Add gamma to eigenvalues
  eigen(tempMat, symmetric = TRUE, only.values = TRUE)$values + gamma
}

# K-Norm mechanism with Rejection Sampling for (sum z_i, sum z_i^2) sensitivity

KNorm <- function(ep, X, rejection_num=100){
  
  # General K-Norm
  # W2 = matrix(runif(rejection_num*2*1,min=-1,max=1),nrow=rejection_num*1,ncol=2)
  # index = which(W2[,2]>=lower(W2[,1]) & W2[,2]<=upper(W2[,1]))
  # 
  # U2 = as.numeric(W2[index[1],])
  # G2 = rgamma(1,shape=3,rate = ep/1)
  # N2 = G2*U2
  
  #l_infty-Norm
  U1 = runif(2,min=-1,max=1)
  G1 = rgamma(1,shape=3,rate = ep/1)
  N1 = G1*U1
  
  z = (X[,-1]+1)/2
  Tz = c(sum(z), sum(z^2))
  
  return(Tz + N1)
}

# Real data analysis settings
dp_info = list(sdp = 0,
               L=-1,
               U=1,
               epsilon=epsilon, #### if it's 0.2, then should go back to lines 165,166 to distribute the budget (90% for betaDP and 10% for zDP)
               q=0.5,
               n=dim(vicMat)[1],#n=10473
               p=1)

hp = list(mu0 = 0,
          sigmaSq0 = sigma_square_1,
          mu1 = 0,
          sigmaSq1 = sigma_square_2,
          aAlpha = shape, ### 9 or 6
          aRate = scale, ### 6 or 4
          bAlpha = shape,
          bRate = scale)


hyper_params = list(hp = hp,
                    stopTime = 10,
                    N=100/no_chains,## at least 100, or ESS larger than 50
                    dim=4)

params_bound = list(beta0 = c(-Inf, Inf),
                    beta1 = c(-Inf, Inf),
                    a = c(0, Inf),
                    b = c(0, Inf))

# Generate sdp of (beta0, beta1, a, b) with use of extended Objective Perturbation and K-Norm mechanism

set.seed(49)
Y = vicMat$Retir
X = cbind(1, xAGEGRP)

# betaDP = ObjPertLog(ep = 0.09, norm = 3,q=0.5, X = X, Y = Y) #(-2.976705,  6.334230)
# zDP = KNorm(ep = 0.01, X = X) #(4996.009, 3174.756)
# sdp = c(betaDP, zDP)

# dp_info$sdp = sdp

sdps <- sapply(1:100, function(i){
  set.seed(i)
  betaDP = ObjPertLog(ep = dp_info$epsilon * 0.9, norm = 3,q=dp_info$q, X = X, Y = Y) #(-2.976705,  6.334230)
  zDP = KNorm(ep = dp_info$epsilon * 0.1, X = X) #(4996.009, 3174.756)

  return(c(betaDP, zDP))
})

set.seed(NULL)
dp_info$sdp = apply(sdps, 1, median)
print(dp_info$sdp)

#####################################################################################################################################

rmodel <- function(i, thetaNew, dp_info){ #OK
  daZ = rbeta(dp_info$n, shape1 = thetaNew[3,i], shape2 = thetaNew[4,i])
  daX = 2*daZ - 1
  daY = sapply(expit(cbind(1, daX) %*% thetaNew[1:2,i]), function(p){rbinom(n=1, size=1, prob=p)})
  
  return(cbind(daY, daX))
}

rmodel.par <- function(newVec, dp_info) {
  # Sample from the beta distribution
  daZ <- rbeta(dp_info$n, shape1 = newVec[3], shape2 = newVec[4])
  
  # Transform daZ to daX
  daX <- 2 * daZ - 1
  
  # Compute probabilities using the logistic function (expit)
  linear_combination <- newVec[1] + newVec[2] * daX  # Avoid explicit cbind
  probabilities <- expit(linear_combination)
  
  # Generate binary outcomes in a vectorized way
  daY <- rbinom(n = dp_info$n, size = 1, prob = probabilities)
  
  # Combine daY and daX into a matrix
  return(cbind(daY, daX))
}

threshold <- function(privacyBudget, Data, dp_info){ #OK
  
  # sdp info
  ept = privacyBudget
  betaDP = dp_info$sdp[1:2]
  zDP = dp_info$sdp[3:4]
  
  # synthetic data info (generated from rmodel)
  cData <- Data
  ySyn <- cData[,1, drop = FALSE] 
  xSyn <- cbind(1, cData[,-1, drop = FALSE])
  
  #beta
  #q=0.85 
  q = ifelse(ept<dp_info$epsilon, 0.01, dp_info$q) #### schedule for q: NEEDS ITS OWN schedule function of t
  m = dim(Data)[2]
  n = dim(Data)[1]
  lambda = m/4
  #epBeta = dp_info$epsilon*0.9
  
  ep1 = ept*0.9
  gamma = lambda/(exp((1-q)*ep1)-1)
  numeratorExponenet = -max(abs(V(ep=ep1,q=q,X=xSyn,Y=ySyn,beta=betaDP)))*ep1 
  denominator = prod(nablaV(ep=ep1,q=q,X=xSyn,beta=betaDP))/(gamma^m)
  logr1 = numeratorExponenet - log(denominator)
  
  #Z
  ep2 = ept*0.1
  zSyn = (xSyn[,-1, drop = FALSE]+1)/2
  Tz = c(sum(zSyn), sum(zSyn^2))
  logr2 = -max(abs(zDP - Tz))*ep2/1 - 0
  
  #log ratio
  r = logr1 + logr2
  return(r)
}

dprior <- function(i, thetaNew, hyper_params){ #OK
  
  logDst = dnorm(x=thetaNew[1,i], mean = hyper_params$hp$mu0, sd = sqrt(hyper_params$hp$sigmaSq0), log = T) +  
            dnorm(x=thetaNew[2,i], mean = hyper_params$hp$mu1, sd = sqrt(hyper_params$hp$sigmaSq1), log = T) +  
              dgamma(x = thetaNew[3,i], shape = hyper_params$hp$aAlpha, rate = hyper_params$hp$aRate, log = T) + 
              dgamma(x = thetaNew[4,i], shape = hyper_params$hp$bAlpha, rate = hyper_params$hp$bRate, log = T)
  
  return(logDst)
}

rprior <- function(hyper_params){ #OK
  
  thetaPrior = matrix(0, nrow = hyper_params$dim, ncol = hyper_params$N)
  
  beta0Prior = t(rnorm(n = hyper_params$N, mean = hyper_params$hp$mu0, sd = sqrt(hyper_params$hp$sigmaSq0)))
  beta1Prior = t(rnorm(n = hyper_params$N, mean = hyper_params$hp$mu1, sd = sqrt(hyper_params$hp$sigmaSq1)))
  aPrior = t(rgamma(n = hyper_params$N, shape = hyper_params$hp$aAlpha, rate = hyper_params$hp$aRate))
  bPrior = t(rgamma(n = hyper_params$N, shape = hyper_params$hp$bAlpha, rate = hyper_params$hp$bRate))

  thetaPrior = rbind(beta0Prior, beta1Prior, aPrior, bPrior)
  
  return(thetaPrior)
}

rperturb <- function(t, index, thetaOld, hyper_params, params_bound){ #OK
  variation = 1/t
  
  beta01abNew = mvrnorm(n=1, mu = c(thetaOld[1:2,index], log(thetaOld[3:4,index])), Sigma = variation*diag(1,4))
  perturbed = c(beta01abNew[1:2], exp(beta01abNew[3:4]))
  
  return(perturbed)
}

dperturb <- function(t, i, thetaNew, thetaOld, hyper_params, params_bound){
  variation = 1/t
  
  logDst = sapply(1:hyper_params$N, function(j){
    dmvnorm(x = thetaNew[1,i], mu = thetaOld[1,j], Sigma = variation*diag(1,1), log = T) +
      dmvnorm(x = thetaNew[2,i], mu = thetaOld[2,j], Sigma = variation*diag(1,1), log = T) +
      dmvnorm(x = log(thetaNew[3,i]), mu = log(thetaOld[3,j]), Sigma = variation*diag(1,1), log = T) +  
      dmvnorm(x = log(thetaNew[4,i]), mu = log(thetaOld[4,j]), Sigma = variation*diag(1,1), log = T)
  })
  
  return(logDst)
}

theta.update <- function(t, thetaOld, weightOld, hyper_params, params_bound){ #OK
  index = sample((1:hyper_params$N), size = 1, replace = T, prob = weightOld)
  pickedParticle = rperturb(t=t, index=index, thetaOld = thetaOld, hyper_params = hyper_params, params_bound = params_bound)
  return(pickedParticle)
}

weight.update <- function(t, i, thetaNew, thetaOld, weightOld, hyper_params, params_bound, dprior){ #OK
  exp(dprior(i=i, thetaNew = thetaNew, hyper_params = hyper_params) - log(sum(exp(log(weightOld) + dperturb(t=t, i=i, thetaNew = thetaNew, thetaOld = thetaOld, hyper_params = hyper_params, params_bound = params_bound)))))
}

schedule <- function(dp_info, hyper_params){ #OK
  sche = seq(by=dp_info$epsilon/hyper_params$stopTime, from=dp_info$epsilon/hyper_params$stopTime, to=dp_info$epsilon)
  #sche = rev(2^(-(0:(hyper_params$stopTime-1))))*dp_info$epsilon
  #sche = 10^(-c(40,30,20:10)/10) #13
  #sche = 10^(-c(100,90,80,70,60,50,40:20)/20) #27
  return(sche)
}

#####################################################################################################################################

ABC.SMC.Exact <- function(dp_info, hyper_params, params_bound){ 
  
  stepResultStore <- list()
  
  ### S1
  
  # outputs a vector
  epSchedule = schedule(dp_info = dp_info, hyper_params = hyper_params)
  
  thetaOld = rprior(hyper_params)
  weightOld = matrix(1/hyper_params$N, nrow = 1, ncol = hyper_params$N)
  attempt = matrix(1, nrow = 1, ncol = hyper_params$N)
  stepResult = rbind(thetaOld, weightOld, attempt)
  
  thetaNew = matrix(1, nrow = hyper_params$dim, ncol = hyper_params$N)
  weightNew = matrix(1, nrow = 1, ncol = hyper_params$N)

    t = 1
  stepResultStore[[t]] = stepResult
 
  repeat {
    thetaANDa <- foreach::foreach(i = 1:hyper_params$N) %dopar% {
      a <- 1
      repeat {
        newVec <- theta.update(t=t, thetaOld=thetaOld, weightOld=weightOld, hyper_params=hyper_params, params_bound=params_bound)
        modelData <- rmodel.par(newVec=newVec, dp_info=dp_info)
        logr <- threshold(privacyBudget=epSchedule[t], Data=modelData, dp_info=dp_info)
        u <- runif(1, min = 0, max = 1)
        if (logr > log(u)) {
          break
        }
        a <- a + 1
      }
      list(newVec = newVec, attempt = a)
    }  
    
    thetaNew = sapply(thetaANDa, function(x) x$newVec)
    weightNew = sapply(1:hyper_params$N, function(i){weight.update(t=t, i=i, thetaNew=thetaNew, thetaOld=thetaOld, weightOld=weightOld, hyper_params=hyper_params, params_bound=params_bound, dprior=dprior)})
    attempt = sapply(thetaANDa, function(x) x$attempt)
    
    thetaOld = thetaNew
    weightOld = weightNew / sum(weightNew)
    
    stepResult = rbind(thetaOld, weightOld, attempt)
    stepResultStore[[t+1]] = stepResult
    
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

res$time = unname(time$toc - time$tic)
res$sdp = dp_info$sdp

if (no_chains > 1) {
  save(res, file = paste("real_data_n", dim(vicMat)[1], "e", epsilon, "q", dp_info$q, "a", shape, "r", scale, "mo", hp$mu0, "so", sigma_square_1, "mt", hp$mu1, "st", sigma_square_2, "chain_id", chain_id, ".RData", sep = ""))
} else {
  save(res, file = paste("real_data_n", dim(vicMat)[1], "e", epsilon, "q", dp_info$q, "a", shape, "r", scale,  "mo", hp$mu0, "so", sigma_square_1, "mt", hp$mu1, "st", sigma_square_2, ".RData", sep = ""))
}
