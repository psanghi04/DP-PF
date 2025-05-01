
###   the model is as follows:
###   x | mu, phi ~ N_p(mu, phi^{-1})
###   y | x, beta, tau ~ N(x*beta, tau^{-1} I)
###   with priors 
###   mu ~ N(theta, Sigma)
###   Phi ~ Wishart_p(d,W)
###   beta | tau ~ N_{p+1}(m,tau^{-1}V^{-1})
###   tau ~ Gamma(a/2, b/2), using the shape, rate parameterization
###   In total, the parameters are (beta, tau, mu, Phi)
###   the hyperparameters are (m, V, a, b, theta, Sigma, d, W)
###   and sufficient statistics are (XtX, XtY, YtY, X_bar, and S=sum((X-X_bar)%*%top(X-X_bar)) )

library(VGAM)###for laplace distribution
library(MASS)### for multivariate normal simulation
library(tictoc)
library(emdbook)
library(LaplacesDemon)
library(truncnorm)
library(optparse)


## clamping + location-scale transformation
## use this function to normalize each covariate x[-,j] and y
normalize <- function(vec, left, right){
  return(2 * (pmin(right, pmax(left, vec)) - left) / (right - left) - 1 );
}

## prior on N
g_log <- function(N){
  return(0)## constant improper prior.
  #return(dpois(N,lambda=100,log=TRUE))
}

#log-probability of transitioning from N to N_.
q_log <-function(N_,N){
  if(N==1 & N_==2){
    return(0)# log(1)
  }
  if(N>1 & abs(N-N_)==1){
    return(log(1/2))
  }
  else{
    print("q is zero")
    return(-Inf)
  }
}

## initialize the Gibbs sampler----
rinit_data <- function(N,params){
  
  ### generate X from model
  x = mvrnorm(n=N,mu=params$mu,Sigma = params$phi_inv)
  X = cbind(1, x);
  
  ### generate Y given Xbeta + noise
  Y <- as.numeric(X %*% params$beta + rnorm(N, 0, s=sqrt(1/params$tau)));
  ### clamp and get diffT1, diffT2, diffT3
  return(data = list(N=N, X = X, Y = Y));
  #keep track of current N as well as current upper bound on N (N_upper)
}

## initialize the Gibbs sampler----
rinit_gibbs <- function(N_guess,hyperparams = hyperparams,params = NULL){
  #N is an initial guess at N. Could be N_DP.
  
  N_upper = 2*N_guess # a crude upper bound for the amount of space we will need.
  ## generate parameters from prior
  if(is.null(params)){
    params = list()
    params$tau = hyperparams$a/hyperparams$b#rgamma(n=1,shape = hyperparams$a/2,rate = hyperparams$b/2)
    params$beta <- hyperparams$m#mvrnorm(n=1,mu=hyperparams$m,Sigma = hyperparams$V_inv/params$tau)
    params$mu <- hyperparams$theta#mvrnorm(n=1,mu=hyperparams$theta,Sigma = ginv(hyperparams$Sigma))
    
    # params$phi <- ginv(hyperparams$W_inv)#rWishart(n=1, df=hyperparams$d, Sigma=ginv(hyperparams$W_inv))[,,1]
    # params$phi_inv <- ginv(params$phi)
    ###
    params$phi <- hyperparams$d
    params$phi_inv <- 1/params$phi
  }
  #beta_mat <- matrix(rnorm((p+1)*N_upper,0,tau_star),nrow=(p+1),ncol=N_upper)
  ### generate X from model
  state = rinit_data(N=N_upper,params)
  state$N = N_guess
  state$N_upper = N_upper
  state$params = params
  state$accMean_change = 0
  state$accMean_jump=0
  return(state);
  #keep track of current N as well as current upper bound on N (N_upper)
}

## compute distance from summary statistics to sdp--------
## only used the first time.
get_diffT <- function(state, clbounds, sdp){
  ### clamp X and Y
  X_cl <- apply(state$X[1:state$N,], 2, function(x) normalize(x, clbounds[1], clbounds[3]));
  X_cl[, 1] <- 1;
  Y_cl <- sapply(state$Y[1:state$N], function(y) normalize(y, clbounds[2], clbounds[4]));
  ### compute summary statistics and distance to noisy sdp
  if(is.null(sdp)){
    state$T1 <- t(X_cl) %*% Y_cl
    state$T2 <- t(Y_cl) %*% Y_cl
    state$T3 <- t(X_cl) %*% X_cl
  }
  else{
    state$T1 <- t(X_cl) %*% Y_cl - sdp[[1]];
    state$T2 <- t(Y_cl) %*% Y_cl - sdp[[2]];
    state$T3 <- t(X_cl) %*% X_cl - sdp[[3]];
  }
  return(state);
}


### these are not used in the DP density, just for the non-private posterior.
###  just run this after the database is updated.
get_suffStats = function(state){
  state$XtX <- t(state$X[1:state$N,])%*%state$X[1:state$N,]
  state$XtY <- t(state$X[1:state$N,])%*%state$Y[1:state$N]
  state$YtY <- t(state$Y[1:state$N])%*%state$Y[1:state$N]
  state$Xbar <- state$XtX[1,2:(p+1)]/state$N
  ###
  state$S <- (state$N-1)*var(state$X[1:state$N,2])
  return(state)
}

## compute density of privacy noise---------
## using Laplace densities
## operations on log-scale for stability
get_logeta <- function(state, N_DP,privacy_params){
  noise_scale = privacy_params$deltaa/ privacy_params$eps
  return(sum(VGAM::dlaplace(state$T1 , 0, noise_scale, TRUE)) + 
           VGAM::dlaplace(state$T2, 0, noise_scale, TRUE) + 
           sum(VGAM::dlaplace(state$T3[lower.tri(state$T3)], 0, noise_scale, TRUE)) + 
           sum(VGAM::dlaplace(diag(state$T3)[2:length(state$T2)], 0, noise_scale, TRUE))+
           VGAM::dlaplace(state$N-N_DP,0,1/privacy_params$ep_N,TRUE));
}


## calculate individual contribution to t(x,y)
get_diffT_individual <- function(x, y, clbounds){
  xcl <- normalize(x, clbounds[1], clbounds[3]);
  xcl[1] <- 1;
  ycl <- normalize(y, clbounds[2], clbounds[4]);
  diffti = list(diffti1 = xcl %o% ycl, 
                diffti2 = t(ycl)%*%ycl, 
                diffti3 = xcl %o% xcl)
  return(diffti);
}


## difference in each individual contribution to t(x,y)
## between current state and proposed state
## ti(xi_star,yi_star) - to(xi, yi)
get_diffeacht <- function(x, x_, y, y_, clbounds){
  xcl <- normalize(x, clbounds[1], clbounds[3]);
  xcl[1] <- 1;
  ycl <- normalize(y, clbounds[2], clbounds[4]);
  xcl_ <- normalize(x_, clbounds[1], clbounds[3]);
  xcl_[1] <- 1;
  ycl_ <- normalize(y_, clbounds[2], clbounds[4]);
  diffti = list(diffti1 = xcl_ %o% ycl_ - xcl %o% ycl, 
                diffti2 = t(ycl_)%*%ycl_ - t(ycl)%*%ycl, 
                diffti3 = xcl_ %o% xcl_ - xcl %o% xcl)
  return(diffti);
}


## independent-MH kernel for xi, yi given xnoti, ynoti, beta, sdp, and other parameters
## propose (xi_star,yi_star) pair from the model f(x,y | theta);
## acceptance ratio is eta(tstar) / eta(t);
update_XiYi_gibbs_change <- function(i,state,hyperparams=hyperparams,privacy_params,N_DP){
  acc <- 0;
  ## propose xi and yi from model
  #print(state$params)
  xi <- c(1,as.numeric(mvrnorm(n=1,mu=state$params$mu,Sigma=state$params$phi_inv)))
  #print(xi)
  yi <- xi %*% state$params$beta + sqrt(1/state$params$tau) * rnorm(1);
  ## compute T*
  diffti <- get_diffeacht(state$X[i,],xi, state$Y[i], yi, privacy_params$clbounds);
  
  state_proposed = state
  state_proposed$X[i,] <- xi;
  state_proposed$Y[i] <- yi;
  state_proposed$T1=state$T1 + diffti$diffti1; ## xty
  state_proposed$T2=state$T2 + diffti$diffti2; ##yty
  state_proposed$T3=state$T3 + diffti$diffti3; ## xtx
  state_proposed$logeta = get_logeta(state_proposed,N_DP,privacy_params)
  #print(state_proposed$logeta)
  
  ## accept and reject step 
  logu <- log(runif(1));
  if(logu < state_proposed$logeta - state$logeta){
    acc <- 1;
    ## accept and change all associated values of x and y 
    state = state_proposed
  }
  state$acc <- acc;
  return(state)
}

## independent-MH kernel for xi, yi given xnoti, ynoti, beta, sdp, and other parameters
## propose (xi_star,yi_star) pair from the model f(x,y | theta);
## acceptance ratio is eta(tstar) / eta(t);
update_XiYi_gibbs_add <- function(state,hyperparams=hyperparams,privacy_params,N_DP){
  acc <- 0;
  ## propose xi and yi from model
  xi <- c(1,as.numeric(mvrnorm(n=1,mu=state$params$mu,Sigma=state$params$phi_inv)))
  yi <- xi %*% state$params$beta + sqrt(1/state$params$tau) * rnorm(1);
  ## compute T*
  diffti <- get_diffT_individual(xi, yi, privacy_params$clbounds);
  if(state$N+1>state$N_upper){
    print("augmented data structure doubled")
    state$X <-rbind(state$X,state$X)
    state$Y <-c(state$Y,state$Y)
    state$N_upper=2*state$N_upper
  }
  state_proposed = state
  state_proposed$X[state$N+1,] <- xi;
  state_proposed$Y[state$N+1] <- yi;
  state_proposed$T1=state$T1 + diffti$diffti1; ## xty
  state_proposed$T2=state$T2 + diffti$diffti2; ##yty
  state_proposed$T3=state$T3 + diffti$diffti3; ## xtx
  state_proposed$N = state$N+1
  state_proposed$logeta = get_logeta(state_proposed,N_DP,privacy_params)
  
  ## accept and reject step 
  logu <- log(runif(1));
  log_accept = state_proposed$logeta - state$logeta+
    q_log(state$N,state$N+1)-q_log(state$N+1,state$N)+
    g_log(state$N+1)-g_log(state$N)
  
  if(logu < log_accept){
    acc <- 1;
    ## accept and change all associated values of x and y 
    state = state_proposed
  }
  state$acc <- acc;
  return(state)
}

update_XiYi_gibbs_delete <- function(i,state,privacy_params,N_DP){
  acc <- 0;
  
  ## compute T*
  diffti <- get_diffT_individual(state$X[i,], state$Y[i], privacy_params$clbounds);
  state_proposed = state
  state_proposed$X[i,] <- state$X[state$N,];
  state_proposed$Y[i] <- state$Y[state$N] ;
  state_proposed$T1=state$T1 - diffti$diffti1; ## xty
  state_proposed$T2=state$T2 - diffti$diffti2; ##yty
  state_proposed$T3=state$T3 - diffti$diffti3; ## xtx
  state_proposed$N = state$N-1
  state_proposed$logeta = get_logeta(state_proposed,N_DP,privacy_params)
  
  ## accept and reject step 
  logu <- log(runif(1));
  log_accept = state_proposed$logeta - state$logeta+
    q_log(state$N,state$N-1)-q_log(state$N-1,state$N)+
    g_log(state$N-1)-g_log(state$N);
  if(logu < log_accept){
    acc <- 1;
    ## accept and change all associated values of x and y 
    state = state_proposed
  }
  state$acc <- acc;
  return(state)
}


## update params given (x,y) using conjugate posterior
# update_params_gibbs <- function(state,hyperparams=hyperparams){
#   ### get sufficient statistics
#   state = get_suffStats(state)
#   m_star = ginv(state$XtX+hyperparams$V_inv)%*% (state$XtY + hyperparams$V_inv%*% hyperparams$m)
#   V_inv_star = state$XtX + hyperparams$V_inv
#   V_star = ginv(V_inv_star)
#   a_star = hyperparams$a + state$N
#   b_star = hyperparams$b + state$YtY + t(hyperparams$m)%*%hyperparams$V_inv%*% hyperparams$m-t(m_star)%*% V_inv_star%*% m_star
#   
#   
#   tau_new = rgamma(n=1,shape = a_star/2,rate = b_star/2)
#   beta_new = mvrnorm(n=1,mu=m_star,Sigma = V_star/tau_new)
#   
#   Sigma_star = ginv(hyperparams$Sigma_inv + state$N*state$params$phi)
#   theta_star = Sigma_star %*%(hyperparams$Sigma_inv %*% hyperparams$theta + state$N*state$params$phi%*%state$Xbar)
#   
#   mu_new = mvrnorm(n=1,mu=theta_star,Sigma = Sigma_star)
#   
#   d_star = hyperparams$d + state$N
#   W_star = ginv(hyperparams$W_inv + state$S + state$N*(mu_new - state$Xbar)%*%t(mu_new - state$Xbar))### threw in inverse?
#   
#   phi_new = rWishart(n=1,df=d_star,Sigma = W_star)[,,1]
#   
#   
#   params_new = list(tau=tau_new,beta=beta_new,mu=mu_new,phi=phi_new,phi_inv = ginv(phi_new))
#   state$params = params_new
#   return(state);
# }

################################################################################################################################################


metrop1 = function(logD,init,nbatch,scale,blen){
  ###   blen doesnt do anything.
  dim = length(init)
  out = list(accept = 0,
             batch = matrix(rep(0,nbatch*dim),nrow=nbatch,ncol=dim))
  U = matrix(runif(nbatch*dim,min=0,max=1),nrow=nbatch,ncol=dim)
  Prop = matrix(rnorm(nbatch*dim,m=0,s=scale),nrow=nbatch,ncol=dim)
  
  for(r in 1:nbatch){
    if(r==1){
      out$batch[1,] =init
      oldLogD = logD(init)
    }
    else
      out$batch[r,] = out$batch[r-1,]
    ###
    
    for(i in 1:dim){
      oldVector = out$batch[r,]
      
      newVector = oldVector
      newVector[i] = newVector[i] + Prop[r,i]
      ###
      if(newVector[2] <= 0 | newVector[3] <= 0){newVector = oldVector}
      ###
      newLogD = logD(newVector)
      
      TestValue = exp((newLogD - oldLogD))
      if(U[r,i]<TestValue){
        out$batch[r,] = newVector
        out$accept = out$accept + 1/(nbatch*dim)
        oldLogD = newLogD
      }
      #print(c("i is", i)) 
    }
    #print(c("r is", r)) 
  }
  return(out)
}

getScale = function(logA){
  scale = .1
  prevBelow=0
  prevAbove = Inf
  count=0
  init = c(0,1,1,0,0) ####
  #init = beta_hat
  
  out = metrop1(logA,init=init,nbatch=200,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
  #out$accept
  # out = metrop(logD,init=t(tail(out$batch,n=1)),nbatch=1000,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
  while((out$accept<.1 | out$accept>.2)&(count<=10)){
    if(out$accept<.1){
      prevAbove = scale 
    }
    else if(out$accept>.2){
      prevBelow = scale 
    }
    
    if(prevAbove<Inf){
      scale = (1/2)*(prevBelow+prevAbove) 
    }
    else{
      scale = 2*scale 
    }
    
    out = metrop1(logA,init=init,nbatch=200,scale=scale,blen=1)# accept around .2 ish # used .15 now .029
    count = count+1
  }
  if(count==11){
    print('scale did not converge')
    return(0)
  }
  return(scale)
}

update_params_gibbs <- function(state,hyperparams=hyperparams){
  ### get sufficient statistics
  state = get_suffStats(state)
  ###
  stateParams = c(state$params$mu, state$params$phi, state$params$tau, state$params$beta)
  ###
  
  logD = function(theta){
    ###
    mu = theta[1]
    phi = theta[2]
    tau = theta[3]
    beta = theta[4:5]
    
    # Log Improper Prior 
    logPrior = tryCatch(
      {
        ### Improper uniform prior on mu and beta
        # result = dweibull(x = tau, scale = 1.25, shape = 2, log = T) + dwishart(Omega = PhiMatrix, nu = 2, S = diag(2), log = T)
        ## t prior on mu and beta
        result = dweibull(x = tau, scale = 1.25, shape = 2, log = T) + log(2) + dt(x = phi, df = 2, log = T) +
          dt(x = mu, df = 2, log = T) + sum(dt(x = beta, df = 2, log = T))
      },
      error = function(e) {
        message('An Error Occurred')
        print(e)
        return(-Inf)
      }, warning = function(w){
        message('A Warning Occurred')
        print(w)
        return(-Inf)
      }, finally = {
        result
      }
    )
    
    
    # Log data likelihood
    logData = 0 
    for(n in 1:state$N){
      xn = state$X[n,]
      yn = state$Y[n]
      logData = tryCatch(
        {
          result = logData + dnorm(x = xn[-1], mean = mu, sd = 1/sqrt(phi), log = T) + dnorm(x = yn, mean = sum(xn*beta), sd = 1/sqrt(tau), log = T) #could use sufficient statistics
        },
        error = function(e) {
          message('An Error Occurred')
          print(e)
          return(-Inf)
        }, warning = function(w){
          message('A Warning Occurred')
          print(w)
          return(-Inf)
        }, finally = {
          result
        }
      )
    }
    
    logPost = logPrior + logData
    return(logPost)
  }
  
  ###
  while (state$scale==0) {
    scale=getScale(logA=logD)
    print(c("getScale completed with scale =",scale))
    state$scale=scale
  }
  ###
  
  nbatch = 10 ### to be tuned
  out = metrop1(logD=logD, init = stateParams, nbatch = nbatch, scale = state$scale, blen = 1)
  params_new = out$batch[nbatch-1,]
  
  mu_new=params_new[1]
  phi_new=params_new[2]
  tau_new=params_new[3]
  beta_new=params_new[4:5]
  params_new = list(tau=tau_new,beta=beta_new,mu=mu_new,phi=phi_new,phi_inv = ginv(phi_new))
  
  state$params = params_new
  return(state);
}


################################################################################################################################################


## One iteration of the Gibbs sampler for linear regression
# swaps decides how many individuals should be swapped
# fullUpdate does a complete update of the database. ignores swaps in this case.
# jumps decides how many individuals should be added/deleted in one "iteration". 
# intuitively, we want loops large enough that we are likely to have made at least one change to the database.
onestep_gibbs <- function(state,hyperparams,privacy_params,N_DP,swaps = 1,jumps=1, fullUpdate = TRUE){
  # if(is.nan(privacy_params$ep_N)==TRUE | privacy_params$ep_N==Inf)
  #   jumps=0
  ### step1
  acc_change <- 0;
  acc_jump <- 0;
  ###  update the database by swapping rows
  if(fullUpdate == TRUE){
    ###  full sweep, like in the original paper
    for(i in 1:state$N){
      state <- update_XiYi_gibbs_change(i,state,hyperparams,privacy_params,N_DP);
      acc_change <- state$acc + acc_change;
    }
    state$accMean_change <- acc_change/state$N;
  }
  else if(swaps>0){
    ### randomly update "swaps" many rows
    for(i in 1:swaps){
      i = sample(size = 1,1:state$N)
      state <- update_XiYi_gibbs_change(i,state,hyperparams,privacy_params,N_DP);
      acc_change <- state$acc + acc_change;
    }
    state$accMean_change <- acc_change/swaps;
  }
  
  ### reversible jump steps. Do "jumps" many times
  # if(jumps>0){
  #   for(index in 1:jumps){
  #     if(state$N==1){
  #       state <- update_XiYi_gibbs_add(state,hyperparams,privacy_params,N_DP);
  #     }
  #     else{
  #       coin <- rbinom(n=1,size=1,prob=1/2)
  #       if(coin==1){
  #         state <- update_XiYi_gibbs_add(state,hyperparams,privacy_params,N_DP);
  #         #print(state$N)
  #       }
  #       else if(coin==0){
  #         i = sample(size=1,1:state$N);
  #         state <- update_XiYi_gibbs_delete(i,state,privacy_params,N_DP);
  #         #print(state$N)
  #       }
  #       else{
  #         print("coin is broken")
  #       }
  #       
  #     }
  #     acc_jump <- state$acc + acc_jump;
  #     
  #   }
  # }
  #state$accMean_jump <- acc_jump / jumps;
  state = get_suffStats(state)
  #print(state$XtY)
  ### better to update the parameters after updating the database.
  state <- update_params_gibbs(state,hyperparams);### step2 of gibbs sampler
  #state$accMean <- (acc_change + acc_jump)
  #print(acc/loops)
  #print(state$params)
  return(state);
}

run_chain = function(sdp, iteration=1,num_cycles=15000,N=1000,ep_N=NaN,eps,N_guess, ####
                     hyperparams=hyperparams,privacy_params = privacy_params,
                     swaps=1,jumps=1,fullUpdate=TRUE){
  
  
  ### true data parameters
  params_true = list(beta = c(0,2),
                     tau = 1,
                     mu = 1,
                     phi = 1,
                     phi_inv = 1)
  
  
  #set.seed(1)### fix the seed for the true (sensitive) data. 
  #data = rinit_data(N,params_true)
  #set.seed(iteration)### use the iteration as the seed for the privatized statistics
  
  
  N_DP =  N + rlaplace(n=1,location=0,scale=1/privacy_params$ep_N)
  
  #data = get_diffT(data,privacy_params$clbounds,sdp=NULL)
  data <- sdp
  # sample the noises
  noise_scale = privacy_params$deltaa/privacy_params$eps
  l1 <- matrix(VGAM::rlaplace(p + 1, 0, noise_scale), ncol = 1);
  l2 <- VGAM::rlaplace(1,0,noise_scale);
  l3 <- matrix(VGAM::rlaplace( ( p + 1) * (p + 1), 0,  noise_scale), ncol = (p + 1));
  l3[lower.tri(l3)] <- t(l3)[lower.tri(l3)];
  l3[1,1] <- 0;
  sdp = list(data$T1+l1,data$T2+l2,data$T3+l3) ### observed DP output, along with N_DP.
  
  #get_logeta(data,sdp)
  
  # need to guess what N should be. Since N_DP is unbiased, we can use it. 
  if(is.na(N_guess))
    N_guess = max(2,ceiling(N_DP))
  ## initiate the chain 
  state <-rinit_gibbs(N_guess,hyperparams)
  state <- get_diffT(state, privacy_params$clbounds, sdp);
  state$logeta <- get_logeta(state,N_DP,privacy_params);
  ###
  #state$scale=0
  state$scale=1
  ###
  
  ## store parameter values, acceptance rate, and log-posterior
  logeta_chain <- double(num_cycles);
  beta_chain <- matrix(NA, nrow = p + 1, ncol = num_cycles);
  tau_chain <- double(num_cycles);
  ###
  mu_chain <- double(num_cycles);
  phi_chain <- double(num_cycles);
  ###
  n_chain <- rep(NA,num_cycles);
  acc_chain_change <- double(num_cycles);
  acc_chain_jump <- double(num_cycles);
  #acc_chain_mean <- double(num_cycles);
  
  
  for(iter in 1 : num_cycles){
    state <- onestep_gibbs(state,hyperparams,privacy_params,N_DP,jumps=1, fullUpdate = TRUE)
    logeta_chain[iter] <- state$logeta;
    
    beta_chain[,iter] <- state$params$beta;
    tau_chain[iter] <- state$params$tau;
    ###
    mu_chain[iter] <- state$params$mu;
    phi_chain[iter] <- as.numeric(state$params$phi)
    ###
    n_chain[iter] <- state$N;
    acc_chain_jump[iter] <- state$accMean_jump;
    acc_chain_change[iter] <- state$accMean_change;
    if(iter %% 100 == 0) cat("iteration", iter, "(",100*iter/num_cycles,"% )\n");
  }
  
  outputs = list(logeta_chain=logeta_chain,
                 tau_chain = tau_chain,
                 beta_chain=beta_chain,
                 mu_chain = mu_chain,
                 phi_chain = phi_chain,
                 n_chain=n_chain,
                 acc_chain_jump=acc_chain_jump,
                 acc_chain_change=acc_chain_change,
                 sdp=sdp,
                 N_DP=N_DP)
  return(outputs)
}

# Define options for command line arguments
option_list <- list(
  make_option(c("-n", "--sample_size"), type="integer", default=200, 
              help="Sample size", metavar="number"),
  make_option(c("-e", "--epsilon"), type="numeric", default=1, 
              help="Privacy budget", metavar="number"),
  make_option(c("-d", "--seed"), type="numeric", default=731, 
              help="Seed for reproducibility of data", metavar="number"),
  make_option(c("-i", "--iter"), type="numeric", default=731, 
              help="Run Number", metavar="number")
)

# Parse options
args <- parse_args(OptionParser(option_list=option_list))

###### Simulation ######
num_cycles=10000
ep_N=10^7 

### Change for different settings
N=args$sample_size #sample size: n
N_guess = N
p=1 #no. of predictors
eps=args$epsilon #epsilon
iter <- args$iter
seed <- args$seed

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
 
 return( c(ur_s1,ur_s2,ur_s3) + rlaplace(((p+2)^2+p)/2, location = 0, scale = deltaa/eps) )
}

# Random sdp
mu = 1
Phi = 1
beta = c(0, 2)
tau = 1
p = 1
deltaa = p^2 + 4*p + 3

set.seed(1)
missingX <- mvrnorm(N, mu = mu, Sigma = 1/Phi) #model matrix generation~N(mu, 1) ##
missingY <- cbind(1, missingX) %*% beta + rnorm(N, sd = sqrt(1/tau)) #model error~N(0, 1)

sdp <- privatized.stat(cbind(missingY, missingX), L=-5, U=5) 

# Fixed sdp
# sdp <- c(0, 0, N/10, 0, N/10)

#set.seed(seed)

# Create the matrices
matrix1 <- matrix(sdp[1:2], ncol = 1)
matrix2 <- matrix(sdp[3], ncol = 1)
matrix3 <- matrix(c(N, sdp[4], sdp[4], sdp[5]), ncol = 2)

# Combine the matrices into a list
sdp <- list("T1" = matrix1, "T2" = matrix2, "T3" = matrix3)

hyperparams = list(p=1,
                   m = rep(0,p+1),
                   V_inv = diag(p+1),
                   a = 2,
                   b = 2,
                   theta = 0,
                   Sigma_inv = 1,
                   d=2)# df needs to be >=p

privacy_params = list(eps = eps,
                      ep_N = ep_N,
                      deltaa = p^2+4*p+3,  ### bounded should have been p^2+4*p+3, not p^2+3*p+3...; unbdd: p^2+4*p+3
                      clbounds=c(-5,-5,5,5))

tic()
chain = run_chain(sdp, iteration=1,num_cycles=num_cycles,N=N,N_guess=NaN,
                  hyperparams=hyperparams,
                  privacy_params = privacy_params,
                  swaps=1,jumps=1,fullUpdate=TRUE)
time <- toc()

chain$time <- unname(time$toc - time$tic)
chain$sdp <- sdp

print("Saving results")
save(chain, file = paste("linear_regression/mcmc/non_conjugate_prior/data/mcmc_nonconjugate_n", N, "e", eps, "seed", seed, "i", iter, ".RData", sep = ""))
print("Saved Results")