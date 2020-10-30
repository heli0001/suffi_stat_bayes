setwd("C:/Users/emmal/Desktop/bayesian sufficient")
rm(list=ls())
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000

clu_func = function(dataid,simdata,m0=30,iter=100,burn=10){
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  u = simdata[[dataid]]$u
  
  ########### MCMC posterior samples
  #
  # propensity score
  gamma_pos = matrix(-99,nrow=iter,ncol=2)
  var_gamma = rep(-99,iter)
  #
  # hyper cluster parameters
  zeta = matrix(-99,nrow=iter,ncol=m)
  #
  # cluster-variable indicator parameters
  delta = matrix(-99,nrow=iter,ncol=m)
  #
  # probablity inclustion
  prob_inc = matrix(-99,nrow=iter,ncol=m)
  
  ## starting value
  # propensity score
  gamma_pos[1,] = c(1,1)
  var_gamma[1] = 100
  zeta[1,] = u*1
  delta[1,] = c(rep(0,m0),rep(1,(m-m0)))
  prob_inc[1,] = rbeta(m,1,1)
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (k in 2:iter){
    # Step 1: update propensity score parameters
    # 1-1 : gamma
    data_gamma = cbind(data_obs$x1,data_obs$x2)
    w_gamma = apply(data_gamma%*%gamma_pos[k-1,]+Istar%*%zeta[k-1,],1,rpg,n=1,h=1)
    k_gamma = data_obs$z-1/2
    L_gamma = k_gamma/w_gamma
    omega_gamma = diag(w_gamma)
    gamma_Sigma = solve(t(data_gamma)%*%omega_gamma%*%data_gamma + 1/var_gamma[k-1]*diag(2))
    gamma_mu = gamma_Sigma%*%(t(data_gamma)%*%omega_gamma%*%(L_gamma-Istar%*%zeta[k-1,]))
    gamma_pos[k,] = mvrnorm(1,gamma_mu,gamma_Sigma)
    
    # 1-2 : var_gamma
    var_gamma[k] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(gamma_pos[k,])%*%gamma_pos[k,]/2)
    
    # 1-3 : zeta, 1,2,...m
    ## check variable inclusion
    ind_0 = which(delta[k-1,]==0) # zero
    ind_n0 = which(delta[k-1,]!=0) # nonzero
    C = Istar[,ind_n0]
    sigma_zeta = diag(1/50,length(ind_n0),length(ind_n0)) # specify 50 as variance
    # update zeta
    data_ps = cbind(data_obs$x1,data_obs$x2) #x1,x2
    w_cl = apply(data_ps%*%gamma_pos[k-1,]+Istar[,ind_n0]%*%zeta[k-1,ind_n0],1,rpg,n=1,h=1)
    k_cl = data_obs$z-1/2
    L_cl = k_cl/w_cl
    omega_cl = diag(w_cl)
    cl_Sigma = solve(sigma_zeta+t(C)%*%omega_cl%*%C)
    cl_mu = cl_Sigma%*%(t(C)%*%omega_cl%*%(L_cl-data_ps%*%gamma_pos[k-1,]))
    zeta[k,ind_n0] = mvrnorm(1,cl_mu,cl_Sigma)
    zeta[k,ind_0] = 0
    
    # 1-4 delta
    for (i in 1:m){
      ## numerator part 
      delta_ti = rep(-99,m)
      delta_ti[i] = 0
      delta_ti[-i] = delta[k-1,-i]
      ind_n0 = which(delta_ti!=0)
      C = Istar[,ind_n0]
      phi_t = diag(1/50,length(ind_n0),length(ind_n0))
      V_t = t(C)%*%omega_cl%*%C + solve(phi_t)
      mu_top = solve(V_t)%*%t(C)%*%omega_cl%*%(L_cl-data_ps%*%gamma_pos[k-1,])
      mut_inv_omega_t_mu_top = t(mu_top)%*%V_t%*%mu_top
      log_top = -1/2*log((det(phi_t)))-1/2*log((det(V_t)))+1/2*mut_inv_omega_t_mu_top
      
      ## denominator part 
      delta_bi = rep(-99,m)
      delta_bi[i] = 1
      delta_bi[-i] = delta[k-1,-i]
      ind_n0 = which(delta_bi!=0)
      C = Istar[,ind_n0]
      phi_t = diag(1/50,length(ind_n0),length(ind_n0))
      V_t = t(C)%*%omega_cl%*%C + solve(phi_t)
      mu_bot = solve(V_t)%*%t(C)%*%omega_cl%*%(L_cl-data_ps%*%gamma_pos[k-1,])
      mut_inv_omega_t_mu_bot = t(mu_bot)%*%V_t%*%mu_bot
      log_bot = -1/2*log((det(phi_t)))-1/2*log((det(V_t)))+1/2*mut_inv_omega_t_mu_bot
      
      ri = exp(log_top-log_bot)
      prob_deltai = 1/(1+(1-prob_inc[k-1,i])*ri/prob_inc[k-1,i])
      delta[k,i] = rbinom(1,1,prob_deltai)
      delta[k-1,i] = delta[k,i] # use the updated one for the following deltas
    }
    
    # step 4 : update hyperparameters
    # 4-1 prob inclusion
    for (i in 1:m){
      beta_a = 1 + delta[k,i]
      beta_b = 1 + 1-delta[k,i]
      prob_inc[k,i] = rbeta(1,beta_a,beta_b)
    }
    
    print(k)
  }
  print(dataid)
  return(list(gamma_pos=gamma_pos,var_gamma=var_gamma,
              zeta=zeta,delta=delta,prob_inc=prob_inc))
}
res = mclapply(1,clu_func,simdata00m50_small,m0=30,iter=100,burn=10,mc.cores = 1)
apply(prob_inc[(burn+1):iter,],2,mean)
apply(zeta[(burn+1):iter,],2,mean)
apply(gamma_pos[(burn+1):iter,],2,mean)
