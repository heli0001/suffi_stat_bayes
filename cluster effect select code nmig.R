setwd("C:/Users/emmal/Desktop/bayesian sufficient")
rm(list=ls())
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=100
burn=10

clu_func = function(dataid,simdata,iter=100,burn=10){
  ptm = proc.time()
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
  #
  # variance of spike and slab distribution
  phi1 = rep(-99,iter) 
  
  ## starting value
  # propensity score
  gamma_pos[1,] = c(1,1)
  var_gamma[1] = 100
  zeta[1,] = u*1
  delta[1,] = rep(1,m)
  prob_inc[1,] = rep(0.5,m)
  phi1[1] = 50
  
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
    rzeta = rep(1,m)
    rzeta[ind_0] = 1/0.00025 ## specify small r
    sigma_zeta = diag(rzeta*1/phi1[k-1],m,m) # specify phi as variance
    # update zeta
    data_ps = cbind(data_obs$x1,data_obs$x2) #x1,x2
    w_cl = apply(data_ps%*%gamma_pos[k,]+Istar%*%zeta[k-1,],1,rpg,n=1,h=1)
    k_cl = data_obs$z-1/2
    L_cl = k_cl/w_cl
    omega_cl = diag(w_cl)
    cl_Sigma = solve(sigma_zeta+t(Istar)%*%omega_cl%*%Istar)
    cl_mu = cl_Sigma%*%(t(Istar)%*%omega_cl%*%(L_cl-data_ps%*%gamma_pos[k,]))
    #zeta[k,] = rmvnorm(1,cl_mu,cl_Sigma,method = "svd")
    zeta[k,] = mvrnorm(1,cl_mu,cl_Sigma)
    
    # 1-4 delta
    pslab = dnorm(zeta[k,],0,sqrt(phi1[k-1]))*prob_inc[k-1,]
    pspike = dnorm(zeta[k,],0,sqrt(0.00025*phi1[k-1]))*(1-prob_inc[k-1,])
    delta[k,] = rbinom(m,1,pslab/(pspike+pslab))
    
    # step 4 : update hyperparameters
    # 4-1 prob inclusion
    beta_a = 1 + delta[k,]
    beta_b = 1 + 1-delta[k,]
    prob_inc[k,] = rbeta(m,beta_a,beta_b)
    # for (i in 1:m){
    #   beta_a = 1 + delta[k,i]
    #   beta_b = 1 + 1-delta[k,i]
    #   prob_inc[k,i] = rbeta(1,beta_a,beta_b)
    # }
    
    # 4-2 slab-spike variance
    b1_phi = zeta[k,]^2%*%rzeta
    phi1[k] = rinvgamma(1,shape=5+m/2,rate=50+b1_phi/2)
    # b1_phi = 0
    # for (i in 1:m){
    #   b1_phi = b1_phi+zeta[k,i]^2*rzeta[i]
    # }
    # phi1[k] = rinvgamma(1,shape=5+m/2,rate=50+b1_phi/2)
    
    print(k)
  }
  
  gamma_pos_mean = apply(gamma_pos[(burn+1):iter,],2,mean)
  zeta_mean = apply(zeta[(burn+1):iter,],2,mean)
  zetaij_mean = Istar%*%zeta_mean
  psfunc_est = function(x,para){
    1/(1+exp(-x[1:2]%*%para-x[3]))
  }
  data_ps = as.matrix(cbind(data_obs$x1,data_obs$x2,zetaij_mean))
  data_obs$w_nmig = apply(data_ps,1,psfunc_est,para=gamma_pos_mean)
  effec_nmig = mean(data_obs$z*data_obs$y_obs/data_obs$w_nmig-(1-data_obs$z)*data_obs$y_obs/(1-data_obs$w_nmig))
  
  print(dataid)
  time_span = proc.time() - ptm
  return(list(gamma_pos=gamma_pos,var_gamma=var_gamma,
              zeta=zeta,delta=delta,prob_inc=prob_inc,phi1=phi1,effec_nmig=effec_nmig,time_span=time_span))
}
clu_res = mclapply(1:2,clu_func,simdatam20,iter=iter,burn=burn,mc.cores = 1)
apply(prob_inc[(burn+1):iter,],2,mean)
apply(zeta[(burn+1):iter,],2,mean)
apply(gamma_pos[(burn+1):iter,],2,mean)
apply(delta[(burn+1):iter,],2,mean)

#### 