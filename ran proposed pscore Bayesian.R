rm(list=ls())
library(MASS)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000
####################### random method ###################
psfunc_ran = function(dataid,simudata,iter=100,burn=0){
  ptm = proc.time()
  set.seed(dataid)
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  N = simudata[[dataid]]$N
  m = simudata[[dataid]]$m
  n = simudata[[dataid]]$n
  u = simudata[[dataid]]$u
  suff_T = simudata[[dataid]]$suff_T
  
  ######## MCMC parameters
  beta_ps_pos = matrix(-99,nrow=iter,ncol=2) # x1,x2
  var_beta = rep(-99,iter)
  alpha_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha
  ####### initialize values
  beta_ps_pos[1,] = rep(1,2)
  var_beta[1] = 100
  alpha_ps_pos[1,] = u
  var_alpha_pos[1] = 1
  data_ps = as.matrix(data_obs[,3:4])
  ## dummy variable/indicator matrix 
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (g in 2:iter){
    ## sample alpha_ps
    data_ps = cbind(data_obs$x1,data_obs$x2) #x1,x2
    w_cl = apply(data_ps%*%beta_ps_pos[g-1,]+Istar%*%alpha_ps_pos[g-1,],1,rpg,n=1,h=1)
    k_cl = data_obs$z-1/2
    L_cl = k_cl/w_cl
    omega_cl = diag(w_cl)
    cl_Sigma = solve(1/var_alpha_pos[g-1]*diag(m)+t(Istar)%*%omega_cl%*%Istar)
    cl_mu = cl_Sigma%*%(t(Istar)%*%omega_cl%*%(L_cl-data_ps%*%beta_ps_pos[g-1,]))
    alpha_ps_pos[g,]  = mvrnorm(1,cl_mu,cl_Sigma)
    
    ## sample var_alpha
    ralpha = sum((alpha_ps_pos[g,])^2)
    shape1 = (0.002+m)/2
    scale1 = (0.002*1 + ralpha)/2
    var_alpha_pos[g] = rinvgamma(1,shape=shape1,rate=scale1)
    
    ## sample beta_ps
    wij = apply(data_ps%*%beta_ps_pos[g-1,]+Istar%*%alpha_ps_pos[g,],1,rpg,n=1,h=1)
    kij = data_obs$z-1/2
    L = kij/wij
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/var_beta[g-1]*diag(2))
    mu_ps = var_ps %*%(t(data_ps)%*%omega%*%(L-Istar%*%alpha_ps_pos[g,]))
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    
    # 1-2 : var_beta
    var_beta[g] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta_ps_pos[g,])%*%beta_ps_pos[g,]/2)
    
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  alpha_ps_pos_mean = apply(alpha_ps_pos[(burn+1):iter,],2,mean)
  alphaij_mean = Istar%*%alpha_ps_pos_mean
  psfunc2 = function(x,para1){
    1/(1+exp(-x[1:2]%*%para1-x[3]))
  }
  data_ps2 = as.matrix(cbind(data_ps,alphaij_mean))
  data$w_ran = apply(data_ps2,1,psfunc2,para1=beta_ps_pos_mean)
  effec_ran = mean(data$z*data$y_obs/data$w_ran-(1-data$z)*data$y_obs/(1-data$w_ran))
  print(dataid)
  time_span = proc.time()-ptm
  return(list(data=data,beta_ps_pos=beta_ps_pos,alpha_ps_pos=alpha_ps_pos,var_alpha_pos=var_alpha_pos,effec_ran=effec_ran,time_span=time_span))
}
raneffect = mclapply(1:20,psfunc_ran,simdatam20,iter=iter,burn=burn,mc.cores = 1)
