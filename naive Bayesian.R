rm(list=ls())
library(MASS)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=11000
burn=1000
####################### naive method ###################
psfunc_naive = function(dataid,simudata,iter=100,burn=0){
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
  ####### initialize values
  beta_ps_pos[1,] = rep(1,2)
  var_beta[1] = 100

  for (g in 2:iter){
    data_ps = cbind(data_obs$x1,data_obs$x2) #x1,x2
    ## sample beta_ps
    wij = apply(data_ps%*%beta_ps_pos[g-1,],1,rpg,n=1,h=1)
    kij = data_obs$z-1/2
    L = kij/wij
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/var_beta[g-1]*diag(2))
    mu_ps = var_ps %*%(t(data_ps)%*%omega%*%L)
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    
    # 1-2 : var_beta
    var_beta[g] = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(beta_ps_pos[g,])%*%beta_ps_pos[g,]/2)
    
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  psfunc2 = function(x,para1){
    1/(1+exp(-x%*%para1))
  }
  data$w_naive = apply(data_ps,1,psfunc2,para1=beta_ps_pos_mean)
  effec_naive = mean(data$z*data$y_obs/data$w_naive-(1-data$z)*data$y_obs/(1-data$w_naive))
  print(dataid)
  time_span = proc.time()-ptm
  return(list(data=data,beta_ps_pos=beta_ps_pos,effec_naive=effec_naive,time_span=time_span))
}
naiveeffect = mclapply(1:20,psfunc_naive,simdata55m10_small,iter=iter,burn=burn,mc.cores = 1)

effect_naive = rep(-99,20)
effect_true = rep(-99,20)
for (i in 1:20){
  effect_naive[i] = naiveeffect[[i]]$effec_naive
  effect_true[i] = mean(simdata55m10_small[[i]]$data$y1-simdata55m10_small[[i]]$data$y0)
}
mean(effec_naive)
mean(effect_true)
mean(effec_naive)-mean(effect_true)
mean((effec_naive-effect_true)^2)
