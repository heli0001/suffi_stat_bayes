setwd("C:/Users/emmal/Desktop/bayesian sufficient")
rm(list=ls())
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

data_gene = function(dataid,N=40){
  set.seed(dataid)
  X = mvrnorm(N,mu=rep(0,9),Sigma=diag(9))
  beta_t = c(2,2,2,0.2,0.2,0.2,0,0,0)
  y = X%*%matrix(beta_t)+rnorm(N,1,1)
  data = as.data.frame(cbind(X,y))
  return (list(data=data,beta_t=beta_t))
}
simdata= mclapply(1,data_gene,N=40,mc.cores=1)

iter=5000
burn=1000

clu_func = function(dataid,simdata,N=40,cp=50/4,iter=100,burn=10){
  set.seed(dataid)
  data = simdata[[dataid]]$data
  ########### MCMC posterior samples
  Beta_pos = matrix(-99,nrow=iter,ncol=9)
  # cluster-variable indicator parameters
  Delta = matrix(-99,nrow=iter,ncol=9)
  #
  # probablity inclustion
  Prob_inc = matrix(-99,nrow=iter,ncol=9)
  
  ## starting value
  # propensity score
  beta_pos = beta_t
  delta = c(rep(1,6),0,0,0)
  prob_inc = rbeta(9,1,1)
  
  like_v = function(delta,i){
    if(all(delta == 0)){
      lval = 0
      return(lval)
    }
    ind_n0 = which(delta != 0)
    C = X[,ind_n0]
    phi_t = diag(1/cp,length(ind_n0),length(ind_n0))
    V_t = t(C)%*%C + solve(phi_t)
    mu = solve(V_t)%*%t(C)%*%y
    lval = -1/2*log((det(phi_t)))-1/2*log((det(V_t)))+1/2*t(mu)%*%V_t%*%mu
    return(lval)
  }
  
  for (k in 1:iter){
    ## sample beta
    ## check variable inclusion
    ind_0 = which(delta==0) # zero
    ind_n0 = which(delta!=0) # nonzero
    Xv = X[,ind_n0]
    beta_pos[ind_0] = 0
    sigma_zeta = diag(1/cp,length(ind_n0),length(ind_n0)) # specify 50/4 as variance
    beta_Sigma = solve(sigma_zeta+t(Xv)%*%Xv)
    beta_mu = beta_Sigma%*%t(Xv)%*%y
    beta_pos[ind_n0] = mvrnorm(1,beta_mu,beta_Sigma)
    
    ## Sample v ##
    for (i in 1:9){
      ## numerator part 
      delta[i]=0
      log_top = like_v(delta,i)
      
      ## denominator part 
      delta[i]=1
      log_bot = like_v(delta,i)
      
      ri = exp(log_top-log_bot)
      prob_deltai = 1/(1+(1-prob_inc[i])*ri/prob_inc[i])
      delta[i] = rbinom(1,1,prob_deltai)
    }
    
    ##  sample prob inclusion
    beta_a = 1 + delta
    beta_b = 1 + 1-delta
    prob_inc = rbeta(9,beta_a,beta_b)
    
    ## save
    Beta_pos[k,] = beta_pos
    Delta[k,] = delta
    Prob_inc[k,] = prob_inc
    
    print(k)
  }
  
  print(dataid)
  return(list(Gamma_pos=Gamma_pos,Var_gamma=Var_gamma,Zeta=Zeta,Delta=Delta,Prob_inc=Prob_inc))
  
}

apply(Prob_inc[(burn+1):iter,],2,mean)
apply(Beta_pos[(burn+1):iter,],2,mean)
apply(Delta[(burn+1):iter,],2,mean)
