load("simdata00m50_small.RData")
load("simdata05m50_small.RData")
load("simdata50m50_small.RData")
load("simdata55m50_small.RData")
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=5000
burn=1000

clu_func = function(dataid,simdata,cp=50/4,iter=iter,burn=burn){
  ptm = proc.time()
  set.seed(dataid)
  
  data = simdata[[dataid]]$data
  data_obs = simdata[[dataid]]$data_obs
  N = dim(data_obs)[1]
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  u = simdata[[dataid]]$u
  
  ########### MCMC posterior samples
  #
  # propensity score
  Gamma_pos = matrix(-99,nrow=iter,ncol=2)
  Var_gamma = rep(-99,iter)
  #
  # hyper cluster parameters
  Zeta = matrix(-99,nrow=iter,ncol=m)
  #
  # cluster-variable indicator parameters
  Delta = matrix(-99,nrow=iter,ncol=m)
  #
  # probablity inclustion
  Prob_inc = matrix(-99,nrow=iter,ncol=m)
  
  ## starting value
  # propensity score
  gamma_pos = c(1,1)
  var_gamma = 100
  zeta = u*1
  delta = rep(1,m)
  prob_inc = rep(0.5,m)
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  like_v = function(delta,i){
    if(all(delta == 0)){
      lval = 0
      return(lval)
    }
    ind_n0 = which(delta != 0)
    C = Istar[,ind_n0]
    phi_t = diag(1/cp,length(ind_n0),length(ind_n0))
    V_t = t(C)%*%C + solve(phi_t)
    mu = solve(V_t)%*%t(C)%*%y
    lval = -1/2*log((det(phi_t)))-1/2*log((det(V_t)))+1/2*t(mu)%*%V_t%*%mu
    return(lval)
  }
  
  for (k in 1:iter){
    # Step 1: update propensity score parameters
    # 1-1 : gamma
    data_gamma = cbind(data_obs$x1,data_obs$x2)
    w_gamma = apply(data_gamma%*%gamma_pos+Istar%*%zeta,1,rpg,n=1,h=1)
    k_gamma = data_obs$z-1/2
    L_gamma = k_gamma/w_gamma
    omega_gamma = diag(w_gamma)
    gamma_Sigma = solve(t(data_gamma)%*%omega_gamma%*%data_gamma + 1/var_gamma*diag(2))
    gamma_mu = gamma_Sigma%*%(t(data_gamma)%*%omega_gamma%*%(L_gamma-Istar%*%zeta))
    gamma_pos = mvrnorm(1,gamma_mu,gamma_Sigma)
    
    # 1-2 : var_gamma
    var_gamma = rinvgamma(1,shape=2/2+0.001,rate=0.001+t(gamma_pos)%*%gamma_pos/2)
    
    # 1-3 : zeta, 1,2,...m
    ## check variable inclusion
    ind_0 = which(delta==0) # zero
    ind_n0 = which(delta!=0) # nonzero
    C = Istar[,ind_n0]
    sigma_zeta = diag(1/50,length(ind_n0),length(ind_n0)) # specify 50 as variance
    # update zeta
    data_ps = cbind(data_obs$x1,data_obs$x2) #x1,x2
    w_cl = apply(data_ps%*%gamma_pos+Istar%*%zeta,1,rpg,n=1,h=1)
    k_cl = data_obs$z-1/2
    L_cl = k_cl/w_cl
    omega_cl = diag(w_cl)
    cl_Sigma = solve(sigma_zeta+t(C)%*%omega_cl%*%C)
    cl_mu = cl_Sigma%*%(t(C)%*%omega_cl%*%(L_cl-data_ps%*%gamma_pos))
    zeta[ind_n0] = mvrnorm(1,cl_mu,cl_Sigma)
    zeta[ind_0] = 0
    
    ## Sample v ##
    for (i in 1:m){
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
    Gamma_pos[k,] = gamma_pos
    Var_gamma[k,] =var_gamma
    Zeta[k,] = zeta
    Delta[k,] = delta
    Prob_inc[k,] = prob_inc
    
    print(k)
  }
  
  print(dataid)
  time_span = proc.time() - ptm
  return(list(Gamma_pos=Gamma_pos,Var_gamma=Var_gamma,Zeta=Zeta,Delta=Delta,Prob_inc=Prob_inc,time_span=time_span))
  
}

apply(Prob_inc[(burn+1):iter,],2,mean)
apply(Beta_pos[(burn+1):iter,],2,mean)
apply(Delta[(burn+1):iter,],2,mean)
