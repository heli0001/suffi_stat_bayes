rm(list=ls())
library(MASS)
library(survival)
library(tictoc)
library(utils)
memory.limit(size = 10000000000000)
library(RcppAlgos)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(parallel)
library(gtools)
########################Part 1: Data generation #######################
rhoxu1=0
rhoxu2=0
rhoyu1=0
rhoyu2=0
m=500
data_gene = function(dataid,rhoxu1 = 0,rhoxu2 = 0,rhoyu1 = 0,rhoyu2 = 0,m = 500){
  n = as.matrix(sample(3:5,m,replace=T,prob=rep(1/3,3)))
  id = (unlist(apply(n,1,seq,from=1,by=1)))
  clu = rep(1:m,n)
  data = as.data.frame(cbind(clu,id))
  N = dim(data)[1]
  inter = rep(1,N)
  data$Int = inter
  data$x1 = rnorm(N,0,1)
  data$x2 = sample(c(-1,0,1),N,replace=T,prob=rep(1/3,3))
  x1mean = rep(-99,m) 
  x2mean = rep(-99,m)
  for (i in 1:m){
    x1mean[i] = mean(data[data$clu==i,]$x1)
    x2mean[i] = mean(data[data$clu==i,]$x2)
  }
  u1 = apply(as.matrix(-rhoxu1*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  u2 = apply(as.matrix(-rhoxu2*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  data$U1 = rep(u1,n)
  data$U2 = rep(u2,n)
  ps0 = function(para){ # includes intercept
    1/(exp(1+para[1]+para[2]+para[3])+exp(1+para[1]+para[2]+para[4])+1)
  }
  ps1 = function(para){
    exp(1+para[1]+para[2]+para[3])/(exp(1+para[1]+para[2]+para[3])+exp(1+para[1]+para[2]+para[4])+1)
  }
  pstrue = matrix(-99,ncol=3,nrow=N)
  pstrue[,1] = apply(data[,4:7],1,ps0)
  pstrue[,2] = apply(data[,4:7],1,ps1)
  pstrue[,3] = 1-pstrue[,1]- pstrue[,2]
  aij = apply(pstrue,1,sample,x=c(0,1,2),size=1,replace=F)
  data$A = aij
  ps_true = ifelse(data$A==0,pstrue[,1],ifelse(data$A==1,pstrue[,2],pstrue[,3]))
  suffi = matrix(-99,nrow=m,ncol=2)
  for (i in 1:m){
    cluA = data[data$clu==i,]$A
    suffi[i,1] = length(cluA[cluA==1])
    suffi[i,2] = length(cluA[cluA==2])
  }
  while (suffi[1,1]==0|suffi[1,2]==0|sum(suffi[1,])==n[1] | sum(suffi[1,])==0){
    sap = apply(as.matrix(pstrue[1:n[1],]),1,sample,x=c(0,1,2),size=1,replace=F)
    data$A[1:n[1]] = sap
    suffi[1,] = c(length(sap[sap==1]),length(sap[sap==2]))
  }
  begin = 1
  total = n[1]
  for (i in 2:m){
    begin = begin + n[i-1]
    total = total+ n[i]
    while (suffi[i,1]==0|suffi[i,2]==0|sum(suffi[i,])==n[i] | sum(suffi[i,])==0){
      sap = apply(as.matrix(pstrue[begin:total,]),1,sample,x=c(0,1,2),size=1,replace=F) 
      data$A[begin:total] = sap
      suffi[i,] = c(length(sap[sap==1]),length(sap[sap==2]))
    }
  }
  data$sufft1 = rep(suffi[,1],n)
  data$sufft2 = rep(suffi[,2],n)
  data$y0 = 1 + data$x1 + data$x2 + rhoyu1*data$U1 +rhoyu2*data$U2+rnorm(N,0,1)
  data$y1 = 1 + data$x1 + data$x2 + 2 + rhoyu1*data$U1 +rhoyu2*data$U2 + rnorm(N,0,1)
  data$y2 = 1 + data$x1 + data$x2 + 4 + rhoyu1*data$U1 +rhoyu2*data$U2 + rnorm(N,0,1)
  data$y_obs = ifelse(data$A==0,data$y0,ifelse(data$A==1,data$y1,data$y2))
  effec_true = c(2,4,2) #effect 01,02,12
  effec_simu = c(mean(data$y1-data$y0),mean(data$y2-data$y0),mean(data$y2-data$y1))
  data_obs = data[,c(1:5,8:10,14)]
  print(dataid)
  return(list(data=data,data_obs = data_obs,m=m,n=n,N=N,u1=u1,u2=u2,suffi=suffi,pstrue=pstrue,ps_true=ps_true,effec_simu=effec_simu))
}
simudata = mclapply(1,data_gene,rhoxu1 = 3,rhoxu2 = 3,rhoyu1 = 3,rhoyu2 = 3,m = 500,mc.cores=1)

psfunc_naive = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  data01 = data_obs[data_obs$A==0|data_obs$A==1,]
  data02 = data_obs[data_obs$A==0|data_obs$A==2,]
  data12 = data_obs[data_obs$A==1|data_obs$A==2,]
  ########### group 0 vs 1 ###########
  N = dim(data01)[1]
  m = length(unique(data01$clu)) ## same as original m
  n = as.vector(table(data01$clu))
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  var_betaps = rep(-99,iter)
  beta_ps_pos[1,] = rnorm(3,0,1)
  var_betaps[1]=100
  data_ps = as.matrix(data01[,3:5])
  for (g in 2:iter){
    set.seed(g)
    ## update beta
    wij = apply(data_ps%*%beta_ps_pos[g-1,],1,rpg,n=1,h=1) 
    kij = data01$A-1/2
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/var_betaps[g-1]*diag(3))
    mu_ps = var_ps %*%t(data_ps)%*%kij
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    ## update variance
    var_betaps[g] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(beta_ps_pos[g,])%*%beta_ps_pos[g,]/2)
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  data01$w_naive = apply(as.matrix(data01[,3:5]),1,psfunc,para=beta_ps_pos_mean)
  effec01_naive = mean(data01$A*data01$y_obs/data01$w_naive-(1-data01$A)*data01$y_obs/(1-data01$w_naive))
  poste01 = list(data01=data01,beta_ps_pos=beta_ps_pos,effec01_naive=effec01_naive)
  ########### group 0 vs 2 ###########
  N = dim(data02)[1]
  data02$A = ifelse(data02$A==2,1,0)
  m = length(unique(data02$clu)) ## same as original m
  n = as.vector(table(data02$clu))
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  var_betaps = rep(-99,iter)
  beta_ps_pos[1,] = rnorm(3,0,1)
  var_betaps[1]=100
  data_ps = as.matrix(data02[,3:5])
  for (g in 2:iter){
    set.seed(g)
    ## update beta
    wij = apply(data_ps%*%beta_ps_pos[g-1,],1,rpg,n=1,h=1) 
    kij = data02$A-1/2
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/var_betaps[g-1]*diag(3))
    mu_ps = var_ps %*%t(data_ps)%*%kij
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    ## update variance
    var_betaps[g] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(beta_ps_pos[g,])%*%beta_ps_pos[g,]/2)
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  data02$w_naive = apply(as.matrix(data02[,3:5]),1,psfunc,para=beta_ps_pos_mean)
  effec02_naive = mean(data02$A*data02$y_obs/data02$w_naive-(1-data02$A)*data02$y_obs/(1-data02$w_naive))
  poste02 = list(data02=data02,beta_ps_pos=beta_ps_pos,effec02_naive=effec02_naive)
  
  
  
  
  
  

  
  
  
  
  
  ########### group 1 vs 2 ###########
  N = dim(data12)[1]
  data12$A = ifelse(data12$A==1,0,1)
  m = length(unique(data12$clu)) ## same as original m
  n = as.vector(table(data12$clu))
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  var_betaps = rep(-99,iter)
  beta_ps_pos[1,] = rnorm(3,0,1)
  var_betaps[1]=100
  data_ps = as.matrix(data12[,3:5])
  for (g in 2:iter){
    set.seed(g)
    ## update beta
    wij = apply(data_ps%*%beta_ps_pos[g-1,],1,rpg,n=1,h=1) 
    kij = data12$A-1/2
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/var_betaps[g-1]*diag(3))
    mu_ps = var_ps %*%t(data_ps)%*%kij
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    ## update variance
    var_betaps[g] = rinvgamma(1,shape=0.001+3/2,rate=0.001+t(beta_ps_pos[g,])%*%beta_ps_pos[g,]/2)
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  data12$w_naive = apply(as.matrix(data12[,3:5]),1,psfunc,para=beta_ps_pos_mean)
  effec12_naive = mean(data12$A*data12$y_obs/data12$w_naive-(1-data12$A)*data12$y_obs/(1-data12$w_naive))
  poste12 = list(data12=data12,beta_ps_pos=beta_ps_pos,effec12_naive=effec12_naive)
  print(dataid)
  return(list(poste01=poste01,poste02=poste02,poste12=poste12))
}
naiveeffect = lapply(1,psfunc_naive,iter=20,burn=0)

psfunc_ran = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  data01 = data_obs[data_obs$A==0|data_obs$A==1,]
  data02 = data_obs[data_obs$A==0|data_obs$A==2,]
  data12 = data_obs[data_obs$A==1|data_obs$A==2,]
  ########### group 0 vs 1 ######################
  N = dim(data01)[1]
  m = length(unique(data01$clu)) ## same as original m
  n = as.vector(table(data01$clu))
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data01[data01$clu==i,]$A)
  }
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  alpha_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha
  beta_ps_pos[1,] = rep(1,3)
  alpha_ps_pos[1,] = mvrnorm(1,mu=rep(0,m),Sigma = diag(m))
  var_alpha_pos[1] = 100
  data_ps = as.matrix(data01[,3:5])
  ## dummy variable/indicator matrix 
  cluind = data01$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (g in 2:iter){
    ## sample alpha_ps
    for (cl in 1:m){
      data_j = as.matrix(data01[data01$clu==cl,3:5])
      ni = dim(data_j)[1]
      waij = apply(data_j%*%beta_ps_pos[g-1,]+alpha_ps_pos[g-1,cl],1,rpg,n=1,h=1)
      kaij = data01[data01$clu==cl,]$A-1/2
      La = kaij/waij
      sig = diag(waij)
      var_a = 1/(1/var_alpha_pos[g-1]+rep(1,ni)%*%sig%*%c(rep(1,ni)))
      mu_a = var_a*(rep(1,ni)%*%sig%*%(La-data_j%*%beta_ps_pos[g-1,]))
      alpha_ps_pos[g,cl] = rnorm(1,mu_a,sqrt(var_a))
    }
    ## sample beta_ps
    wij = apply(data_ps%*%beta_ps_pos[g-1,]+Istar%*%alpha_ps_pos[g,],1,rpg,n=1,h=1)
    kij = data01$A-1/2
    L = kij/wij
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/100*diag(3))
    mu_ps = var_ps %*%(t(data_ps)%*%omega%*%(L-Istar%*%alpha_ps_pos[g,]))
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    
    ## sample var_alpha
    ralpha = sum((alpha_ps_pos[g,])^2)
    shape1 = (0.002+m)/2
    scale1 = (0.002*1 + ralpha)/2
    var_alpha_pos[g] = rinvgamma(1,shape=shape1,rate=scale1)
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  alpha_ps_pos_mean = apply(alpha_ps_pos,2,mean)
  alphaij_mean = Istar%*%alpha_ps_pos_mean
  psfunc2 = function(x,para1){
    exp(x[1:3]%*%para1+x[4])/(1+exp(x[1:3]%*%para1+x[4]))
  }
  data_ps2 = as.matrix(cbind(data_ps,alphaij_mean))
  data01$w_ran = apply(data_ps2,1,psfunc2,para1=beta_ps_pos_mean)
  effec01_ran = mean(data01$A*data01$y_obs/data01$w_ran-(1-data01$A)*data01$y_obs/(1-data01$w_ran))
  
  poste01 = list(data01=data01,beta_ps_pos=beta_ps_pos,alpha_ps_pos=alpha_ps_pos,var_alpha_pos=var_alpha_pos,effec01_ran=effec01_ran)
  ########### group 0 vs 2 ######################
  N = dim(data02)[1]
  data02$A = ifelse(data02$A==2,1,0)
  m = length(unique(data02$clu)) ## same as original m
  n = as.vector(table(data02$clu))
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data02[data02$clu==i,]$A)
  }
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  alpha_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha
  beta_ps_pos[1,] = rep(1,3)
  alpha_ps_pos[1,] = mvrnorm(1,mu=rep(0,m),Sigma = diag(m))
  var_alpha_pos[1] = 100
  data_ps = as.matrix(data02[,3:5])
  ## dummy variable/indicator matrix 
  cluind = data02$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (g in 2:iter){
    ## sample alpha_ps
    for (cl in 1:m){
      data_j = as.matrix(data02[data02$clu==cl,3:5])
      ni = dim(data_j)[1]
      waij = apply(data_j%*%beta_ps_pos[g-1,]+alpha_ps_pos[g-1,cl],1,rpg,n=1,h=1)
      kaij = data02[data02$clu==cl,]$A-1/2
      La = kaij/waij
      sig = diag(waij)
      var_a = 1/(1/var_alpha_pos[g-1]+rep(1,ni)%*%sig%*%c(rep(1,ni)))
      mu_a = var_a*(rep(1,ni)%*%sig%*%(La-data_j%*%beta_ps_pos[g-1,]))
      alpha_ps_pos[g,cl] = rnorm(1,mu_a,sqrt(var_a))
    }
    ## sample beta_ps
    wij = apply(data_ps%*%beta_ps_pos[g-1,]+Istar%*%alpha_ps_pos[g,],1,rpg,n=1,h=1)
    kij = data02$A-1/2
    L = kij/wij
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/100*diag(3))
    mu_ps = var_ps %*%(t(data_ps)%*%omega%*%(L-Istar%*%alpha_ps_pos[g,]))
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    
    ## sample var_alpha
    ralpha = sum((alpha_ps_pos[g,])^2)
    shape1 = (0.002+m)/2
    scale1 = (0.002*1 + ralpha)/2
    var_alpha_pos[g] = rinvgamma(1,shape=shape1,rate=scale1)
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  alpha_ps_pos_mean = apply(alpha_ps_pos,2,mean)
  alphaij_mean = Istar%*%alpha_ps_pos_mean
  psfunc2 = function(x,para1){
    exp(x[1:3]%*%para1+x[4])/(1+exp(x[1:3]%*%para1+x[4]))
  }
  data_ps2 = as.matrix(cbind(data_ps,alphaij_mean))
  data02$w_ran = apply(data_ps2,1,psfunc2,para1=beta_ps_pos_mean)
  effec02_ran = mean(data02$A*data02$y_obs/data02$w_ran-(1-data02$A)*data02$y_obs/(1-data02$w_ran))
  
  poste02 = list(data02=data02,beta_ps_pos=beta_ps_pos,alpha_ps_pos=alpha_ps_pos,var_alpha_pos=var_alpha_pos,effec02_ran=effec02_ran)
  ########### group 1 vs 2 ######################
  N = dim(data12)[1]
  data12$A = ifelse(data12$A==1,0,1)
  m = length(unique(data12$clu)) ## same as original m
  n = as.vector(table(data12$clu))
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data12[data12$clu==i,]$A)
  }
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  alpha_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha
  beta_ps_pos[1,] = rep(1,3)
  alpha_ps_pos[1,] = mvrnorm(1,mu=rep(0,m),Sigma = diag(m))
  var_alpha_pos[1] = 100
  data_ps = as.matrix(data12[,3:5])
  ## dummy variable/indicator matrix 
  cluind = data12$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (g in 2:iter){
    ## sample alpha_ps
    for (cl in 1:m){
      data_j = as.matrix(data12[data12$clu==cl,3:5])
      ni = dim(data_j)[1]
      waij = apply(data_j%*%beta_ps_pos[g-1,]+alpha_ps_pos[g-1,cl],1,rpg,n=1,h=1)
      kaij = data12[data12$clu==cl,]$A-1/2
      La = kaij/waij
      sig = diag(waij)
      var_a = 1/(1/var_alpha_pos[g-1]+rep(1,ni)%*%sig%*%c(rep(1,ni)))
      mu_a = var_a*(rep(1,ni)%*%sig%*%(La-data_j%*%beta_ps_pos[g-1,]))
      alpha_ps_pos[g,cl] = rnorm(1,mu_a,sqrt(var_a))
    }
    ## sample beta_ps
    wij = apply(data_ps%*%beta_ps_pos[g-1,]+Istar%*%alpha_ps_pos[g,],1,rpg,n=1,h=1)
    kij = data12$A-1/2
    L = kij/wij
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/100*diag(3))
    mu_ps = var_ps %*%(t(data_ps)%*%omega%*%(L-Istar%*%alpha_ps_pos[g,]))
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
    
    ## sample var_alpha
    ralpha = sum((alpha_ps_pos[g,])^2)
    shape1 = (0.002+m)/2
    scale1 = (0.002*1 + ralpha)/2
    var_alpha_pos[g] = rinvgamma(1,shape=shape1,rate=scale1)
    print(g)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  alpha_ps_pos_mean = apply(alpha_ps_pos,2,mean)
  alphaij_mean = Istar%*%alpha_ps_pos_mean
  psfunc2 = function(x,para1){
    exp(x[1:3]%*%para1+x[4])/(1+exp(x[1:3]%*%para1+x[4]))
  }
  data_ps2 = as.matrix(cbind(data_ps,alphaij_mean))
  data12$w_ran = apply(data_ps2,1,psfunc2,para1=beta_ps_pos_mean)
  effec12_ran = mean(data12$A*data12$y_obs/data12$w_ran-(1-data12$A)*data12$y_obs/(1-data12$w_ran))
  
  
  poste12 = list(data12=data12,beta_ps_pos=beta_ps_pos,alpha_ps_pos=alpha_ps_pos,var_alpha_pos=var_alpha_pos,effec12_ran=effec12_ran)
  
  print(dataid)
  return(list(poste01=poste01,poste02=poste02,poste12=poste12))
  }
raneffect = lapply(1,psfunc_ran,iter=20,burn=0)

psfunc_suffi = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  data01 = data_obs[data_obs$A==0|data_obs$A==1,]
  data02 = data_obs[data_obs$A==0|data_obs$A==2,]
  data12 = data_obs[data_obs$A==1|data_obs$A==2,]
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  suffi_method = function(beta_est,data,m=m,n=n,suff_T=suff_T){
    dem=rep(-99,m)
    prop_ps = 0
    for (cls in 1:m){
      ## denominator part
      data_clu = as.matrix(data[data$clu==cls,3:5])
      ni=n[cls]
      suffT =suff_T[cls]
      rdenom = matrix(-99,nrow=suffT+1,ncol=ni+1)
      rdenom[lower.tri(rdenom)]=0 ## lower part =0
      rdenom[1,]=1
      for (rowid in 2:(suffT+1)){
        for (colid in 2:(ni+1)){
          rdenom[rowid,colid]=rdenom[rowid,colid-1] + exp(data_clu[colid-1,]%*%beta_est)*rdenom[rowid-1,colid-1]
        }
      }
      dem[cls] = rdenom[suffT+1,ni+1]
      for (sub in 1:n[cls]){
        ## numerator part
        rnum = matrix(-99,nrow=suffT+1,ncol=ni+1)
        rnum[lower.tri(rnum)]=0 ## lower part =0
        rnum[1,]=1
        for (rowid in 2:(suffT+1)){
          for (colid in 2:(ni+1)){
            rnum[rowid,colid]=rnum[rowid,colid-1] + exp(data_clu[colid-1,]%*%beta_est)*(1-as.numeric(colid-1==sub))*rnum[rowid-1,colid-1]
          }
        }
        if (data[data$clu==cls,]$A[sub]==0){
          nume= rnum[suffT+1,ni+1]
        }else{
          nume= dem[cls]-rnum[suffT+1,ni+1]
        }
        
        prob = nume/dem[cls]
        prop_ps = rbind(prop_ps,prob)
      }
    }
    ps_suffi=prop_ps[-1]
    return (ps_suffi)
  }
  #full joint log posterior of beta|y for the logistic regression model:
  posteriorfunc=function(param,data_obs,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=1){
    psall1_suffi = suffi_method(param,data_obs,m=m,n=n,suff_T=suff_T)
    ps_suffi = ifelse(data_obs$A==1,psall1_suffi,1-psall1_suffi)
    # like = 0
    # for (i in 1:N){
    #   like=like+dbinom(data_obs$A[i],size=1,prob=ps_suffi[i],log=TRUE) 
    # }
    like = sum(dbinom(data_obs$A,1,ps_suffi,log=TRUE))
    prior = sum(dnorm(param,prior.mn,prior.sd,log=TRUE))
    return(like+prior)
  }
  
  ############### group 0 vs 1 #########################
  N = dim(data01)[1]
  m = length(unique(data01$clu)) ## same as original m
  n = as.vector(table(data01$clu))
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data01[data01$clu==i,]$A)
  }
  # Start sampling
  #Initial values:
  beta = rnorm(3,0,1)
  beta_ps_pos = matrix(-99,iter,3)
  acc = rep(0,3)
  cur_log_post = posteriorfunc(beta,data01,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=10)
  for(g in 1:iter){
    #Update beta using MH sampling:
    for(j in 1:3){
      # Draw candidate:
      canbeta = beta
      canbeta[j] = rnorm(1,beta[j],1)
      can_log_post = posteriorfunc(canbeta,data01,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=10)
      # Compute acceptance ratio:
      R = (can_log_post-cur_log_post)
      if(log(runif(1))<R){
        beta = canbeta
        cur_log_post = can_log_post
        acc[j] = acc[j]+1
      }
    }
    print(g)
    print(acc)
    beta_ps_pos[g,]=beta
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  ps_suffi_mean = suffi_method(beta_ps_pos_mean,data01,m=m,n=n,suff_T=suff_T)
  data01$w_suffi = ifelse(data01$A==1,ps_suffi_mean,1-ps_suffi_mean)
  effec01_suffi = mean(data01$A*data01$y_obs/data01$w_suffi-(1-data01$A)*data01$y_obs/(1-data01$w_suffi))
  poste01 = list(data01=data01,beta_ps_pos=beta_ps_pos,effec01_suffi=effec01_suffi)
  
  ############### group 0 vs 2 #########################
  N = dim(data02)[1]
  data02$A = ifelse(data02$A==2,1,0)
  m = length(unique(data02$clu)) ## same as original m
  n = as.vector(table(data02$clu))
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data02[data02$clu==i,]$A)
  }
  # Start sampling
  #Initial values:
  beta = rnorm(3,0,1)
  beta_ps_pos = matrix(-99,iter,3)
  acc = rep(0,3)
  cur_log_post = posteriorfunc(beta,data02,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=10)
  for(g in 1:iter){
    #Update beta using MH sampling:
    for(j in 1:3){
      # Draw candidate:
      canbeta = beta
      canbeta[j] = rnorm(1,beta[j],1)
      can_log_post = posteriorfunc(canbeta,data02,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=10)
      # Compute acceptance ratio:
      R = (can_log_post-cur_log_post)
      if(log(runif(1))<R){
        beta = canbeta
        cur_log_post = can_log_post
        acc[j] = acc[j]+1
      }
    }
    print(g)
    print(acc)
    beta_ps_pos[g,]=beta
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  ps_suffi_mean = suffi_method(beta_ps_pos_mean,data02,m=m,n=n,suff_T=suff_T)
  data02$w_suffi = ifelse(data02$A==1,ps_suffi_mean,1-ps_suffi_mean)
  effec02_suffi = mean(data02$A*data02$y_obs/data02$w_suffi-(1-data02$A)*data02$y_obs/(1-data02$w_suffi))
  poste02 = list(data02=data02,beta_ps_pos=beta_ps_pos,effec02_suffi=effec02_suffi)
  
  ############### group 1 vs 2 #########################
  N = dim(data12)[1]
  data12$A = ifelse(data12$A==2,1,0)
  m = length(unique(data12$clu)) ## same as original m
  n = as.vector(table(data12$clu))
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data12[data12$clu==i,]$A)
  }
  # Start sampling
  #Initial values:
  beta = rnorm(3,0,1)
  beta_ps_pos = matrix(-99,iter,3)
  acc = rep(0,3)
  cur_log_post = posteriorfunc(beta,data12,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=10)
  for(g in 1:iter){
    #Update beta using MH sampling:
    for(j in 1:3){
      # Draw candidate:
      canbeta = beta
      canbeta[j] = rnorm(1,beta[j],1)
      can_log_post = posteriorfunc(canbeta,data12,m=m,n=n,suff_T=suff_T,prior.mn=0,prior.sd=10)
      # Compute acceptance ratio:
      R = (can_log_post-cur_log_post)
      if(log(runif(1))<R){
        beta = canbeta
        cur_log_post = can_log_post
        acc[j] = acc[j]+1
      }
    }
    print(g)
    print(acc)
    beta_ps_pos[g,]=beta
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  ps_suffi_mean = suffi_method(beta_ps_pos_mean,data12,m=m,n=n,suff_T=suff_T)
  data12$w_suffi = ifelse(data12$A==1,ps_suffi_mean,1-ps_suffi_mean)
  effec12_suffi = mean(data12$A*data12$y_obs/data12$w_suffi-(1-data12$A)*data12$y_obs/(1-data12$w_suffi))
  poste12 = list(data12=data12,beta_ps_pos=beta_ps_pos,effec12_suffi=effec12_suffi)
  
  print(dataid)
  return(list(poste01=poste01,poste02=poste02,poste12=poste12))
  
}
suffieffect = lapply(1,psfunc_suffi,iter=100,burn=10)
