rm(list=ls())
library(MASS)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
########################Part 1: Data generation #######################
data_gene = function(dataid,rhoxu = 0,rhoyu = 0,m = 500){
  n = as.matrix(sample(5:20,m,replace=T,prob=rep(1/16,16)))
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
  u = apply(as.matrix(-rhoxu*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  data$U = rep(u,n)
  ps=function(para){# intercept, x1,x2,U
    exp(para[1]+para[2]+para[3]+para[4])/(1+exp(para[1]+para[2]+para[3]+para[4]))
  }
  ps_true = apply(data[,3:6],1,ps)
  aij = apply(as.matrix(ps_true),1,rbinom,n=1,size=1)
  data$A = aij
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data[data$clu==i,]$A)
  }
  while (suff_T[1]==n[1] | suff_T[1]==0){
    sap = apply(as.matrix(ps_true[1:n[1]]),1,rbinom,n=1,size=1)
    data$A[1:n[1]] = sap
    suff_T[1] = sum(sap)
  }
  begin = 1
  total = n[1]
  for (i in 2:m){
    begin = begin + n[i-1]
    total = total+ n[i]
    while (suff_T[i]==n[i] | suff_T[i]==0){
      sap = apply(as.matrix(ps_true[begin:total]),1,rbinom,n=1,size=1) 
      data$A[begin:total] = sap
      suff_T[i] = sum(sap)
    }
  }
  data$sufft = rep(suff_T,n)
  data$pstrue = ps_true
  data$y0 = 1 + data$x1 + data$x2 + rhoyu*data$U+ rnorm(N,0,1)
  data$y1 = 1 + data$x1 + data$x2 + 2 + rhoyu*data$U + rnorm(N,0,1)
  data$y_obs = data$y1*data$A + data$y0*(1-data$A)
  effec_simu = mean(data$y1-data$y0)
  effec_naiv = mean(data$A*data$y_obs-(1-data$A)*data$y_obs)
  data_obs = data[,c(1:5,7,8,12)]
  return(list(data=data,data_obs=data_obs,N=N,m=m,n=n,u=u,suff_T = suff_T,effec_simu=effec_simu,effec_naiv=effec_naiv))
}
simudata = lapply(c(1:5),data_gene,rhoxu = 0,rhoyu = 5,m = 20)

data = simudata[[dataid]]$data
data_obs = simudata[[dataid]]$data_obs
N = simudata[[dataid]]$N
m = simudata[[dataid]]$m
n = simudata[[dataid]]$n
u = simudata[[dataid]]$u
suff_T = simudata[[dataid]]$suff_T

###################### part 2: estimated pscore usign sufficient statistic
beta_est=c(1,1,1)
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
ps_suffi2=prop_ps[-1]
data$suffi_weight = ifelse(data$A==1,ps_suffi2,1-ps_suffi2)