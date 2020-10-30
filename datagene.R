setwd("C:/Users/emmal/Desktop/bayesian sufficient")
rm(list=ls())
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
#data_gene = function(dataid,m=50,m0=30,rhoxu=5,rhoyu=0,effect_true=2){
data_gene = function(dataid,m=500,rhoxu=0,rhoyu=0,effect_true=2){
  set.seed(dataid+11111)
  #n = as.matrix(sample(10:20,m,replace=T,prob=rep(1/11,11)))
  n = as.matrix(sample(2:5,m,replace=T,prob=rep(1/4,4)))
  id = unlist(apply(n,1,seq,from=1,by=1))
  clu = rep(1:m,n)
  data = as.data.frame(cbind(clu,id))
  N = dim(data)[1]
  data$x1 = rnorm(N,0,1)
  data$x2 = sample(c(-1,0,1),N,replace=T,prob=rep(1/3,3))
  x1mean = rep(-99,m) 
  x2mean = rep(-99,m)
  for (i in 1:m){
    x1mean[i] = mean(data[data$clu==i,]$x1)
    x2mean[i] = mean(data[data$clu==i,]$x2)
  }
  u = apply(as.matrix(-rhoxu*(x1mean+x2mean)),1,rnorm,n=1,sd=1)
  #u[1:m0]=0
  data$U = rep(u,n)
  psfunc=function(x,para){
    1/(1+exp(-x%*%para))
  }
  data$pstrue = apply(as.matrix(data[,3:5]),1,psfunc,para=c(1,1,1))
  
  ### remove too small/large pscore
  data = data[data$pstrue>=0.05 & data$pstrue<=0.95,]
  m = length(unique(data$clu))
  n = as.numeric(table(data$clu)) 
  n_ind = which(n>1)
  m_vec = c(1:m)[n_ind]
  data = data[data$clu %in%m_vec,]
  N = dim(data)[1]
  u = unique(data$U)
  m = length(unique(data$clu))
  n = as.numeric(table(data$clu)) 
  data$clu = rep(1:m,n)
  data$id = unlist(apply(matrix(n),1,seq,from=1,by=1))
  
  zij = apply(matrix(data$pstrue),1,rbinom,n=1,size=1)
  data$z = zij
  suff_T = rep(-99,m)
  for (i in 1:m){
    suff_T[i] = sum(data[data$clu==i,]$z)
  }
  ##### check each cluster, must have treatment and control both ####
  while (suff_T[1]==n[1] | suff_T[1]==0){
    sap = apply(as.matrix(data$pstrue[1:n[1]]),1,rbinom,n=1,size=1)
    data$z[1:n[1]] = sap
    suff_T[1] = sum(sap)
  }
  begin = 1
  total = n[1]
  for (i in 2:m){
    begin = begin + n[i-1]
    total = total+ n[i]
    while (suff_T[i]==n[i] | suff_T[i]==0){
      sap = apply(as.matrix(data$pstrue[begin:total]),1,rbinom,n=1,size=1) 
      data$z[begin:total] = sap
      suff_T[i] = sum(sap)
    }
  }
  data$sufft = rep(suff_T,n)
  data$y0 = 1 + data$x1 + data$x2 + rnorm(N,0,1)
  data$y1 = 1 + data$x1 + data$x2 + effect_true + rhoyu*data$U + rnorm(N,0,1)
  data$y_obs = data$y1*data$z + data$y0*(1-data$z)
  
  data_obs = data[,c(1:4,7,8,11)] # include int
  print(dataid)
  return (list(data=data,data_obs=data_obs,u=u,m=m,n=n,suff_T=suff_T))
}
simdatam20= mclapply(1:100,data_gene,m=500,rhoxu=5,rhoyu=0,effect_true=2,mc.cores=1)

sum=0
for (i in 1:2){
  n = simdatam20[[i]]$n
  suff = simdatam20[[i]]$suff_T
  if (sum(n==suff)!=0){
    sum=sum+1
  }
}
