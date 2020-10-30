setwd("C:/Users/emmal/Desktop/bayesian sufficient")
rm(list=ls())
library(survival)
library(MASS)
library(parallel)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(dplyr)

psfunc=function(x,para){
  1/(1+exp(-x%*%para))
}

data_gene = function(dataid,m=50,m0=30,rhoxu=0,rhoyu=0,effect_true=2){
  set.seed(dataid)
  data = data.frame(x1=double(),x2=double(),U=double(),pstrue=double(),clu=integer(),
                    z = integer(),sufft=integer(),y0=double(),y1=double(),y_obs=double())
  clu_num = 1
  while (clu_num<=m0){
    n = sample(10:20,1,replace=T,prob=rep(1/11,11))
    x1 = rnorm(n,0,1)
    x2 = sample(c(-1,0,1),n,replace=T,prob=rep(1/3,3))
    data_sub = as.data.frame(cbind(x1,x2))
    data_sub$U = 0
    data_sub$pstrue = apply(as.matrix(data_sub),1,psfunc,para=c(1,1,1))
    data_sub = data_sub[data_sub$pstrue>=0.05 & data_sub$pstrue<=0.95,]
    if (dim(data_sub)[1]<=1){
      next
    }else{
      data_sub$clu = clu_num
      data_sub$z = apply(matrix(data_sub$pstrue),1,rbinom,n=1,size=1)
      suff_T = sum(data_sub$z)
      while (suff_T==dim(data_sub)[1] | suff_T==0){
        data_sub$z= apply(matrix(data_sub$pstrue),1,rbinom,n=1,size=1)
        suff_T = sum(data_sub$z)
      }
      data_sub$sufft = suff_T
      data_sub$y0 = 1 + data_sub$x1 + data_sub$x2 + rnorm(dim(data_sub)[1],0,1)
      data_sub$y1 = 1 + data_sub$x1 + data_sub$x2 + effect_true + rhoyu*data_sub$U + rnorm(dim(data_sub)[1],0,1)
      data_sub$y_obs = data_sub$y1*data_sub$z + data_sub$y0*(1-data_sub$z)
      data = rbind(data,data_sub)
      clu_num = clu_num+1
    }
  }
  
  while (clu_num<=m){
    n = sample(10:20,1,replace=T,prob=rep(1/11,11))
    x1 = rnorm(n,0,1)
    x2 = sample(c(-1,0,1),n,replace=T,prob=rep(1/3,3))
    u = rnorm(1,mean=-rhoxu*(mean(x1)+mean(x2)),sd=1)
    data_sub = as.data.frame(cbind(x1,x2))
    data_sub$U = u
    data_sub$pstrue = apply(as.matrix(data_sub),1,psfunc,para=c(1,1,1))
    data_sub = data_sub[data_sub$pstrue>=0.05 & data_sub$pstrue<=0.95,]
    if (dim(data_sub)[1]<=1){
      next
    }else{
      data_sub$clu = clu_num
      data_sub$z = apply(matrix(data_sub$pstrue),1,rbinom,n=1,size=1)
      suff_T = sum(data_sub$z)
      while (suff_T==dim(data_sub)[1] | suff_T==0){
        data_sub$z= apply(matrix(data_sub$pstrue),1,rbinom,n=1,size=1)
        suff_T = sum(data_sub$z)
      }
      data_sub$sufft = suff_T
      data_sub$y0 = 1 + data_sub$x1 + data_sub$x2 + rnorm(dim(data_sub)[1],0,1)
      data_sub$y1 = 1 + data_sub$x1 + data_sub$x2 + effect_true + rhoyu*data_sub$U + rnorm(dim(data_sub)[1],0,1)
      data_sub$y_obs = data_sub$y1*data_sub$z + data_sub$y0*(1-data_sub$z)
      data = rbind(data,data_sub)
      clu_num = clu_num+1
    }
  }
  n = rep(-99,m)
  suff_T = rep(-99,m)
  u = rep(-99,m)
  for (i in 1:m){
    n[i] = dim(data[data$clu==i,])[1]
    suff_T[i] = sum(data[data$clu==i,]$z)
    u[i] = data[data$clu==i,]$U[1]
  }
  id = unlist(apply(matrix(n),1,seq,from=1,by=1))
  data_obs = as.data.frame(cbind(data$clu,id,data$x1,data$x2,data$z,data$sufft,data$y_obs))
  colnames(data_obs) = c("clu","id","x1","x2","z","sufft","y_obs")
  ####### simualted ATE
  simu_true = mean(data$y1-data$y0)
  print (dataid)
  return (list(data=data,data_obs=data_obs,u=u,m=m,n=n,suff_T=suff_T,simu_true=simu_true))
}
simdatam100= mclapply(1:2,data_gene,m=50,m0=30,rhoxu=5,rhoyu=0,effect_true=2,mc.cores=1)

