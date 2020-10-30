rm(list=ls())
library(MASS)
library(survival)
library(tictoc)
library(Matrix)
library(mclogit)
library(RcppAlgos)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
library(parallel)
library(gtools)
rhoxu1=0
rhoxu2=0
rhoyu1=0
rhoyu2=0
m=500
data_gene = function(dataid,rhoxu1 = 0,rhoxu2 = 0,rhoyu1 = 0,rhoyu2 = 0,m = 100){
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
  ps0 = function(para){
    1/(exp(para[1]+para[2]+para[3])+exp(para[1]+para[2]+para[4])+1)
  }
  ps1 = function(para){
    exp(para[1]+para[2]+para[3])/(exp(para[1]+para[2]+para[3])+exp(para[1]+para[2]+para[4])+1)
  }
  pstrue = matrix(-99,ncol=3,nrow=N)
  pstrue[,1] = apply(data[,4:7],1,ps0)
  pstrue[,2] = apply(data[,4:7],1,ps1)
  pstrue[,3] = 1-pstrue[,1]- pstrue[,2]
  aij = apply(pstrue,1,sample,x=c(0,1,2),size=1,replace=F)
  data$A = aij
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
  
  ###### remove the sufficient statistics ==3 cluster, No.of each group can not exceed n/2+1 within each cluster
  cl_rm1 = which(suffi[,1]==3)
  set1 = setdiff(c(1:500),cl_rm1)
  cl_rm2 = which(suffi[,2]==3)
  set2 = setdiff(set1,cl_rm2)
  cl_rm3 = which(n-(suffi[,1]+suffi[,2])==3)
  set3 = setdiff(set2,cl_rm3)
  data = data[data$clu %in% set3,] ## updated dataset
  # m=length(set3)
  # n = n[set3]
  # N=dim(data)[1]
  # suffi = suffi[set3,]
  effec_true = c(2,4,2) #effect 01,02,12
  effec_simu = c(mean(data$y1-data$y0),mean(data$y2-data$y0),mean(data$y2-data$y1))
  data_obs = data[,c(1:5,8:10,14)]
  print(dataid)
  return(list(set3=set3,data=data,data_obs = data_obs,m=m,n=n,N=N,u1=u1,u2=u2,suffi=suffi,pstrue=pstrue,effec_simu=effec_simu))
}
simudata = mclapply(1:3,data_gene,rhoxu1 = 0,rhoxu2 = 0,rhoyu1 = 0,rhoyu2 = 0,m = 500,mc.cores=1)

# if multinomial, could use mclogit for conditional MLE
#library(mclogit)
#multlog_summ = mclogit(cbind(clu, A)~x1+x2, data=data)
#beta_mle = multlog_summ$coefficients
data = simudata[[1]]$data
data_obs = simudata[[1]]$data_obs
N=simudata[[1]]$N
n = simudata[[1]]$n
m = simudata[[1]]$m
u1 = simudata[[1]]$u1
u2 = simudata[[1]]$u2
suffi=simudata[[1]]$suffi
beta1_mle=c(1,1)
beta2_mle=c(1,1)
prop_ps = 0
for (cl in set3){
  data_c = as.matrix(data[data$clu==cl,c(4:5,8)]) ## x1,x2,A
  ########### calculate denominator
  data_den = data_c[,1:2]
  Vbar = permuteGeneral(c(0,1,2), freqs = c(n[cl]-sum(suffi[cl,]),suffi[cl,1],suffi[cl,2]))
  sum_den = rep(-99,dim(Vbar)[1])
  sum_abar = rep(-99,n[cl])
  for (k in 1:dim(Vbar)[1]){
    abar = Vbar[k,]
    Iabar = matrix(-99,nrow=n[cl],ncol=3)
    for (q in 1:n[cl]){
      Iabar[q,1] = ifelse(abar[q]==1,1,0)
      Iabar[q,2] = ifelse(abar[q]==2,1,0)
      Iabar[q,3] = ifelse(abar[q]==0,1,0)
    }
    # sum_abar = sum(Iabar%*%d_need)
    for (p in 1:n[cl]){
      sum_abar[p] = Iabar[p,]%*%c(data_den[p,]%*%beta1_mle,data_den[p,]%*%beta2_mle,1)
    }
    sum_den[k] = exp(sum(sum_abar))
  }
  for (sub in 1:n[cl]){
    data_num = data_c[-sub,1:2]
    Vij = Vbar[Vbar[,sub]==data_c[sub,3],]
    ########## calculate numerator
    #xij = data_c[sub,1:2]
    #d_need = c(xij%*%beta1_mle,xij%*%beta2_mle,1)
    sum_num = rep(-99,dim(Vij)[1])
    sum_astar = rep (-99,n[cl]-1)
    for (k in 1:dim(Vij)[1]){
      astar_j = Vij[k,-sub]
      Iasta = matrix(-99,nrow=n[cl]-1,ncol=3)
      for (q in 1:(n[cl]-1)){
        Iasta[q,1] = ifelse(astar_j[q]==1,1,0)
        Iasta[q,2] = ifelse(astar_j[q]==2,1,0)
        Iasta[q,3] = ifelse(astar_j[q]==0,1,0)
      }
      # sum_astar = sum(Iasta%*%d_need)
      for (p in 1:(n[cl]-1)){
        sum_astar[p] = Iasta[p,]%*%c(data_num[p,]%*%beta1_mle,data_num[p,]%*%beta2_mle,1)
      }
      sum_num[k] = exp(sum(sum_astar))
    }
    Isub = rep(-99,3)
    Isub[1] = ifelse(data_c[sub,3]==1,1,0)
    Isub[2] = ifelse(data_c[sub,3]==2,1,0)
    Isub[3] = ifelse(data_c[sub,3]==0,1,0)
    # num_sub = Isub%*%d_need
    num_sub = Isub%*%c(data_c[sub,1:2]%*%beta1_mle,data_c[sub,1:2]%*%beta2_mle,1)
    num = exp(num_sub)*sum(sum_num)
    
    prob = num/sum(sum_den)
    prop_ps = rbind(prop_ps,prob)  
  }
}
ps_suffi = prop_ps[-1,]

############# to calcuate the probability for other groups ###########
updateA = c()
for(i in set3){
  renew = sample(data[data$clu==i,]$A)
  while (sum(renew==data[data$clu==i,]$A)!=0){
    renew = sample(data[data$clu==i,]$A)
  }
  print(i)
  updateA=c(updateA,renew)
}

datanew = data
datanew$A = updateA
prop_ps = 0
for (cl in set3){
  data_c = as.matrix(datanew[datanew$clu==cl,c(4:5,8)]) ## x1,x2,A
  ########### calculate denominator
  data_den = data_c[,1:2]
  Vbar = permuteGeneral(c(0,1,2), freqs = c(n[cl]-sum(suffi[cl,]),suffi[cl,1],suffi[cl,2]))
  sum_den = rep(-99,dim(Vbar)[1])
  sum_abar = rep(-99,n[cl])
  for (k in 1:dim(Vbar)[1]){
    abar = Vbar[k,]
    Iabar = matrix(-99,nrow=n[cl],ncol=3)
    for (q in 1:n[cl]){
      Iabar[q,1] = ifelse(abar[q]==1,1,0)
      Iabar[q,2] = ifelse(abar[q]==2,1,0)
      Iabar[q,3] = ifelse(abar[q]==0,1,0)
    }
    # sum_abar = sum(Iabar%*%d_need)
    for (p in 1:n[cl]){
      sum_abar[p] = Iabar[p,]%*%c(data_den[p,]%*%beta1_mle,data_den[p,]%*%beta2_mle,1)
    }
    sum_den[k] = exp(sum(sum_abar))
  }
  for (sub in 1:n[cl]){
    data_num = data_c[-sub,1:2]
    Vij = Vbar[Vbar[,sub]==data_c[sub,3],]
    ########## calculate numerator
    #xij = data_c[sub,1:2]
    #d_need = c(xij%*%beta1_mle,xij%*%beta2_mle,1)
    sum_num = rep(-99,dim(Vij)[1])
    sum_astar = rep (-99,n[cl]-1)
    for (k in 1:dim(Vij)[1]){
      astar_j = Vij[k,-sub]
      Iasta = matrix(-99,nrow=n[cl]-1,ncol=3)
      for (q in 1:(n[cl]-1)){
        Iasta[q,1] = ifelse(astar_j[q]==1,1,0)
        Iasta[q,2] = ifelse(astar_j[q]==2,1,0)
        Iasta[q,3] = ifelse(astar_j[q]==0,1,0)
      }
      # sum_astar = sum(Iasta%*%d_need)
      for (p in 1:(n[cl]-1)){
        sum_astar[p] = Iasta[p,]%*%c(data_num[p,]%*%beta1_mle,data_num[p,]%*%beta2_mle,1)
      }
      sum_num[k] = exp(sum(sum_astar))
    }
    Isub = rep(-99,3)
    Isub[1] = ifelse(data_c[sub,3]==1,1,0)
    Isub[2] = ifelse(data_c[sub,3]==2,1,0)
    Isub[3] = ifelse(data_c[sub,3]==0,1,0)
    # num_sub = Isub%*%d_need
    num_sub = Isub%*%c(data_c[sub,1:2]%*%beta1_mle,data_c[sub,1:2]%*%beta2_mle,1)
    num = exp(num_sub)*sum(sum_num)
    
    prob = num/sum(sum_den)
    prop_ps = rbind(prop_ps,prob)  
  }
}
ps_suffi2 = prop_ps[-1,]
ps_suffi3 = 1-ps_suffi-ps_suffi2
thirdg = apply(as.matrix(cbind(data$A,datanew$A)),1,setdiff,x=c(0,1,2))
groups = cbind(data$A,datanew$A,thirdg)
ordersind = t(apply(groups,1,order))
probs = cbind(ps_suffi,ps_suffi2,ps_suffi3)

prob012 = matrix(-99,nrow=dim(probs)[1],ncol=3)
for (i in 1:dim(probs)[1]){
  prob012[i,] = probs[i,ordersind[i,]]
}


