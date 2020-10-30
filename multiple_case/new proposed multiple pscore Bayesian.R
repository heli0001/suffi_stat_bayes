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
  effec_true = c(2,4,2) #effect 01,02,12
  effec_simu = c(mean(data$y1-data$y0),mean(data$y2-data$y0),mean(data$y2-data$y1))
  data_obs = data[,c(1:5,8:10,14)]
  print(dataid)
  return(list(data=data,data_obs = data_obs,m=m,n=n,N=N,u1=u1,u2=u2,suffi=suffi,pstrue=pstrue,effec_simu=effec_simu))
}
simudata = mclapply(1:3,data_gene,rhoxu1 = 0,rhoxu2 = 0,rhoyu1 = 0,rhoyu2 = 0,m = 500,mc.cores=1)

psfunc_naive = function(dataid){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  m = simudata[[dataid]]$m
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  u1 = simudata[[dataid]]$u1
  u2 = simudata[[dataid]]$u2
  suffi = simudata[[dataid]]$suffi
  y0_naive = mean(data[data$A==0,]$y_obs)
  y1_naive = mean(data[data$A==1,]$y_obs)
  y2_naive = mean(data[data$A==2,]$y_obs)
  effect_naive = c(y1_naive-y0_naive,y2_naive-y0_naive,y2_naive-y1_naive)
  print(dataid)
  return(effect_naive=effect_naive)
}

psfunc_fix = function(dataid,iter=100,burn=0){
   data = simudata[[dataid]]$data
   data_obs = simudata[[dataid]]$data_obs
   m = simudata[[dataid]]$m
   N = simudata[[dataid]]$N
   n = simudata[[dataid]]$n
   u1 = simudata[[dataid]]$u1
   u2 = simudata[[dataid]]$u2
   suffi = simudata[[dataid]]$suffi
   psfunc = function (x,para){ ## directly use the pscore parameters are 1s
     exp(x%*%para)
   }
   beta1_ps_pos = matrix(-99,nrow=iter,ncol=m+3) # intercept,x1,x2, m dummy variables for cluster
   beta2_ps_pos = matrix(-99,nrow=iter,ncol=m+3) # intercept,x1,x2, m dummy variables for cluster
   mean_prior = rep(1,m+3)
   beta1_ps_pos[1,] = c(rep(1,3),u1)
   beta2_ps_pos[1,] = c(rep(1,3),u2)
   ## dummy variable/indicator matrix 
   cluind = data_obs$clu
   Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
   Istar[cbind(seq_along(cluind), cluind)] = 1
   data_ps = as.matrix(cbind(data_obs[,3:5],Istar))
   for (g in 2:iter){
     ### update beta1
     c1 = log(1+exp(data_ps%*%beta2_ps_pos[g-1,]))
     wi1 = apply(as.matrix(data_ps%*%beta1_ps_pos[g-1,]-c1),1,rpg,n=1,h=1) 
     zi1 = ifelse(data_obs$A==1,1,0)
     ki1 = zi1-1/2
     omega1 = diag(wi1)
     var_ps1 = solve(t(data_ps)%*%omega1%*%data_ps + 1/5*diag(m+3))
     mu_ps1 = var_ps1 %*%(t(data_ps)%*%(ki1+omega1%*%c1)+1/5*diag(m+3)%*%mean_prior)
     #mu_ps1 = var_ps %*%(t(data_ps)%*%kij+1/100*diag(m+3)%*%mean_prior)
     beta1_ps_pos[g,] = mvrnorm(1,mu_ps1,var_ps1)
     
     ### update beta2
     c2 = log(1+exp(data_ps%*%beta1_ps_pos[g,]))
     wi2 = apply(as.matrix(data_ps%*%beta2_ps_pos[g-1,]-c2),1,rpg,n=1,h=1) 
     zi2 = ifelse(data_obs$A==2,1,0)
     ki2 = zi2-1/2
     omega2 = diag(wi2)
     var_ps2 = solve(t(data_ps)%*%omega2%*%data_ps + 1/5*diag(m+3))
     mu_ps2 = var_ps2 %*%(t(data_ps)%*%(ki2+omega2%*%c2)+1/5*diag(m+3)%*%mean_prior)
     beta2_ps_pos[g,] = mvrnorm(1,mu_ps2,var_ps2)
   }
   beta1_ps_pos_mean = apply(beta1_ps_pos[(burn+1):iter,],2,mean)
   beta2_ps_pos_mean = apply(beta2_ps_pos[(burn+1):iter,],2,mean)
   w_fix_ps1 = apply(data_ps,1,psfunc,para=beta1_ps_pos_mean)/(1+apply(data_ps,1,psfunc,para=beta1_ps_pos_mean)+apply(data_ps,1,psfunc,para=beta2_ps_pos_mean))
   w_fix_ps2 = apply(data_ps,1,psfunc,para=beta2_ps_pos_mean)/(1+apply(data_ps,1,psfunc,para=beta1_ps_pos_mean)+apply(data_ps,1,psfunc,para=beta2_ps_pos_mean))
   w_fix_ps0 = 1-w_fix_ps1-w_fix_ps2
   data$w_fix = ifelse(data$A==0,w_fix_ps0,ifelse(data$A==1,w_fix_ps1,w_fix_ps2))
   data0 = data[data$A==0,]
   data1 = data[data$A==1,]
   data2 = data[data$A==2,]
   y0_fix = sum(data0$y_obs/data0$w_fix)/sum(1/data0$w_fix)
   y1_fix = sum(data1$y_obs/data1$w_fix)/sum(1/data1$w_fix)
   y2_fix = sum(data2$y_obs/data2$w_fix)/sum(1/data2$w_fix)
   effect_fix = c(y1_fix-y0_fix,y2_fix-y0_fix,y2_fix-y1_fix)
   print(dataid)
   return(list(data=data,beta1_ps_pos=beta1_ps_pos,beta2_ps_pos=beta2_ps_pos,effect_fix=effect_fix))
 }
fixeffet = mclapply(1:10,psfunc_fix,iter=10,burn=0,mc.cores=1)

psfunc_ran = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  m = simudata[[dataid]]$m
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  u1 = simudata[[dataid]]$u1
  u2 = simudata[[dataid]]$u2
  suffi = simudata[[dataid]]$suffi
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta1_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  alpha1_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var1_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha1
  
  beta2_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  alpha2_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var2_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha2
  
  beta1_ps_pos[1,] = rep(1,3)
  alpha1_ps_pos[1,] = u1
  var1_alpha_pos[1] = 1
  
  beta2_ps_pos[1,] = rep(1,3)
  alpha2_ps_pos[1,] = u2
  var2_alpha_pos[1] = 1
  
  data_ps = as.matrix(data_obs[,3:5])
  ## dummy variable/indicator matrix 
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (g in 2:iter){
    ## sample alpha1_ps
    for (cl in 1:m){
      data_j = as.matrix(data_obs[data_obs$clu==cl,3:5])
      ni = dim(data_j)[1]
      alpha_cl2 = rep(alpha2_ps_pos[g-1,cl],n[cl])
      c1 = log(1+exp(data_j%*%beta2_ps_pos[g-1,]+alpha_cl2))
      waij1 = apply(data_j%*%beta1_ps_pos[g-1,]+alpha1_ps_pos[g-1,cl]-c1,1,rpg,n=1,h=1)
      z1 = ifelse(data_obs[data_obs$clu==cl,]$A==1,1,0)
      kaij1 = z1-1/2
      La1 = kaij1/waij1
      sig1 = diag(waij1)
      var_a1 = 1/(1/var1_alpha_pos[g-1]+rep(1,ni)%*%sig1%*%c(rep(1,ni)))
      mu_a1 = var_a1*(rep(1,ni)%*%sig1%*%(La1-data_j%*%beta1_ps_pos[g-1,]+c1))
      alpha1_ps_pos[g,cl] = rnorm(1,mu_a1,sqrt(var_a1))
    }
    ## sample alpha2_ps
    for (cl in 1:m){
      data_j = as.matrix(data_obs[data_obs$clu==cl,3:5])
      ni = dim(data_j)[1]
      alpha_cl1 = rep(alpha1_ps_pos[g,cl],n[cl])
      c2 = log(1+exp(data_j%*%beta1_ps_pos[g-1,]+alpha_cl1))
      waij2 = apply(data_j%*%beta2_ps_pos[g-1,]+alpha2_ps_pos[g-1,cl]-c2,1,rpg,n=1,h=1)
      z2 = ifelse(data_obs[data_obs$clu==cl,]$A==2,1,0)
      kaij2 = z2-1/2
      La2 = kaij2/waij2
      sig2 = diag(waij2)
      var_a2 = 1/(1/var2_alpha_pos[g-1]+rep(1,ni)%*%sig2%*%c(rep(1,ni)))
      mu_a2 = var_a2*(rep(1,ni)%*%sig2%*%(La2-data_j%*%beta2_ps_pos[g-1,]+c2))
      alpha2_ps_pos[g,cl] = rnorm(1,mu_a2,sqrt(var_a2))
    }
    
    
    ## sample beta1_ps
    cb1 = log(1+exp(data_ps%*%beta2_ps_pos[g-1,]+Istar%*%alpha2_ps_pos[g,]))
    wij1 = apply(data_ps%*%beta1_ps_pos[g-1,]+Istar%*%alpha1_ps_pos[g,]-cb1,1,rpg,n=1,h=1)
    zb1 = ifelse(data_obs$A==1,1,0)
    kij1 = zb1-1/2
    L1 = kij1/wij1
    omega1 = diag(wij1)
    var_ps1 = solve(t(data_ps)%*%omega1%*%data_ps + 1/100*diag(3))
    mu_ps1 = var_ps1 %*%(t(data_ps)%*%omega1%*%(L1-Istar%*%alpha1_ps_pos[g,]+cb1))
    beta1_ps_pos[g,] = mvrnorm(1,mu_ps1,var_ps1)
    
    ## sample beta2_ps
    cb2 = log(1+exp(data_ps%*%beta1_ps_pos[g,]+Istar%*%alpha1_ps_pos[g,]))
    wij2 = apply(data_ps%*%beta2_ps_pos[g-1,]+Istar%*%alpha2_ps_pos[g,]-cb2,1,rpg,n=1,h=1)
    zb2 = ifelse(data_obs$A==2,1,0)
    kij2 = zb2-1/2
    L2 = kij2/wij2
    omega2 = diag(wij2)
    var_ps2 = solve(t(data_ps)%*%omega2%*%data_ps + 1/100*diag(3))
    mu_ps2 = var_ps2 %*%(t(data_ps)%*%omega2%*%(L2-Istar%*%alpha2_ps_pos[g,]+cb2))
    beta2_ps_pos[g,] = mvrnorm(1,mu_ps2,var_ps2)
    
    
    ## sample var1_alpha
    ralpha1 = sum((alpha1_ps_pos[g,])^2)
    shape1 = (0.002+m)/2
    scale1 = (0.002*1 + ralpha1)/2
    var1_alpha_pos[g] = rinvgamma(1,shape=shape1,rate=scale1)
    
    ## sample var2_alpha
    ralpha2 = sum((alpha2_ps_pos[g,])^2)
    shape2 = (0.002+m)/2
    scale2 = (0.002*1 + ralpha2)/2
    var2_alpha_pos[g] = rinvgamma(1,shape=shape2,rate=scale2)
    
    print(g)
  }
  beta1_ps_pos_mean = apply(beta1_ps_pos[(burn+1):iter,],2,mean)
  beta2_ps_pos_mean = apply(beta2_ps_pos[(burn+1):iter,],2,mean)
  alpha1_ps_pos_mean = apply(alpha1_ps_pos[(burn+1):iter,],2,mean)
  alpha2_ps_pos_mean = apply(alpha2_ps_pos[(burn+1):iter,],2,mean)
  alphaij1_mean = Istar%*%alpha1_ps_pos_mean
  alphaij2_mean = Istar%*%alpha2_ps_pos_mean
  psfunc_ran = function(x,para1){
    exp(x[1:3]%*%para1+x[4])
  }
  data_ps2 = as.matrix(cbind(data_ps,alphaij1_mean))
  data_ps3 = as.matrix(cbind(data_ps,alphaij2_mean))
  w_ran1 = apply(data_ps2,1,psfunc_ran,para1=beta1_ps_pos_mean)/(1+apply(data_ps2,1,psfunc_ran,para1=beta1_ps_pos_mean)+apply(data_ps3,1,psfunc_ran,para1=beta2_ps_pos_mean))
  w_ran2 = apply(data_ps3,1,psfunc_ran,para1=beta2_ps_pos_mean)/(1+apply(data_ps2,1,psfunc_ran,para1=beta1_ps_pos_mean)+apply(data_ps3,1,psfunc_ran,para1=beta2_ps_pos_mean))
  w_ran0 = 1-w_ran1-w_ran2
  data$w_ran = ifelse(data$A==0,w_ran0,ifelse(data$A==1,w_ran1,w_ran2))
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  y0_ran = sum(data0$y_obs/data0$w_ran)/sum(1/data0$w_ran)
  y1_ran = sum(data1$y_obs/data1$w_ran)/sum(1/data1$w_ran)
  y2_ran = sum(data2$y_obs/data2$w_ran)/sum(1/data2$w_ran)
  effect_ran = c(y1_ran-y0_ran,y2_ran-y0_ran,y2_ran-y1_ran)
  print(dataid)
  return(list(data=data,beta1_ps_pos=beta1_ps_pos,beta2_ps_pos=beta2_ps_pos,alpha1_ps_pos=alpha1_ps_pos,alpha2_ps_pos=alpha2_ps_pos,var1_alpha_pos=var1_alpha_pos,var2_alpha_pos=var2_alpha_pos,effect_ran=effect_ran))
}
raneffect = mclapply(1:10,psfunc_ran,iter=10,burn=0,mc.cores=1)

psfunc_suffi = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  m = simudata[[dataid]]$m
  N = simudata[[dataid]]$N
  n = simudata[[dataid]]$n
  u1 = simudata[[dataid]]$u1
  u2 = simudata[[dataid]]$u2
  suffi = simudata[[dataid]]$suffi
  suffi_method = function(beta1_mle,beta2_mle,data){
    prop_ps = 0
    for (cl in 1:m){
      data_c = as.matrix(data[data$clu==cl,c(3,4:5,8)]) ## intercept, x1,x2,A
      data_den = data_c[,1:3]
      Vbar = permuteGeneral(c(0,1,2), freqs = c(n[cl]-sum(suffi[cl,]),suffi[cl,1],suffi[cl,2]))
      nr = nrow(Vbar)
      aaa = lapply(split(Vbar,rep(1:ceiling(nr/10000), each=10000, length.out=nr)),matrix, ncol=n[cl])
      ########### calculate denominator
      den_chuncks = function(dataid){
        sum_den = rep(-99,dim(aaa[[dataid]])[1])
        sum_abar = rep(-99,n[cl])
        for (k in 1:dim(aaa[[dataid]])[1]){
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
        print(dataid)
       return (sum_den) 
      }
      tic()
      bbb = mclapply(1:length(aaa),den_chuncks,mc.cores = 1)
      toc()
      denom = sum(unlist(bbb))
      for (sub in 1:n[cl]){
        data_num = data_c[-sub,1:3]
        # allposs = permutations(n=3,r=n[cl],v=c(0,1,2),repeats.allowed=T) couldn't save all
        # Vbar = allposs
        # Vij = allposs[allposs[,sub]==data_c[sub,4],]
        # keep.crit = rep(-99,dim(Vij)[1])
        # keep.crit2 = rep(-99,dim(Vbar)[1])
        # for (i in 1:dim(Vij)[1]){
        #   t = Vij[i,]
        #   keep.crit[i] = ifelse(length(t[t==1])==suffi[cl,1] & length(t[t==2])==suffi[cl,2],"keep","discard")
        # }
        # ind = which(keep.crit=="keep")
        # Vij = Vij[ind,]
        # for (i in 1:dim(Vbar)[1]){
        #   t = Vbar[i,]
        #   keep.crit2[i] = ifelse(length(t[t==1])==suffi[cl,1] & length(t[t==2])==suffi[cl,2],"keep","discard")
        # }
        # ind2 = which(keep.crit2=="keep")
        # Vbar = Vbar[ind2,]
        Vij = Vbar[Vbar[,sub]==data_c[sub,4],]
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
        Isub[1] = ifelse(data_c[sub,4]==1,1,0)
        Isub[2] = ifelse(data_c[sub,4]==2,1,0)
        Isub[3] = ifelse(data_c[sub,4]==0,1,0)
        # num_sub = Isub%*%d_need
        num_sub = Isub%*%c(data_c[sub,1:3]%*%beta1_mle,data_c[sub,1:3]%*%beta2_mle,1)
        num = exp(num_sub)*sum(sum_num)
        
        prob = num/denom
        prop_ps = rbind(prop_ps,prob)  
      }
    }
    ps_suffi = prop_ps[-1,]
    return (ps_suffi)
  }
  #full joint log posterior of beta|y for the logistic regression model:
  posteriorfunc=function(param1,param2,prior.mn=0,prior.sd=1){
    psall1_suffi = suffi_method(param1,param2,data)
    # ps_suffi = ifelse(data_obs$A==1,psall1_suffi,1-psall1_suffi)
    # like = 0
    # for (i in 1:N){
    #   like=like+dbinom(data_obs$A[i],size=1,prob=ps_suffi[i],log=TRUE) 
    # }
    like = sum(dbinom(1,1,psall1_suffi,log=TRUE))
    prior1 = sum(dnorm(param1,prior.mn,prior.sd,log=TRUE))
    prior2 = sum(dnorm(param2,prior.mn,prior.sd,log=TRUE))
    return(like+prior1+prior2)
  }
  
  # Start sampling
  #Initial values:
  beta1 = rnorm(3,0,1)
  beta2 = rnorm(3,0,1)
  # beta = c(1,1,1)
  beta1_ps_pos = matrix(-99,iter,3)
  beta2_ps_pos = matrix(-99,iter,3)
  acc1 = rep(0,3)
  acc2 = rep(0,3)
  cur_log_post = posteriorfunc(beta1,beta2,prior.mn=0,prior.sd=10)
  for(g in 1:iter){
    #Update beta1 using MH sampling:
    for(j in 1:3){
      # Draw candidate:
      canbeta1 = beta1
      canbeta1[j] = rnorm(1,beta1[j],1)
      can1_log_post = posteriorfunc(canbeta1,beta2,prior.mn=0,prior.sd=10)
      # Compute acceptance ratio:
      R = exp(can1_log_post-cur_log_post)
      if(runif(1)<R){
        beta1 = canbeta1
        cur_log_post = can1_log_post
        acc1[j] = acc1[j]+1
      }
    }
    #Update beta2 using MH sampling:
    for (j in 1:3){
      canbeta2 = beta2
      canbeta2[j] = rnorm(1,beta2[j],1)
      can2_log_post = posteriorfunc(beta1,canbeta2,prior.mn=0,prior.sd=10)
      # Compute acceptance ratio:
      R = exp(can2_log_post-cur_log_post)
      if(runif(1)<R){
        beta2 = canbeta2
        cur_log_post = can2_log_post
        acc2[j] = acc2[j]+1
      }
    }
    print(g)
    print(acc1)
    print(acc2)
    beta1_ps_pos[g,]=beta1
    beta2_ps_pos[g,]=beta2
  }
  beta1_ps_pos_mean = apply(beta1_ps_pos[(burn+1):iter,],2,mean)
  beta2_ps_pos_mean = apply(beta2_ps_pos[(burn+1):iter,],2,mean)
  ps_suffi_mean = suffi_method(beta1_ps_pos_mean,beta2_ps_pos_mean,data)
  data$ps_weight = ps_suffi_mean
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  ybar_0 = sum(data0$y_obs/data0$ps_weight)/sum(1/data0$ps_weight)
  ybar_1 = sum(data1$y_obs/data1$ps_weight)/sum(1/data1$ps_weight)
  ybar_2 = sum(data2$y_obs/data2$ps_weight)/sum(1/data2$ps_weight)
  effect_suffi = c(ybar_1-ybar_0,ybar_2-ybar_0,ybar_2-ybar_1)
  print(dataid)
  return(data=data,beta1_ps_pos=beta1_ps_pos,beta2_ps_pos=beta2_ps_pos,effect_suffi=effect_suffi)
}
suffieffect = mclapply(1:1,psfunc_suffi,iter=100,burn=0,mc.cores=1)

effect_fix = function(dataid,fixeffect0000,iter=1000,burn=100){
  data = simudata0000[[dataid]]$data
  data_obs = simudata0000[[dataid]]$data_obs
  N = simudata0000[[dataid]]$N
  u1 = simudata0000[[dataid]]$u1
  u2 = simudata0000[[dataid]]$u2
  suffi = simudata0000[[dataid]]$suffi
  m = simudata0000[[dataid]]$m
  n = simudata0000[[dataid]]$n
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)
  }
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  data_ps = as.matrix(cbind(data_obs[,3:5],Istar))
  beta1_ps_pos = fixeffect0000[[dataid]]$beta1_ps_pos
  beta2_ps_pos = fixeffect0000[[dataid]]$beta2_ps_pos
  #########################################
  beta1_ps_pos_mean = apply(beta1_ps_pos[(burn+1):iter,],2,mean)
  beta2_ps_pos_mean = apply(beta2_ps_pos[(burn+1):iter,],2,mean)
  w_fix_ps1 = apply(data_ps,1,psfunc,para=beta1_ps_pos_mean)/(1+apply(data_ps,1,psfunc,para=beta1_ps_pos_mean)+apply(data_ps,1,psfunc,para=beta2_ps_pos_mean))
  w_fix_ps2 = apply(data_ps,1,psfunc,para=beta2_ps_pos_mean)/(1+apply(data_ps,1,psfunc,para=beta1_ps_pos_mean)+apply(data_ps,1,psfunc,para=beta2_ps_pos_mean))
  w_fix_ps0 = 1-w_fix_ps1-w_fix_ps2
  data$w_fix = ifelse(data$A==0,w_fix_ps0,ifelse(data$A==1,w_fix_ps1,w_fix_ps2))
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  y0_fix = sum(data0$y_obs/data0$w_fix)/sum(1/data0$w_fix)
  y1_fix = sum(data1$y_obs/data1$w_fix)/sum(1/data1$w_fix)
  y2_fix = sum(data2$y_obs/data2$w_fix)/sum(1/data2$w_fix)
  effect_fix = c(y1_fix-y0_fix,y2_fix-y0_fix,y2_fix-y1_fix)
  print(dataid)
  return(effect_fix=effect_fix)
}
fix_effect = mclapply(1:10,effect_fix,fixeffect0000,iter=1000,burn=100)
fff = apply(matrix(unlist(fix_effect),byrow=T,ncol=3),2,mean)

effect_ran = function(dataid,raneffect0000,iter=1000,burn=100){
  data = simudata0000[[dataid]]$data
  data_obs = simudata0000[[dataid]]$data_obs
  N = simudata0000[[dataid]]$N
  u1 = simudata0000[[dataid]]$u1
  u2 = simudata0000[[dataid]]$u2
  suffi = simudata0000[[dataid]]$suffi
  m = simudata0000[[dataid]]$m
  n = simudata0000[[dataid]]$n
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  data_ps = as.matrix(data_obs[,3:5])
  ## dummy variable/indicator matrix 
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  beta1_ps_pos =  raneffect0000[[dataid]]$beta1_ps_pos
  beta2_ps_pos =  raneffect0000[[dataid]]$beta2_ps_pos
  alpha1_ps_pos =  raneffect0000[[dataid]]$alpha1_ps_pos
  alpha2_ps_pos =  raneffect0000[[dataid]]$alpha2_ps_pos
  ##############################################################
  beta1_ps_pos_mean = apply(beta1_ps_pos[(burn+1):iter,],2,mean)
  beta2_ps_pos_mean = apply(beta2_ps_pos[(burn+1):iter,],2,mean)
  alpha1_ps_pos_mean = apply(alpha1_ps_pos[(burn+1):iter,],2,mean)
  alpha2_ps_pos_mean = apply(alpha2_ps_pos[(burn+1):iter,],2,mean)
  alphaij1_mean = Istar%*%alpha1_ps_pos_mean
  alphaij2_mean = Istar%*%alpha2_ps_pos_mean
  psfunc_ran = function(x,para1){
    exp(x[1:3]%*%para1+x[4])
  }
  data_ps2 = as.matrix(cbind(data_ps,alphaij1_mean))
  data_ps3 = as.matrix(cbind(data_ps,alphaij2_mean))
  w_ran1 = apply(data_ps2,1,psfunc_ran,para1=beta1_ps_pos_mean)/(1+apply(data_ps2,1,psfunc_ran,para1=beta1_ps_pos_mean)+apply(data_ps3,1,psfunc_ran,para1=beta2_ps_pos_mean))
  w_ran2 = apply(data_ps3,1,psfunc_ran,para1=beta2_ps_pos_mean)/(1+apply(data_ps2,1,psfunc_ran,para1=beta1_ps_pos_mean)+apply(data_ps3,1,psfunc_ran,para1=beta2_ps_pos_mean))
  w_ran0 = 1-w_ran1-w_ran2
  data$w_ran = ifelse(data$A==0,w_ran0,ifelse(data$A==1,w_ran1,w_ran2))
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  y0_ran = sum(data0$y_obs/data0$w_ran)/sum(1/data0$w_ran)
  y1_ran = sum(data1$y_obs/data1$w_ran)/sum(1/data1$w_ran)
  y2_ran = sum(data2$y_obs/data2$w_ran)/sum(1/data2$w_ran)
  effect_ran = c(y1_ran-y0_ran,y2_ran-y0_ran,y2_ran-y1_ran)
  print(dataid)
  return(effect_ran=effect_ran)
}
ran_effect = mclapply(1:10,effect_ran,raneffect0000,iter=1000,burn=100)
rrr = apply(matrix(unlist(ran_effect),byrow=T,ncol=3),2,mean)

effect_suffi = function(dataid,suffieffect0000,iter=1000,burn=100){
  data = suffieffect0000[[dataid]]$data
  m = simudata0000[[dataid]]$m
  n = simudata0000[[dataid]]$n
  suffi = simudata0000[[dataid]]$suffi
  beta1_ps_pos = suffieffect0000[[dataid]]$beta1_ps_pos
  beta2_ps_pos = suffieffect0000[[dataid]]$beta2_ps_pos
  suffi_method = function(beta1_mle,beta2_mle,data){
    prop_ps = 0
    for (cl in 1:m){
      data_c = as.matrix(data[data$clu==cl,c(3,4:5,8)]) ## intercept, x1,x2,A
      for (sub in 1:n[cl]){
        data_num = data_c[-sub,1:3]
        data_den = data_c[,1:3]
        allposs = permutations(n=3,r=n[cl],v=c(0,1,2),repeats.allowed=T)
        Vbar = allposs
        Vij = allposs[allposs[,sub]==data_c[sub,4],]
        keep.crit = rep(-99,dim(Vij)[1])
        keep.crit2 = rep(-99,dim(Vbar)[1])
        for (i in 1:dim(Vij)[1]){
          t = Vij[i,]
          keep.crit[i] = ifelse(length(t[t==1])==suffi[cl,1] & length(t[t==2])==suffi[cl,2],"keep","discard")
        }
        ind = which(keep.crit=="keep")
        Vij = Vij[ind,]
        for (i in 1:dim(Vbar)[1]){
          t = Vbar[i,]
          keep.crit2[i] = ifelse(length(t[t==1])==suffi[cl,1] & length(t[t==2])==suffi[cl,2],"keep","discard")
        }
        ind2 = which(keep.crit2=="keep")
        Vbar = Vbar[ind2,]
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
        Isub[1] = ifelse(data_c[sub,4]==1,1,0)
        Isub[2] = ifelse(data_c[sub,4]==2,1,0)
        Isub[3] = ifelse(data_c[sub,4]==0,1,0)
        # num_sub = Isub%*%d_need
        num_sub = Isub%*%c(data_c[sub,1:3]%*%beta1_mle,data_c[sub,1:3]%*%beta2_mle,1)
        num = exp(num_sub)*sum(sum_num)
        ########### calculate denominator
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
        prob = num/sum(sum_den)
        prop_ps = rbind(prop_ps,prob)  
      }
    }
    ps_suffi = prop_ps[-1,]
    return (ps_suffi)
  }
  ################################################################
  beta1_ps_pos_mean = apply(beta1_ps_pos[(burn+1):iter,],2,mean)
  beta2_ps_pos_mean = apply(beta2_ps_pos[(burn+1):iter,],2,mean)
  ps_suffi_mean = suffi_method(beta1_ps_pos_mean,beta2_ps_pos_mean,data)
  data$ps_weight = ps_suffi_mean
  data0 = data[data$A==0,]
  data1 = data[data$A==1,]
  data2 = data[data$A==2,]
  ybar_0 = sum(data0$y_obs/data0$ps_weight)/sum(1/data0$ps_weight)
  ybar_1 = sum(data1$y_obs/data1$ps_weight)/sum(1/data1$ps_weight)
  ybar_2 = sum(data2$y_obs/data2$ps_weight)/sum(1/data2$ps_weight)
  effect_suffi = c(ybar_1-ybar_0,ybar_2-ybar_0,ybar_2-ybar_1)
  print(dataid)
  return(effect_suffi=effect_suffi)
}
suffi_effect = mclapply(1:10,effect_suffi,suffieffect0000,iter=1000,burn=100)
sss = apply(matrix(unlist(suffi_effect),byrow=T,ncol=3),2,mean)

effect_naive =matrix(-99,nrow=10,ncol=3)
for (i in 1:10){
  data = simudata0000[[i]]$data
  y0_naive = mean(data[data$A==0,]$y_obs)
  y1_naive = mean(data[data$A==1,]$y_obs)
  y2_naive = mean(data[data$A==2,]$y_obs)
  effect_naive[i,] = c(y1_naive-y0_naive,y2_naive-y0_naive,y2_naive-y1_naive)
}
nnn = apply(effect_naive,2,mean)

effect_simu = matrix(-99,nrow=10,ncol=3)
for (i in 1:10){
  effect_simu[i,] = simudata0000[[i]]$effec_simu
}
eee = apply(effect_simu,2,mean)
