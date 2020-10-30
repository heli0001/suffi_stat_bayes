rm(list=ls())
library(MASS)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)
########################Part 1: Data generation #######################
data_gene = function(dataid,rhoxu = 0,rhoyu = 0,m = 500){
  n = as.matrix(sample(2:20,m,replace=T,prob=rep(1/19,19)))
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
simudata = lapply(c(1:5),data_gene,rhoxu = 0,rhoyu = 0,m = 20)

psfunc_fix = function(dataid,iter=100,burn=0){
   data = simudata[[dataid]]$data
   data_obs = simudata[[dataid]]$data_obs
   N = simudata[[dataid]]$N
   m = simudata[[dataid]]$m
   n = simudata[[dataid]]$n
   u = simudata[[dataid]]$u
   suff_T = simudata[[dataid]]$suff_T
   psfunc = function (x,para){ ## directly use the pscore parameters are 1s
     exp(x%*%para)/(1+exp(x%*%para))
   }
   beta_ps_pos = matrix(-99,nrow=iter,ncol=m+3) # intercept,x1,x2, m dummy variables for cluster
   mean_prior = rep(1,m+3)
   beta_ps_pos[1,] = c(rep(1,3),u)
   ## dummy variable/indicator matrix 
   cluind = data_obs$clu
   Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
   Istar[cbind(seq_along(cluind), cluind)] = 1
   data_ps = as.matrix(cbind(data_obs[,3:5],Istar))
   for (g in 2:iter){
     wij = apply(data_ps%*%beta_ps_pos[g-1,],1,rpg,n=1,h=1) 
     kij = data_obs$A-1/2
     omega = diag(wij)
     var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/100*diag(m+3))
     # mu_ps = var_ps %*%t(data_ps)%*%kij
     mu_ps = var_ps %*%(t(data_ps)%*%kij+1/100*diag(m+3)%*%mean_prior)
     beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
     print(g)
   }
   beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
   data$w_fix = apply(data_ps,1,psfunc,para=beta_ps_pos_mean)
   effec_fix = mean(data$A*data$y_obs/data$w_fix-(1-data$A)*data$y_obs/(1-data$w_fix))
   print(dataid)
   return(list(data=data,beta_ps_pos=beta_ps_pos,effec_fix=effec_fix))
 }
fixeffet = lapply(c(1:5),psfunc_fix,iter=20,burn=0)

psfunc_ran = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  N = simudata[[dataid]]$N
  m = simudata[[dataid]]$m
  n = simudata[[dataid]]$n
  u = simudata[[dataid]]$u
  suff_T = simudata[[dataid]]$suff_T
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=3) # intercept,x1,x2
  alpha_ps_pos = matrix(-99,nrow=iter,ncol=m)
  var_alpha_pos = rep(-99,iter) ## hyperparameter for random effect alpha
  beta_ps_pos[1,] = rep(1,3)
  alpha_ps_pos[1,] = u
  var_alpha_pos[1] = 1
  data_ps = as.matrix(data_obs[,3:5])
  ## dummy variable/indicator matrix 
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  for (g in 2:iter){
    ## sample alpha_ps
    for (cl in 1:m){
      data_j = as.matrix(data_obs[data_obs$clu==cl,3:5])
      ni = dim(data_j)[1]
      waij = apply(data_j%*%beta_ps_pos[g-1,]+alpha_ps_pos[g-1,cl],1,rpg,n=1,h=1)
      kaij = data_obs[data_obs$clu==cl,]$A-1/2
      La = kaij/waij
      sig = diag(waij)
      var_a = 1/(1/var_alpha_pos[g-1]+rep(1,ni)%*%sig%*%c(rep(1,ni)))
      mu_a = var_a*(rep(1,ni)%*%sig%*%(La-data_j%*%beta_ps_pos[g-1,]))
      alpha_ps_pos[g,cl] = rnorm(1,mu_a,sqrt(var_a))
    }
    ## sample beta_ps
    wij = apply(data_ps%*%beta_ps_pos[g-1,]+Istar%*%alpha_ps_pos[g,],1,rpg,n=1,h=1)
    kij = data_obs$A-1/2
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
  data$w_ran = apply(data_ps2,1,psfunc2,para1=beta_ps_pos_mean)
  effec_ran = mean(data$A*data$y_obs/data$w_ran-(1-data$A)*data$y_obs/(1-data$w_ran))
  print(dataid)
  return(list(data=data,beta_ps_pos=beta_ps_pos,alpha_ps_pos=alpha_ps_pos,var_alpha_pos=var_alpha_pos,effec_ran=effec_ran))
}
raneffect = lapply(c(1:5),psfunc_ran,iter=20,burn=0)

psfunc_suffi = function(dataid,iter=100,burn=0){
  data = simudata[[dataid]]$data
  data_obs = simudata[[dataid]]$data_obs
  N = simudata[[dataid]]$N
  m = simudata[[dataid]]$m
  n = simudata[[dataid]]$n
  u = simudata[[dataid]]$u
  suff_T = simudata[[dataid]]$suff_T
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
  posteriorfunc=function(param,prior.mn=0,prior.sd=1){
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
  
  # Start sampling
  #Initial values:
  # beta = rnorm(3,0,1)
  beta = c(1,1,1)
  beta_ps_pos = matrix(-99,iter,3)
  acc = rep(0,3)
  cur_log_post = posteriorfunc(beta,prior.mn=0,prior.sd=10)
  for(g in 1:iter){
    #Update beta using MH sampling:
    for(j in 1:3){
      # Draw candidate:
      canbeta = beta
      canbeta[j] = rnorm(1,beta[j],1)
      can_log_post = posteriorfunc(canbeta,prior.mn=0,prior.sd=10)
      # Compute acceptance ratio:
      R = exp(can_log_post-cur_log_post)
      if(runif(1)<R){
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
  ps_suffi_mean = suffi_method(beta_ps_pos_mean,data_obs,m=m,n=n,suff_T=suff_T)
  data$w_suffi = ifelse(data_obs$A==1,ps_suffi_mean,1-ps_suffi_mean)
  effec_suffi = mean(data$A*data$y_obs/data$w_suffi-(1-data$A)*data$y_obs/(1-data$w_suffi))
  print(dataid)
  return(list(data=data,beta_ps_pos=beta_ps_pos,effec_suffi=effec_suffi))
}
suffieffect = lapply(c(1:5),psfunc_suffi,iter=20,burn=0)
