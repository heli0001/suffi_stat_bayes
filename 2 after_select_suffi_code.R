ps_est = function(x,para){
  1/(1+exp(-x[1:2]%*%para-x[3]))
}


binary_suffi_func=function(cls,data,m,beta_est,suff_T,n){
    prop_ps = c()
    dem=rep(-99,m)
    ## denominator part
    data_clu = as.matrix(data[data$clu==cls,3:4])
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
      if (data[data$clu==cls,]$z[sub]==0){
        nume= rnum[suffT+1,ni+1]
      }else{
        nume= dem[cls]-rnum[suffT+1,ni+1]
      }
      prob = nume/dem[cls]
      prop_ps = c(prop_ps,prob)
    }
    return(prop_ps)
  }

iter=11000
burn=1000

effect_func = function(dataid,simdata,clu_res,iter=iter,burn=burn){
  ptm = proc.time()
  set.seed(dataid)
  data_obs = simdata[[dataid]]$data_obs
  m = simdata[[dataid]]$m
  n = simdata[[dataid]]$n
  suff_T = simdata[[dataid]]$suff_T
  gamma_pos = clu_res[[dataid]]$gamma_pos
  delta = clu_res[[dataid]]$delta
  zeta = clu_res[[dataid]]$zeta
  
  
  ## dummy variable/indicator matrix for cluster
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  
  
  ###### based on NMIG for significant clusters, use sufficient method to calculate pscore, 
  ###### others still use the logistic one
  gamma_est = apply(gamma_pos[(burn+1):iter,],2,mean)
  inclu_prob = apply(delta[(burn+1):iter,],2,mean)
  
  # nonsig clusters
  non_inclu = which(inclu_prob<0.5)
  data_ps = cbind(data_obs,Istar)
  data_nsig = data_ps[data_ps$clu %in% non_inclu,]
  zeta_est =  apply(zeta[(burn+1):iter,],2,mean)[non_inclu]
  data_need = cbind(data_nsig$x1,data_nsig$x2,as.matrix(data_nsig[,-c(1:7)][,non_inclu])%*%zeta_est)
  data_nsig$ps_est = apply(data_need,1,ps_est,para=gamma_est)
  effect1 = data_nsig$z*data_nsig$y_obs/data_nsig$ps_est-(1-data_nsig$z)*data_nsig$y_obs/(1-data_nsig$ps_est)
  
  # sig clusters
  inclu = which(inclu_prob>=0.5)
  data_sig = data_obs[data_obs$clu %in% inclu,]
  data_sig$ps_est = unlist(lapply(inclu,binary_suffi_func,data_obs,m,gamma_est,suff_T,n))
  effect2 = data_sig$z*data_sig$y_obs/data_sig$ps_est-(1-data_sig$z)*data_sig$y_obs/(data_sig$ps_est)
  
  #################### calculate using IPW ###########################
  effect_nmig_suffi = mean(c(effect1,effect2))  
  print(dataid)
  time_span = proc.time() - ptm
  return(list(inclu=inclu,effect_nmig_suffi=effect_nmig_suffi,time_span=time_span))
}
effect_res = mclapply(1:2,effect_func,simdatam20,clu_res,iter=iter,burn=burn,mc.cores = 1)
