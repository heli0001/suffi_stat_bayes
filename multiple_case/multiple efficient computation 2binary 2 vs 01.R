tic()
############# Run based on 2 vs 01 ##################
data_new1=data
data_new1$A_relabel1=ifelse(data_new1$A==2,1,0) ## group 2=1,others(0&1)=0
suffi_new = suffi[,2] ## number of group 0 within each cluster
data_new1$suffi_new = rep(suffi_new,n)

############# sufficient statistic based estimation
library(survival)
clog_summ = clogit(A_relabel1~x1+x2+strata(clu), data=data_new1) ## conditional logistic MLE
beta_est = clog_summ$coefficients

# logs1 = glm(A_relabel1~1+x1+x2,family = "binomial",data=data_new1)
# beta_est=coef(logs1)
dem=rep(-99,m)
prop_ps = 0
for (cls in 1:m){
  ## denominator part
  data_clu = as.matrix(data_new1[data_new1$clu==cls,4:5])
  ni=n[cls]
  suffT =suffi_new[cls]
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
    if (data_new1[data_new1$clu==cls,]$A_relabel1[sub]==0){
      nume= rnum[suffT+1,ni+1]
    }else{
      nume= dem[cls]-rnum[suffT+1,ni+1]
    }
    
    prob = nume/dem[cls]
    prop_ps = rbind(prop_ps,prob)
  }
}
ps_suffi_A=prop_ps[-1]

############## rerun it again only based on 0&1 groups ############
data_bc = data[data$A==1|data$A==0,]
data_bc$A_relabel2=ifelse(data_bc$A==1,1,0) ## group 1=1,group0=0
suffi_new2 = rep(-99,m) ## number of group 1 within each cluster
for (i in 1:m){
  suffi_new2[i] = sum(data_bc[data_bc$clu==i,]$A_relabel2)
}
n_new = table(data_bc$clu)
data_bc$suffi_new2 = rep(suffi_new2,n_new)

############# sufficient statistic based estimation 
clog_summ2 = clogit(A_relabel2~x1+x2+strata(clu), data=data_bc) ## conditional logistic MLE
beta_est2 = clog_summ2$coefficients
#logs2 = glm(A_relabel2~1+x1+x2,family = "binomial",data=data_bc)
#beta_est2=coef(logs2)
dem2=rep(-99,m)
prop_ps2 = 0
for (cls in 1:m){
  ## denominator part
  data_clu = as.matrix(data_bc[data_bc$clu==cls,4:5])
  ni=n_new[cls]
  suffT =suffi_new2[cls]
  rdenom = matrix(-99,nrow=suffT+1,ncol=ni+1)
  rdenom[lower.tri(rdenom)]=0 ## lower part =0
  rdenom[1,]=1
  for (rowid in 2:(suffT+1)){
    for (colid in 2:(ni+1)){
      rdenom[rowid,colid]=rdenom[rowid,colid-1] + exp(data_clu[colid-1,]%*%beta_est2)*rdenom[rowid-1,colid-1]
    }
  }
  dem2[cls] = rdenom[suffT+1,ni+1]
  for (sub in 1:n_new[cls]){
    ## numerator part
    rnum = matrix(-99,nrow=suffT+1,ncol=ni+1)
    rnum[lower.tri(rnum)]=0 ## lower part =0
    rnum[1,]=1
    for (rowid in 2:(suffT+1)){
      for (colid in 2:(ni+1)){
        rnum[rowid,colid]=rnum[rowid,colid-1] + exp(data_clu[colid-1,]%*%beta_est2)*(1-as.numeric(colid-1==sub))*rnum[rowid-1,colid-1]
      }
    }
    if (data_bc[data_bc$clu==cls,]$A_relabel2[sub]==0){
      nume= rnum[suffT+1,ni+1]
    }else{
      nume= dem2[cls]-rnum[suffT+1,ni+1]
    }
    
    prob = nume/dem2[cls]
    prop_ps2 = rbind(prop_ps2,prob)
  }
}
ps_suffi_BC=prop_ps2[-1]
index_bc = which(data_new1$A!=2) ## index of units in group 0 and 1
ps_suffi4=ps_suffi_A
ps_suffi4[index_bc] = ps_suffi_A[index_bc]*ps_suffi_BC ## update pscore for group1 and 2
toc()
plot(ps_suffi,ps_suffi4)

aaa = cbind(mean((ps_suffi-ps_suffi2)^2),mean((ps_suffi-ps_suffi3)^2),mean((ps_suffi-ps_suffi4)^2))
colnames(aaa) = c("0 vs 12","1 vs 02","2 vs 01")

summ_x1 = rbind(summary(data[data$A==0,]$x1),
summary(data[data$A==1,]$x1),
summary(data[data$A==2,]$x1))

summ_x2 = rbind(summary(data[data$A==0,]$x2),
summary(data[data$A==1,]$x2),
summary(data[data$A==2,]$x2))