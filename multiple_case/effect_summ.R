load("simudata05.RData")
library(MASS)
library(parallel)
library(hier.part)
library(invgamma)
library(truncnorm)
library(BayesLogit)
library(mvtnorm)

iter=10
burn=0
dataid=100

psfunc_fix = function(dataid,iter=100,burn=0){
  data = simudata05[[dataid]]$data
  data_obs = simudata05[[dataid]]$data_obs
  N = simudata05[[dataid]]$N
  m = simudata05[[dataid]]$m
  n = simudata05[[dataid]]$n
  u = simudata05[[dataid]]$u
  suff_T = simudata05[[dataid]]$suff_T
  psfunc = function (x,para){ ## directly use the pscore parameters are 1s
    exp(x%*%para)/(1+exp(x%*%para))
  }
  beta_ps_pos = matrix(-99,nrow=iter,ncol=m+3) # intercept,x1,x2, m dummy variables for cluster
  # beta_ps_pos[1,] = rep(1,m+3)
  beta_ps_pos[1,] = c(rep(1,3),u)
  ## dummy variable/indicator matrix 
  cluind = data_obs$clu
  Istar = matrix(0, nrow=length(cluind), ncol=length(unique(cluind)))
  Istar[cbind(seq_along(cluind), cluind)] = 1
  # data_ps = fastDummies::dummy_cols(data_obs[,1])
  data_ps = as.matrix(cbind(data_obs[,3:5],Istar))
  for (g in 2:iter){
    wij = apply(data_ps%*%beta_ps_pos[g-1,],1,rpg,n=1,h=1) 
    kij = data_obs$A-1/2
    omega = diag(wij)
    var_ps = solve(t(data_ps)%*%omega%*%data_ps + 1/100*diag(m+3))
    mu_ps = var_ps %*%t(data_ps)%*%kij
    beta_ps_pos[g,] = mvrnorm(1,mu_ps,var_ps)
  }
  beta_ps_pos_mean = apply(beta_ps_pos[(burn+1):iter,],2,mean)
  data$w_fix = apply(data_ps,1,psfunc,para=beta_ps_pos_mean)
  effec_fix = mean(data$A*data$y_obs/data$w_fix-(1-data$A)*data$y_obs/(1-data$w_fix))
  print(dataid)
  return(list(data=data,beta_ps_pos=beta_ps_pos,effec_fix=effec_fix))
}
fixeffet05 = mclapply(1:dataid,psfunc_fix,iter=iter,burn=burn,mc.cores=20)
save(fixeffet05,file="effect05_under_fix.RData")

load("simudata50.RData")
load("effect50_under_fix.RData")
load("effect50_under_ran.RData")
load("effect50_under_suffi.RData")

data_gene = 10
effect_naive_est = rep(-99,data_gene)
effect_true = rep(-99,data_gene)
effect_fix_est = rep(-99,data_gene)
effect_ran_est = rep(-99,data_gene)
effect_suffi_est = rep(-99,data_gene)
for (i in 1:data_gene){
  effect_naive_est[i] = simudata50[[i]]$effec_naiv
  effect_true[i]=simudata50[[i]]$effec_simu
  effect_fix_est[i]=fixeffet50[[i]]$effec_fix
  effect_ran_est[i]=raneffet50[[i]]$effec_ran
  effect_suffi_est[i]=suffieffect50[[i]]$effec_suffi
}
sqrt(mean(effect_fix_est-effect_true)^2)
sqrt(mean(effect_ran_est-effect_true)^2)
sqrt(mean(effect_suffi_est-effect_true)^2)

simu_summ = cbind(mean(effect_true),sd(effect_true),0)
naive_summ = cbind(mean(effect_naive_est),sd(effect_naive_est),sqrt(mean(effect_naive_est-effect_true)^2))
fix_summ = cbind(mean(effect_fix_est),sd(effect_fix_est),sqrt(mean(effect_fix_est-effect_true)^2))
ran_summ = cbind(mean(effect_ran_est),sd(effect_ran_est),sqrt(mean(effect_ran_est-effect_true)^2))
suffi_summ = cbind(mean(effect_suffi_est),sd(effect_suffi_est),sqrt(mean(effect_suffi_est-effect_true)^2))

effect_summ50 = rbind(simu_summ,naive_summ,fix_summ,ran_summ,suffi_summ)
colnames(effect_summ50) = c("mean","sd","RSME")
rownames(effect_summ50) = c("simu","naive","fix","ran","suffi")

load("simudata00.RData")
load("effect_under_fix.RData")
load("effect_under_ran.RData")
load("effect_under_suffi.RData")
data_gene = 10
effect_naive_est = rep(-99,data_gene)
effect_true = rep(-99,data_gene)
effect_fix_est = rep(-99,data_gene)
effect_ran_est = rep(-99,data_gene)
effect_suffi_est = rep(-99,data_gene)
for (i in 1:data_gene){
  effect_naive_est[i] = simudata00[[i]]$effec_naiv
  effect_true[i]=simudata00[[i]]$effec_simu
  effect_fix_est[i]=fixeffet[[i]]$effec_fix
  effect_ran_est[i]=raneffet[[i]]$effec_ran
  effect_suffi_est[i]=suffieffect[[i]]$effec_suffi
}
sqrt(mean(effect_fix_est-effect_true)^2)
sqrt(mean(effect_ran_est-effect_true)^2)
sqrt(mean(effect_suffi_est-effect_true)^2)

simu_summ = cbind(mean(effect_true),sd(effect_true),0)
naive_summ = cbind(mean(effect_naive_est),sd(effect_naive_est),sqrt(mean(effect_naive_est-effect_true)^2))
fix_summ = cbind(mean(effect_fix_est),sd(effect_fix_est),sqrt(mean(effect_fix_est-effect_true)^2))
ran_summ = cbind(mean(effect_ran_est),sd(effect_ran_est),sqrt(mean(effect_ran_est-effect_true)^2))
suffi_summ = cbind(mean(effect_suffi_est),sd(effect_suffi_est),sqrt(mean(effect_suffi_est-effect_true)^2))

effect_summ00 = rbind(simu_summ,naive_summ,fix_summ,ran_summ,suffi_summ)
colnames(effect_summ00) = c("mean","sd","RSME")
rownames(effect_summ00) = c("simu","naive","fix","ran","suffi")

############### fixed effect ####################
load("effectfix00.RData")
load("effectfix05.RData")
load("effectfix50.RData")
load("effectfix55.RData")
load("simudata00.RData")
load("simudata05.RData")
load("simudata50.RData")
load("simudata55.RData")
data_gene = 10
effect_fix_est00 = rep(-99,data_gene)
effect_fix_est05 = rep(-99,data_gene)
effect_fix_est50 = rep(-99,data_gene)
effect_fix_est55 = rep(-99,data_gene)
effect_true00 = rep(-99,data_gene)
effect_true05 = rep(-99,data_gene)
effect_true50 = rep(-99,data_gene)
effect_true55 = rep(-99,data_gene)

for (i in 1:data_gene){
  effect_true00[i]=simudata00[[i]]$effec_simu
  effect_true05[i]=simudata05[[i]]$effec_simu
  effect_true50[i]=simudata50[[i]]$effec_simu
  effect_true55[i]=simudata55[[i]]$effec_simu
  
  effect_fix_est00[i]=fixeffet00[[i]]$effec_fix
  effect_fix_est05[i]=fixeffet05[[i]]$effec_fix
  effect_fix_est50[i]=fixeffet50[[i]]$effec_fix
  effect_fix_est55[i]=fixeffet55[[i]]$effec_fix
}
fix00_summ = cbind(mean(effect_fix_est00),sd(effect_fix_est00),sqrt(mean(effect_fix_est00-effect_true00)^2))
fix05_summ = cbind(mean(effect_fix_est05),sd(effect_fix_est05),sqrt(mean(effect_fix_est05-effect_true05)^2))
fix50_summ = cbind(mean(effect_fix_est50),sd(effect_fix_est50),sqrt(mean(effect_fix_est50-effect_true50)^2))
fix55_summ = cbind(mean(effect_fix_est55),sd(effect_fix_est55),sqrt(mean(effect_fix_est55-effect_true55)^2))
rbind(fix00_summ,fix05_summ,fix50_summ,fix55_summ)
