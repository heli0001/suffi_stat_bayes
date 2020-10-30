load("C:/Users/emmal/Desktop/effect_nmig_res00m100_small.RData")
load("C:/Users/emmal/Desktop/effect_nmig_suffi_res00m100_small.RData")
load("C:/Users/emmal/Desktop/effect_suffi00m100_small.RData")
load("C:/Users/emmal/Desktop/ran_res00m100_small.RData")
load("C:/Users/emmal/Desktop/naive_res00m100_small.RData")
load("C:/Users/emmal/Desktop/simdata00m100_small.RData")

m = simdata00m100_small[[1]]$m
datagene=100
effect_suffi = rep(-99,datagene)
effect_nmig = rep(-99,datagene)
effect_nmig_suffi = rep(-99,datagene)
effect_naive = rep(-99,datagene)
effect_ran = rep(-99,datagene)
effect_true = rep(-99,datagene)
for (i in 1:datagene){
  effect_suffi[i] = effect_suffi00m100_small[[i]]$effect_suffi
  #effect_nmig[i] = effect_nmig_res00m100_small[[i]]$effect_nmig
  #effect_nmig_suffi[i] = effect_nmig_suffi_res00m100_small[[i]]$effect_nmig_suffi
  #effect_naive[i] = naive_res00m100_small[[i]]$effec_naive
  #effect_ran[i] = ran_res00m100_small[[i]]$effec_ran
  effect_true[i] = simdata00m100_small[[i]]$simu_true
}
res_se = rbind(sd(effect_naive)/sqrt(datagene),
sd(effect_ran)/sqrt(datagene),
sd(effect_nmig)/sqrt(datagene),
sd(effect_nmig_suffi)/sqrt(datagene),
sd(effect_suffi)/sqrt(datagene))


res_mse = rbind(mean((effect_naive-effect_true)^2),
mean((effect_ran-effect_true)^2),
mean((effect_nmig-effect_true)^2), 
mean((effect_nmig_suffi-effect_true)^2),
mean((effect_suffi-effect_true)^2))


res_bias = rbind(mean(effect_naive)-mean(effect_true),
mean(effect_ran)-mean(effect_true),
mean(effect_nmig)-mean(effect_true),
mean(effect_nmig_suffi)-mean(effect_true),
mean(effect_suffi)-mean(effect_true))



## time span
time_naive = rep(-99,datagene)
time_ran = rep(-99,datagene)
time_suffi = rep(-99,datagene)
time_nmig_suffi = rep(-99,datagene)
time_nmig = rep(-99,datagene)
for (i in 1:datagene){
  #time_naive[i] = naive_res00m100_small[[i]]$time_span[1]
  #time_ran[i] = ran_res00m100_small[[i]]$time_span[1]
  time_suffi[i] = effect_suffi00m100_small[[i]]$time_span[1]
  #time_nmig_suffi[i] = effect_nmig_suffi_res00m100_small[[i]]$time_span[1]
  #time_nmig[i] = effect_nmig_res00m100_small[[i]]$time_span[1]
}

time_clu = rep(-99,datagene)
for (i in 1:datagene){
  time_clu[i] = clu_res00m100_small[[i]]$time_span[1]
}


## 00 elapse: 
res_time = rbind(mean(time_naive),
mean(time_ran),
mean(time_nmig)+mean(time_clu),
mean(time_nmig_suffi)+mean(time_clu),
mean(time_suffi))


######## selected cluster effect ########
sel_clu = rep(-99,datagene)
for (i in 1:datagene){
  sel_clu[i] = length(effect_nmig_suffi_res00m100_small[[i]]$inclu)/m
}
res_select=mean(sel_clu)

ls = list(res_time=res_time,res_bias=res_bias,res_se=res_se,res_mse=res_mse,res_select=res_select)
save(ls,file="ls.RData")
