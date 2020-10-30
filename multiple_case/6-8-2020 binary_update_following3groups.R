################### efficient way to calculate binary case ##############
SumAllComb <- function(m1, n, x1) {
 # x1 is a vector 
  # n = length(x1) 
  if (m1==0) {
    return(1)
  }else if (m1>n){
    return(0)
  }else if (m1<0) {
    return(0)
  }else if (n==0 & m1>0) {
    return(0)
  }else {
    return(SumAllComb(m1,n-1,x1) + SumAllComb(m1-1,n-1,x1)*x1[n])
  }
}
icpwfunc <- function(exb1, a1, ni){
  # calculate icpw in a culster
  # exb1 is a vector of (exp(X_i1%*%beta1),...,exp(X_ini%*%beta1)) w.r.t. the beta1 in category 1
  # ni is sample size of cluster i
  # a1 is vector of (A_i1,..,A_ini)
  t1 <- sum(a1)
  if(t1==ni | ni==1 | t1==0){
    return(rep(1,ni))
  }else{
    denom <- SumAllComb(m1=t1, n=ni, x1=exb1)
    numer <- sapply(c(1:ni), FUN=function(j){SumAllComb(x1=exb1[-j],m1=(t1-a1[j]),n=(ni-1))})
    numer[a1==1] = numer[a1==1]*exb1[a1==1]
    ps_suffi_cl = numer/denom
    return(ps_suffi_cl)
  }
}

clog_summ = clogit(A~x1+x2+strata(clu), data=data)
beta_est = clog_summ$coefficients
exb1 = exp(as.matrix(data[,4:5])%*%beta_est)
a1 = data$A
data$exb1 = exb1
data$a1 = a1

icpw_all = function(mind,data){
  data_cl1 = data[data$clu==mind,]
  ps_suffi_cl1 = icpwfunc(data_cl1$exb1,data_cl1$a1,n[mind])
  print(mind)
  return(ps_suffi_cl1)
}

ps_suffi2 = unlist(lapply(c(1:m),icpw_all,data))
