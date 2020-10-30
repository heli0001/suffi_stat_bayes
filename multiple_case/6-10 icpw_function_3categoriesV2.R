SumAllComb <- function(m1, m2, n, x1, x2) {
  # this function is to get sum of combinations of (m1,m2) elements from x1 and x2
  # x1 and x2 are vector 
  # n = length(x1) = length(x2)
  mat <- array(0, c(m1+2,m2+2,n+1))
    mat[2,2,] = 1
    for (l in seq(n)) {
      for (i in seq(0,length=min(m1,l)+1)) {
        for (j in seq(0,length=min(l-i,m2)+1)) {
          if ((i>0) | (j>0)) {
            mat[i+2,j+2,l+1] <- mat[i+2,j+2,l] + mat[i+1,j+2,l]*x1[l] + mat[i+2,j+1,l]*x2[l]
          }
        }
      }
    }
    return(mat[m1+2,m2+2,n+1])
  }

icpwfunc <- function(exb1, exb2, a1, a2, ni){
  # calculate icpw in a culster
  # exb1 is a vector of (exp(X_i1%*%beta1),...,exp(X_ini%*%beta1)) w.r.t. the beta1 in category 1
  # exb2 is a vector of (exp(X_i1%*%beta2),...,exp(X_ini%*%beta2)) w.r.t. the beta2 in category 2
  # ni is sample size of cluster i
  # a is vector of (A_i1,..,A_ini)
  t1 <- sum(a1)
  t2 <- sum(a2)
  if(t1==ni | t2==ni | ni==1 | t1==0 | t2==0){
    return(rep(1,ni))
  }else{
    denom <- SumAllComb(m1=t1, m2=t2, n=ni, x1=exb1, x2=exb2)
    numer <- sapply(c(1:ni), FUN=function(j){SumAllComb(x1=exb1[-j],x2=exb2[-j],m1=(t1-a1[j]),m2=(t2-a2[j]),n=(ni-1))})
    numer1[a1==1] <- numer[a1==1]*exb[a1==1]
    numer2[a2==1] <- numer[a2==1]*exb[a2==1]
    ps1 <- denom/numer1
    ps2 <- denom/numer2
    ps_list <- list(ps1,ps2)
    return(ps_list)
  }
}