########################################################################
#                                                                      #
# This file prints the MLEs of yield and height-of-burst, their        #
# asymptotic covariance matrix, and 95% confidence intervals for new   #
# event yield and height-of-burst.                                     #
#                                                                      #
# © 2023. Triad National Security, LLC. All rights reserved.           #
# This program was produced under U.S. Government contract             #
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which   #
# is operated by Triad National Security, LLC for the U.S. Department  #
# of Energy/National Nuclear Security Administration. All rights in    #
# the program are reserved by Triad National Security, LLC, and the    #
# U.S. Department of Energy/National Nuclear Security Administration.  #
# The Government is granted for itself and others acting on its behalf #
# a nonexclusive, paid-up, irrevocable worldwide license in this       #
# material to reproduce, prepare derivative works, distribute copies   #
# to the public, perform publicly and display publicly, and to permit  #
# others to do so.                                                     #
#                                                                      #
########################################################################

require(Matrix)

# Maximum Likelihood Summaries

mle_0 = p_cal$tmle_0
mle_0[1] = exp(mle_0[1])/10^6
print(mle_0)

err_0 = 2*sqrt(diag(p_cal$Sigma_mle_0$II_nev))
err_0[1] = exp(err_0[1])
print(err_0)
print(signif((err_0[1]-1)*100,2))

if( length(mle_0) > 1 ){
  Var = p_cal$Sigma_mle_0$II_nev
  D = diag(1/sqrt(diag(p_cal$Sigma_mle_0$II_nev)))
  Corr = D %*% Var %*% D
  print(Corr)
}

# Confidence Intervals
z_alpha = qnorm(.975)
theta0_it = p_cal$mle
if( exists("II_nev_it",where=p_cal$Sigma_mle_0,inherits=FALSE) ){
  lb_0 = theta0_it -
         z_alpha*sqrt(diag(p_cal$Sigma_mle_0$II_nev_it))
  ub_0 = theta0_it +
         z_alpha*sqrt(diag(p_cal$Sigma_mle_0$II_nev_it))
} else {
  lb_0 = theta0_it -
         z_alpha*sqrt(diag(p_cal$Sigma_mle_0$II_nev))
  ub_0 = theta0_it +
         z_alpha*sqrt(diag(p_cal$Sigma_mle_0$II_nev))
}
if( exists("itransform",where=p_cal,inherits=FALSE) ){
  if( p_cal$itransform ){
    bmat = rbind(lb_0,ub_0)
    if( p_cal$ntheta0 > 1 ){
      bcall = "expand.grid("
      for( rr in 1:(p_cal$ntheta0-1) ){
        bcall = paste(bcall,"bmat[,",rr,"],",sep="")
      }
      bcall = paste(bcall,"bmat[,p_cal$ntheta0])",sep="")
      bmat = eval(parse(text=bcall))
    }
    tbmat = NULL
    for( rr in 1:nrow(bmat) ){
      tbmat = rbind(tbmat,p_cal$tau(bmat[rr,],pc=p_cal))
    }
    lb_0 = apply(tbmat,2,min)
    ub_0 = apply(tbmat,2,max)
  }
}
if( exists("itheta0_bounds",where=p_cal,inherits=FALSE) ){
  lb_0 = p_cal$transform(lb_0, pc=p_cal)
  ub_0 = p_cal$transform(ub_0, pc=p_cal)
}
lb_0[1] = exp(lb_0[1])/10^6
ub_0[1] = exp(ub_0[1])/10^6
print(lb_0)
print(ub_0)
