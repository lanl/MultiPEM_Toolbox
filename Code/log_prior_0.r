########################################################################
#                                                                      #
# This file contains code for calculating the log-prior density of the #
# new event parameters.                                                #
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

lprior_0 = function(x, pc)
{
  # log-prior density
  lp = 0

  # extract inference parameters for new event
  if( exists("itransform",where=pc,inherits=FALSE) ){
    if( pc$itransform ){
      lp = lp + pc$log_absdet_j_tau(x, pc=pc)
      x = pc$tau(x, pc=pc)
    }
  }
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    ith0_bds = pc$itheta0_bounds
    if( length(ith0_bds[[1]]) > 0 ){
      jt = sum(log(pc$dnotExp(x[ith0_bds[[1]]])))
      lp = lp + jt
    }
    if( length(ith0_bds[[2]]) > 0 ){
      jt = sum(log(pc$dnotExp(x[ith0_bds[[2]]])))
      lp = lp + jt
    }
    if( length(ith0_bds[[3]]) > 0 ){
      tau = pc$notExp(x[ith0_bds[[3]]])
      jt = pc$sum_theta0_logrange +
           sum(log(pc$dnotExp(x[ith0_bds[[3]]]))) -
           2*sum(log(1+tau))
      lp = lp + jt
    }
    x = pc$transform(x, pc=pc)
  }
  if( "lp_theta0" %in% names(pc) ){
    Arg = "(x,pc)"
    lp_theta0_call = paste("pc$flp$",pc$lp_theta0$f,Arg,sep="")
    # evaluate log-prior for new event inference parameters
    lp = lp + eval(parse(text=lp_theta0_call))
  }
  return(as.numeric(lp))
}
