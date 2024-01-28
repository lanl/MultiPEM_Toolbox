########################################################################
#                                                                      #
# This file contains code for transforming new event parameters to     #
# half-intervals or intervals and the associated inverse transforms.   #
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

transform <- function(x, pc=p_cal)
{
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    ith0_bds = pc$itheta0_bounds
    if( length(ith0_bds[[1]]) > 0 ){
      x[ith0_bds[[1]]] = pc$theta0_bounds[ith0_bds[[1]],1] +
                         pc$notExp(x[ith0_bds[[1]]])
    }
    if( length(ith0_bds[[2]]) > 0 ){
      x[ith0_bds[[2]]] = pc$theta0_bounds[ith0_bds[[2]],2] -
                         pc$notExp(x[ith0_bds[[2]]])
    }
    if( length(ith0_bds[[3]]) > 0 ){
      tau = pc$notExp(x[ith0_bds[[3]]])
      x[ith0_bds[[3]]] = pc$theta0_bounds[ith0_bds[[3]],1] +
                         pc$theta0_range/(1+1/tau)
    }
  }
  return(x)
}

notExp <- function(x)
{
  f <- x
  ind <- x > 1
  f[ind] <- exp(1) * (x[ind]^2 + 1)/2
  ind <- (x <= 1) & (x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  x[ind] <- -x[ind]
  f[ind] <- exp(1) * (x[ind]^2 + 1)/2
  f[ind] <- 1/f[ind]
  f
}

dnotExp <- function(x)
{
  f <- x
  ind <- x > 1
  f[ind] <- exp(1) * x[ind]
  ind <- (x <= 1) & (x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  x[ind] <- -x[ind]
  f[ind] <- 4 * x[ind] / exp(1) / (x[ind]^2 + 1)^2
  f
}

d2notExp <- function(x)
{
  f <- x
  ind <- x > 1
  f[ind] <- exp(1)
  ind <- (x <= 1) & (x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  x[ind] <- -x[ind]
  f[ind] <- (12 * x[ind]^2 - 4) / exp(1) / (x[ind]^2 + 1)^3
  f
}

inv_transform <- function(x, pc=p_cal)
{
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    ith0_bds = pc$itheta0_bounds
    if( length(ith0_bds[[1]]) > 0 ){
      x[ith0_bds[[1]]] = x[ith0_bds[[1]]] -
                         pc$theta0_bounds[ith0_bds[[1]],1]
      x[ith0_bds[[1]]] = pc$notLog(x[ith0_bds[[1]]])
    }
    if( length(ith0_bds[[2]]) > 0 ){
      x[ith0_bds[[2]]] = pc$theta0_bounds[ith0_bds[[2]],2] -
                         x[ith0_bds[[2]]]
      x[ith0_bds[[2]]] = pc$notLog(x[ith0_bds[[2]]])
    }
    if( length(ith0_bds[[3]]) > 0 ){
      x[ith0_bds[[3]]] = (x[ith0_bds[[3]]] -
                         pc$theta0_bounds[ith0_bds[[3]],1]) /
                         (pc$theta0_bounds[ith0_bds[[3]],2] -
                         x[ith0_bds[[3]]])
      x[ith0_bds[[3]]] = pc$notLog(x[ith0_bds[[3]]])
    }
  }
  return(x)
}

notLog <- function(x) 
{
  f <- x
  ind <- x > exp(1)
  f[ind] <- sqrt(2 * x[ind]/exp(1) - 1)
  ind <- !ind & x > exp(-1)
  f[ind] <- log(x[ind])
  ind <- x <= exp(-1)
  x[ind] <- 1/x[ind]
  f[ind] <- sqrt(2 * x[ind]/exp(1) - 1)
  f[ind] <- -f[ind]
  f
}
