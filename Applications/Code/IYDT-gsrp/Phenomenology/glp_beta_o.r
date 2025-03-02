########################################################################
#                                                                      #
# This file contains code for computing the gradient of the log-prior  #
# distribution of the optical forward model coefficient parameters     #
# (beta).                                                              #
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

#
# Optical log-prior distribution
#  beta   : vector containing all forward model coefficients
#     p   : list of fixed parameters
#

lq_o <- function(beta,p)
{
  g <- rep(0,length(beta))
  t1 <- p$dnotExp(beta[3])
  g[3] <- -t1/p$notExp(beta[3])+p$d2notExp(beta[3])/t1
  t2 <- p$dnotExp(beta[4])
  g[4] <- -t2/p$notExp(beta[4])+p$d2notExp(beta[4])/t2
  return(g)
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
