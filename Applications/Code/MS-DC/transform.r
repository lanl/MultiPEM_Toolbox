########################################################################
#                                                                      #
# This file contains code for a user-provided transform of new event   #
# parameters and the associated jacobian and inverse transform.        #
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

# transform
tau <- function(x, pc=p_cal)
{
  x[2] = x[2]*exp(pc$tpars$yield_scaling*x[1])
  return(x)
}

# Jacobian
j_tau <- function(x, pc=p_cal)
{
  Jac = matrix(0,pc$ntheta0,pc$ntheta0)
  Jac[1,1] = 1
  w_sc = exp(pc$tpars$yield_scaling*x[1])
  Jac[2,1] = pc$tpars$yield_scaling*x[2]*w_sc
  Jac[2,2] = w_sc
  return(Jac)
}

# log abs(det(Jacobian))
log_absdet_j_tau <- function(x, pc=p_cal)
{
  return(pc$tpars$yield_scaling*x[1])
}

# gradient log abs(det(Jacobian))
dlog_absdet_j_tau <- function(x, pc=p_cal)
{
  return(c(pc$tpars$yield_scaling,0))
}

# inverse transform
inv_tau <- function(x, pc=p_cal)
{
  x[2] = x[2]/exp(pc$tpars$yield_scaling*x[1])
  return(x)
}
