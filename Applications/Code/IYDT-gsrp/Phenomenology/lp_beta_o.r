########################################################################
#                                                                      #
# This file contains code for computing the log-prior distribution of  #
# the optical forward model coefficient parameters (beta).             #
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

lp_o <- function(beta,p)
{
  t1 = log(p$notExp(beta[3]))+log(p$notExp(beta[4]))
  return(-t1+log(p$dnotExp(beta[3]))+log(p$dnotExp(beta[4])))
}
