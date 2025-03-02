########################################################################
#                                                                      #
# This file contains code for computing the log-prior distribution of  #
# the theta_0 parameter for the new event.                             #
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
# log yield and HOB/DOB log-prior distribution
#  theta0 : (log yield, HOB/DOB)
#       p : list of fixed parameters
#

lp_0 <- function(theta0,p)
{
  w <- theta0[1]; h <- theta0[2];
  lp0 <- dnorm(w,mean=p$pi_w_mu,sd=p$pi_w_sd,log=TRUE)
  lp0 <- lp0+dnorm(h,mean=p$pi_h_mu,sd=p$pi_h_sd,log=TRUE)
  return(lp0)
}
