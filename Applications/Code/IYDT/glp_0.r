########################################################################
#                                                                      #
# This file contains code for computing the gradient of the log-prior  #
# distribution of the theta_0 parameter for the new event.             #
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

lq_0 <- function(theta0,p)
{
  g <- rep(0,length(theta0))
  w <- theta0[1]; h <- theta0[2];
  g[1] <- (p$pi_w_mu - w)/p$pi_w_sd^2
  g[2] <- (p$pi_h_mu - h)/p$pi_h_sd^2
  return(g)
}
