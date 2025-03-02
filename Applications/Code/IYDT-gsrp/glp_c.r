########################################################################
#                                                                      #
# This file contains code for computing the gradient of the log-prior  #
# distribution of the calibration inference parameters.                #
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
# log chemical-to-nuclear equivalency log-prior distribution
#  calp : (log C2N)
#     p : list of fixed parameters
#

lq_c <- function(calp,p)
{
  g <- rep(0,length(calp))
  c <- calp[1]
  g[1] <- (p$pi_c_mu - c)/p$pi_c_sd^2
  return(g)
}
