########################################################################
#                                                                      #
# This file contains code for calculating the log-likelihood of the    #
# errors-in-variables yield parameters.                                #
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

ll_eiv = function(x, pc)
{
  # log-likelihood due to errors-in-variables constraints
  resid = pc$eiv_w - x
  ll = -0.5*sum((resid/pc$eiv_w_sd)^2)
  ll = ll - pc$nsource*log(2*pi)/2

  return(ll)
}
