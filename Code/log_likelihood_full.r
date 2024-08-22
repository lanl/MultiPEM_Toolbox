########################################################################
#                                                                      #
# This file contains code for calculating the full log-likelihood      #
# of the calibration parameters, including (optionally) errors-in-     #
# variables for calibration data yields. New event data can also be    #
# optionally incorporated, allowing for joint inference of new event   #
# and calibration parameters based on calibration and new event data.  #
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

ll_full = function(x, pc=p_cal)
{
  ll = pc$ll_cal(x, pc)

  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
      x = x[-(1:pc$ntheta0)]
    }
    if( pc$ncalp > 0 ){
      x = x[-(1:pc$ncalp)]
    }
    x = x[1:pc$nsource]
    ll = ll + pc$ll_eiv(x, pc)
  }

  return(ll)
}
