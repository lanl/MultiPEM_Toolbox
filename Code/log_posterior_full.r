########################################################################
#                                                                      #
# This file contains code for calculating the full log-posterior       #
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

lpost_full = function(x, pc=p_cal)
{
  pc$ll_full(x, pc) + pc$lprior(x, pc)
}
