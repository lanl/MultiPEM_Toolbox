########################################################################
#                                                                      #
# This file contains code for computing default level 1 and level 2    #
# variance component coefficient matrices for new event data.          #
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

calc_zmat_0 = function(pc=p_cal)
{
  # Variance component coefficient matrices
  if( pc$pvc_1 > 0 ){
    for( hh in 1:pc$H ){
      nsource = pc$h[[hh]]$nsource
      if( any(pc$h[[hh]]$pvc_1 > 0) ){
        # Source
        pc$h[[hh]]$Z1 = c(pc$h[[hh]]$Z1,vector("list",1))
        pc$h[[hh]]$Z1[[nsource]] = vector("list",pc$h[[hh]]$Rh)
        for( rr in 1:pc$h[[hh]]$Rh ){
          if( pc$h[[hh]]$pvc_1[rr] > 0 &&
              pc$h[[hh]]$n[[nsource]][rr] > 0 ){
            pc$h[[hh]]$Z1[[nsource]][[rr]] =
            Matrix(rep(1,pc$h[[hh]]$n[[nsource]][rr]),ncol=1)
          }
        }
      }
    }
  }
  if( pc$pvc_2 > 0 ){
    for( hh in 1:pc$H ){
      nsource = pc$h[[hh]]$nsource_groups
      if( any(pc$h[[hh]]$pvc_2 > 0) ){
        # Path
        pc$h[[hh]]$Z2 = c(pc$h[[hh]]$Z2,vector("list",1))
        pc$h[[hh]]$Z2[[nsource]] = vector("list",pc$h[[hh]]$Rh)
        for( rr in 1:pc$h[[hh]]$Rh ){
          if( pc$h[[hh]]$pvc_2[rr] > 0 &&
              pc$h[[hh]]$ng[[nsource]][rr] > 0 ){
            pc$h[[hh]]$Z2[[nsource]][[rr]] =
            Matrix(0,pc$h[[hh]]$ng[[nsource]][rr],
                     pc$h[[hh]]$nplev[nsource,rr],doDiag=FALSE)
            for( ss in 1:pc$h[[hh]]$nplev[nsource,rr] ){
              ir = p_cal$nh[[hh]]$i[[nsource]]$r[[rr]]$p[[ss]]
              pc$h[[hh]]$Z2[[nsource]][[rr]][ir,ss] = 1
            }
          }
        }
      }
    }
  }
  return(pc)
}
