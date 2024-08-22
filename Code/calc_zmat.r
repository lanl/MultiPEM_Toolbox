########################################################################
#                                                                      #
# This file contains code for computing default level 1 and level 2    #
# variance component coefficient matrices.                             #
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

calc_zmat = function(pc=p_cal)
{
  # Variance component coefficient matrices
  if( pc$pvc_1 > 0 ){
    for( hh in 1:pc$H ){
      if( any(pc$h[[hh]]$pvc_1 > 0) ){
        # Level 1
        pc$h[[hh]]$Z1 = vector("list",pc$h[[hh]]$nsource)
        for( qq in 1:pc$h[[hh]]$nsource ){
          pc$h[[hh]]$Z1[[qq]] = vector("list",pc$h[[hh]]$Rh)
          for( rr in 1:pc$h[[hh]]$Rh ){
            if( pc$h[[hh]]$pvc_1[rr] > 0 &&
                pc$h[[hh]]$n[[qq]][rr] > 0 ){
              pc$h[[hh]]$Z1[[qq]][[rr]] =
              Matrix(rep(1,pc$h[[hh]]$n[[qq]][rr]),ncol=1)
            }
          }
        }
      }
      if( pc$pvc_2 > 0 ){
        if( any(pc$h[[hh]]$pvc_1 > 0 & pc$h[[hh]]$pvc_2 > 0) ){
          # Level 2
          pc$h[[hh]]$Z2 = vector("list",pc$h[[hh]]$nsource)
          for( qq in 1:pc$h[[hh]]$nsource ){
            pc$h[[hh]]$Z2[[qq]] = vector("list",pc$h[[hh]]$Rh)
            for( rr in 1:pc$h[[hh]]$Rh ){
              if( pc$h[[hh]]$pvc_1[rr] > 0 ){
                if( pc$h[[hh]]$pvc_2[rr] > 0 &&
                    pc$h[[hh]]$n[[qq]][rr] > 0 ){
                  pc$h[[hh]]$Z2[[qq]][[rr]] =
                  Matrix(0,pc$h[[hh]]$n[[qq]][rr],
                           pc$h[[hh]]$nplev[qq,rr],doDiag=FALSE)
                  for( ss in 1:pc$h[[hh]]$nplev[qq,rr] ){
                    st_nhijr = 0
                    if( ss > 1 ){
                      st_nhijr =
                        sum(pc$nh[[hh]]$i[[qq]]$r[[rr]][1:(ss-1)])
                    }
                    pc$h[[hh]]$Z2[[qq]][[rr]][st_nhijr+
                      (1:pc$nh[[hh]]$i[[qq]]$r[[rr]][ss]),ss] =
                    Matrix(rep(1,
                               pc$nh[[hh]]$i[[qq]]$r[[rr]][ss],ncol=1))
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  return(pc)
}
