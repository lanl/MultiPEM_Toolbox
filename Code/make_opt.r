########################################################################
#                                                                      #
# This file contains code for organizing the parameter values provided #
# to the function.                                                     #
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

make_opt = function(xfin, pc)
{
  # list of parameter values
  opt = list()

  # collect new event inference parameters
  if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
    opt$theta0 = xfin[1:pc$ntheta0]
    xfin = xfin[-(1:pc$ntheta0)]
  }

  # collect calibration inference parameters
  if( pc$ncalp > 0 ){
    xfin_calp = xfin[1:pc$ncalp]
    names(xfin_calp) = pc$cal_par_names
    opt$calp = xfin_calp
    xfin = xfin[-(1:pc$ncalp)]
  }

  # collect errors-in-variables inference parameters
  if( pc$nsource > 0 ){
    xfin_eiv = xfin[1:pc$nsource]
    names(xfin_eiv) = pc$eiv_sources
    opt$w_eiv = xfin_eiv
    xfin = xfin[-(1:pc$nsource)]
  }

  # collect common coefficients by response
  if( pc$pbeta > 0 ){
    opt$beta = NULL
    for( hh in 1:pc$H ){
      pbeta = sum(pc$h[[hh]]$pbeta)
      if( pbeta > 0 ){
        opt$beta = c(opt$beta, xfin[1:pbeta])
        xfin = xfin[-(1:pbeta)]
      }
    }
  }

  # collect emplacement condition dependent coefficients by response
  if( pc$ptbeta > 0 ){
    opt$tbeta = NULL
    for( hh in 1:pc$H ){
      if( "ptbeta" %in% names(pc$h[[hh]]) ){
        ptbeta = sum(pc$h[[hh]]$ptbeta)
      } else { ptbeta = 0 }
      if( ptbeta > 0 ){
        opt$tbeta = c(opt$tbeta, xfin[1:ptbeta])
        xfin = xfin[-(1:ptbeta)]
      }
    }
  }

  # collect variance components by response
  if( pc$pvc_1 > 0 ){
    opt$vc_1 = NULL
    for( hh in 1:pc$H ){
      pvc_1 = sum(pc$h[[hh]]$pvc_1)
      if( pvc_1 > 0 ){
        opt$vc_1 = c(opt$vc_1, xfin[1:pvc_1])
        xfin = xfin[-(1:pvc_1)]
      }
    }
  }
  if( pc$pvc_2 > 0 ){
    opt$vc_2 = NULL
    for( hh in 1:pc$H ){
      pvc_2 = sum(pc$h[[hh]]$pvc_2)
      if( pvc_2 > 0 ){
        opt$vc_2 = c(opt$vc_2, xfin[1:pvc_2])
        xfin = xfin[-(1:pvc_2)]
      }
    }
  }

  # collect observational error covariance parameters
  opt$eps = NULL
  for( hh in 1:pc$H ){
    ncpar = pc$h[[hh]]$Rh*(pc$h[[hh]]$Rh+1)/2
    opt$eps = c(opt$eps, xfin[1:ncpar])
    xfin = xfin[-(1:ncpar)]
  }

  return(opt)
}
