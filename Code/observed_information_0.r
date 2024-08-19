########################################################################
#                                                                      #
# This file contains code for calculating the (inverse) observed       #
# information matrix of the new event inference parameters, and        #
# (if relevent) the calibration inference parameters.                  #
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

obs_info_0 = function(p_cal,opt,imle=TRUE,t_cal=NULL)
{
  #
  # FUNCTION INPUTS
  #

  # p_cal: environment storing all objects needed in characterization
  #        calculations
  # opt: object resulting from call to optim()
  # imle: flag for opt from likelihood maximization (TRUE/FALSE)
  # t_cal: object used if bounds supplied to optimization

  #
  # END FUNCTION INPUTS
  #

  Cmat = list()
  Cmat$acov_0 = 0
  if( p_cal$ncalp > 0 ){ Cmat$acov_cal = 0 }

  Hess = -opt$hessian
  if( exists("opt_B",where=p_cal,inherits=FALSE) && p_cal$opt_B ){
    if( is.null(t_cal) ){
      stop("List t_cal must be provided for bounded optimization.")
    } 
    x = p_cal$inv_transform(opt$par[1:p_cal$ntheta0], pc=t_cal)
    if( exists("itransform",where=p_cal,inherits=FALSE) ){
      if( p_cal$itransform ){
        xx = p_cal$inv_tau(x, pc=p_cal)
        Tf = p_cal$j_tau(xx, pc=p_cal)
      }
    }
    dx = rep(1,p_cal$ntheta0)
    itheta0_bounds = t_cal$itheta0_bounds
    if( length(itheta0_bounds[[1]]) > 0 ){
      dx[itheta0_bounds[[1]]] = p_cal$dnotExp(x[itheta0_bounds[[1]]])
    }
    if( length(itheta0_bounds[[2]]) > 0 ){
      dx[itheta0_bounds[[2]]] = -p_cal$dnotExp(x[itheta0_bounds[[2]]])
    }
    if( length(itheta0_bounds[[3]]) > 0 ){
      ttau = p_cal$notExp(x[itheta0_bounds[[3]]])
      dx[itheta0_bounds[[3]]] = p_cal$dnotExp(x[itheta0_bounds[[3]]])*
                                p_cal$theta0_range/(1+ttau)^2
    }
    Tf_b = diag(nrow(Hess))
    diag(Tf_b)[1:p_cal$ntheta0] = dx
    Hess = t(Tf_b) %*% Hess %*% Tf_b
    if( exists("itransform",where=p_cal,inherits=FALSE) ){
      if( p_cal$itransform ){ Hess = t(Tf) %*% Hess %*% Tf }
    }
  }
  Hess = forceSymmetric(Hess)
  cHessCatch = p_cal$tryCatch.W.E(chol(Hess))
  delta = 1.e-17
  idelta = FALSE
  while( !is(cHessCatch$value,"Matrix") ){
    idelta = TRUE
    delta = 10*delta
    dHess = Hess
    diag(dHess) = diag(dHess) + delta
    cHessCatch = p_cal$tryCatch.W.E(chol(dHess))
  }
  cHess = cHessCatch$value
  IHess = chol2inv(cHess)
  if( imle ){
    Cmat$II_nev_it = IHess[1:p_cal$ntheta0, 1:p_cal$ntheta0]
  } else {
    p_cal$IHess = IHess
  }
  if( idelta ){
    print(paste("Perturbation added to Hessian diagonals: ",delta,
                sep=""))
  }

  if( imle ){
    II_nev = Cmat$II_nev_it
    x = opt$par[1:p_cal$ntheta0]
    if( exists("itransform",where=p_cal,inherits=FALSE) ){
      if( p_cal$itransform ){
        Tf = p_cal$j_tau(x, pc=p_cal)
        II_nev = Tf %*% II_nev %*% t(Tf)
        x = p_cal$tau(x, pc=p_cal)
      }
    }
    itheta0_bounds = vector("list",3)
    if( exists("itheta0_bounds",where=p_cal,inherits=FALSE) ){
      itheta0_bounds=p_cal$itheta0_bounds
    } else if( !is.null(t_cal) ){
      itheta0_bounds=t_cal$itheta0_bounds
    }   
    dx = rep(1,p_cal$ntheta0)
    if( length(p_cal$itheta0_bounds[[1]]) > 0 ){
      dx[p_cal$itheta0_bounds[[1]]] =
        p_cal$dnotExp(x[p_cal$itheta0_bounds[[1]]])
    }
    if( length(p_cal$itheta0_bounds[[2]]) > 0 ){
      dx[p_cal$itheta0_bounds[[2]]] =
        -p_cal$dnotExp(x[p_cal$itheta0_bounds[[2]]])
    }
    if( length(p_cal$itheta0_bounds[[3]]) > 0 ){
      tau = p_cal$notExp(x[p_cal$itheta0_bounds[[3]]])
      dx[p_cal$itheta0_bounds[[3]]] =
        p_cal$dnotExp(x[p_cal$itheta0_bounds[[3]]])*
        p_cal$theta0_range/(1+tau)^2
    }
    Tf_b = diag(p_cal$ntheta0)
    diag(Tf_b) = dx
    Cmat$II_nev = Tf_b %*% II_nev %*% t(Tf_b)
    Cmat$acov_0 = 1

    if( p_cal$ncalp > 0 ){
      Cmat$II_calp = IHess[p_cal$ntheta0+1:p_cal$ncalp,
                           p_cal$ntheta0+1:p_cal$ncalp]
      Cmat$acov_cal = 1
    }
  }

  if( imle ){ return(Cmat)
  } else { return(p_cal) }
}
