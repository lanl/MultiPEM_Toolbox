########################################################################
#                                                                      #
# This file contains code for calculating the log-likelihood of the    #
# calibration parameters, with an option to include errors-in-         #
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

ll_cal = function(x, pc)
{
  # use R Matrix package
  require(Matrix)

  # extract inference parameters for new event
  if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
    theta0 = x[1:pc$ntheta0]
    if( exists("itransform",where=pc,inherits=FALSE) ){
      if( pc$itransform ){ theta0 = pc$tau(theta0, pc=pc) }
    }
    if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
      theta0 = pc$transform(theta0, pc=pc)
    }
    x = x[-(1:pc$ntheta0)]
  }
  # extract calibration inference parameters
  if( pc$ncalp > 0 ){
    calp = x[1:pc$ncalp]
    x = x[-(1:pc$ncalp)]
  }
  # extract errors-in-variables yield parameters
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    w_eiv = x[1:pc$nsource]
    x = x[-(1:pc$nsource)]
  }
  # extract common forward model parameters
  if( pc$pbeta > 0 ){
    beta_all = x[1:pc$pbeta]
    x = x[-(1:pc$pbeta)]
  }
  # extract forward model parameters dependent on
  # emplacement condition
  if( pc$ptbeta > 0 ){
    tbeta_all = x[1:pc$ptbeta]
    x = x[-(1:pc$ptbeta)]
  }
  # extract source variance components
  if( pc$pvc_1 > 0 ){
    vc1_all = x[1:pc$pvc_1]
    x = x[-(1:pc$pvc_1)]
  }
  # extract path variance components
  if( pc$pvc_2 > 0 ){
    vc2_all = x[1:pc$pvc_2]
    x = x[-(1:pc$pvc_2)]
  }

  # log-likelihood due to calibration data and (optionally)
  # new event data
  ll = 0

  # compute residuals by emplacement condition
  # compute covariance matrices by source
  for(hh in 1:pc$H){
    # named objects for phenomenology "hh"
    pnames = names(pc$h[[hh]])

    # set up list "pm" with known quantities passed to
    # phenomenology "hh" forward model
    pm = new.env(hash=TRUE)
    if( exists("notExp",where=pc,inherits=FALSE) ){
      pm$notExp = pc$notExp
    }
    if( "llpars" %in% pnames ){
      for( na in names(pc$h[[hh]]$llpars) ){
        eval(parse(text=paste("pm$",na," = ","pc$h[[hh]]$llpars$",
                              na,sep="")))
      }
    }

    # extract forward model parameters for phenomenology "hh"
    if( pc$pbeta > 0 && any(pc$h[[hh]]$pbeta > 0) ){
      pbeta = sum(pc$h[[hh]]$pbeta)
      beta = beta_all[1:pbeta]
      beta_all = beta_all[-(1:pbeta)]
    } else { pbeta = 0 }
    if( pc$ptbeta > 0 && "ptbeta" %in% pnames ){
      ptbeta = sum(pc$h[[hh]]$ptbeta)
      tbeta = tbeta_all[1:ptbeta]
      tbeta_all = tbeta_all[-(1:ptbeta)]
    } else { ptbeta = 0 }

    # number of responses for phenomenology "hh"
    Rh = pc$h[[hh]]$Rh

    # construct source variance component
    # covariance matrices
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      Sigma_hr1 = vector("list",Rh)
      for( rr in 1:Rh ){
        pvc_1 = pc$h[[hh]]$pvc_1[rr]
        if( pvc_1 > 0 ){
          Sigma_hr1[[rr]] = Diagonal(pvc_1,exp(vc1_all[1:pvc_1]))
          vc1_all = vc1_all[-(1:pvc_1)]
        }
      }
    }
    # construct path variance component
    # covariance matrices
    if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
      Sigma_hr2 = vector("list",Rh)
      for( rr in 1:Rh ){
        pvc_2 = pc$h[[hh]]$pvc_2[rr]
        if( pvc_2 > 0 ){
          Sigma_hr2[[rr]] = Diagonal(pvc_2,exp(vc2_all[1:pvc_2]))
          vc2_all = vc2_all[-(1:pvc_2)]
        }
      }
    }
    # construct observational error
    # covariance matrices
    L_h = Diagonal(Rh,exp(x[1:Rh]))
    x = x[-(1:Rh)]
    if( Rh > 1 ){
      L_h[upper.tri(L_h)] = x[1:choose(Rh,2)]
      x = x[-(1:choose(Rh,2))]
    }
    Sigma_h = t(L_h) %*% L_h

    # setup for calibration inference parameters
    if( "cal_par_names" %in% names(pc$h[[hh]]) ){
      csub = which(pc$cal_par_names %in% pc$h[[hh]]$cal_par_names)
      if( length(csub) > 0 ){
        Cp = calp[csub]
        pm$cal_par_names = pc$h[[hh]]$cal_par_names
        pm$ncalp = length(pm$cal_par_names)
      }
    } else { pm$cal = FALSE }

    # iterate over sources in phenomenology "hh"
    for( gg in 1:pc$h[[hh]]$nsource_groups ){
      sc = 0
      tsc = length(pc$h[[hh]]$Source_Groups[[gg]])
      resid = NULL
      Omega = vector("list",tsc)
      for( ii in pc$h[[hh]]$Source_Groups[[gg]] ){
        sc = sc + 1
        # number of responses for source "ii"
        n_hi = pc$h[[hh]]$n[[ii]]
        n_hi_tot = sum(n_hi)

        # setup argument for forward model call
        Arg = "(c(beta_t"
        if( exists("cal_par_names",where=pm,inherits=FALSE) ){
          if( !("nev" %in% pnames) || !pc$h[[hh]]$nev[ii] ){
            pm$cal = TRUE
            Arg = paste(Arg,",Cp",sep="")
          } else { pm$cal = FALSE }
        }
        if( "eiv" %in% pnames ){
          if( !is.null(pc$h[[hh]]$eiv[[ii]]) ){
            W = w_eiv[pc$h[[hh]]$eiv[[ii]]]
            pm$theta_names = "W"
            Arg = paste(Arg,",W",sep="")
          } else {
            if( exists("theta_names",where=pm,inherits=FALSE) ){
              rm(theta_names, envir=pm)
            }
          }
        }
        if( "nev" %in% pnames && pc$h[[hh]]$nev[ii] ){
          if( "itheta0" %in% pnames ){
            Arg = paste(Arg,",theta0[pc$h[[hh]]$itheta0]",sep="")
          } else { Arg = paste(Arg,",theta0",sep="") }
        }
        Arg=paste(Arg,"),pm)",sep="")

        # named parameters in forward model call
        if( "theta_names" %in% pnames &&
            !is.null(pc$h[[hh]]$theta_names[[ii]]) ){
          pm$theta_names = pc$h[[hh]]$theta_names[[ii]]
        }

        # extract forward model parameters for emplacement condition
        # "tt" associated with source "ii"
        if( ptbeta > 0 ){
          for( rr in 1:Rh ){
            if( n_hi[rr] > 0 ){
              tt = as.numeric(as.character(
                   pc$h[[hh]]$X[[ii]][[rr]]$Type[1]))
              break
            } else { next }
          }
          st_betat = 0
          if( tt > 1 ){ st_betat = sum(pc$h[[hh]]$ptbeta[1:(tt-1)]) }
          if( pc$h[[hh]]$ptbeta[tt] > 0 ){
            betat = tbeta[st_betat+(1:pc$h[[hh]]$ptbeta[tt])]
          }
        } else { tt = 1 }

        # covariance matrices
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          Xi_hi1 = vector("list",Rh)
        }

        # iterate over responses "rr"
        for(rr in 1:Rh){
          if( n_hi[rr] > 0 ){
            # covariate matrix for source "ii"
            pm$X = pc$h[[hh]]$X[[ii]][[rr]]

            # extract forward model parameters for response "rr"
            if( pbeta > 0 && pc$h[[hh]]$pbeta[rr] > 0 ){
              st_beta = 0
              if( rr > 1 ){ 
                st_beta = sum(pc$h[[hh]]$pbeta[1:(rr-1)])
              }
              betar = beta[st_beta+(1:pc$h[[hh]]$pbeta[rr])]
            } else { betar = NULL }
            if( ptbeta > 0 && pc$h[[hh]]$pbetat[[tt]][rr] > 0 ){
              st_betatr = 0
              if( rr > 1 ){
                st_betatr = sum(pc$h[[hh]]$pbetat[[tt]][1:(rr-1)])
              }
              betatr = betat[st_betatr+(1:pc$h[[hh]]$pbetat[[tt]][rr])]
            } else { betatr = NULL }
            if( is.null(betatr) ){ beta_t = betar }
            if( is.null(betar) ){ beta_t = betatr }
            if( !is.null(betar) && !is.null(betatr) ){
              beta_t = numeric(pc$h[[hh]]$pbeta[rr]+
                               pc$h[[hh]]$pbetat[[tt]][rr])
              beta_t[pc$h[[hh]]$ibetar[[(tt-1)*Rh+rr]]] = betar
              beta_t[pc$h[[hh]]$ibetatr[[(tt-1)*Rh+rr]]] = betatr
            }
            pm$pbeta = length(beta_t)
            if( "iResponse" %in% pnames ){
              pm$iresp = pc$h[[hh]]$iResponse[rr]
            }

            # calculate forward model
            fcall = paste("pc$ffm$",pc$h[[hh]]$f[rr],Arg,sep="")
            yhat = eval(parse(text=fcall))
            if( any(is.nan(yhat)) ){ return(-Inf) }

            # calculate residual vector
            resid = c(resid,pc$h[[hh]]$Y[[ii]][[rr]] - yhat)

            # calculate components of model covariance matrix
            if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
              if( pc$h[[hh]]$pvc_1[rr] > 0 &&
                  !is.null(pc$h[[hh]]$Z1[[ii]][[rr]]) ){
                Xi_hi1[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]] %*%
                               Sigma_hr1[[rr]] %*%
                               t(pc$h[[hh]]$Z1[[ii]][[rr]])
              } else {
                Xi_hi1[[rr]] = Diagonal(n_hi[rr],0)
              }
            }
          }
        }

        # calculate model covariance matrix
        Omega[[sc]] = Matrix(0,n_hi_tot,n_hi_tot,sparse=FALSE,
                             doDiag=FALSE)
        for( r1 in 1:Rh ){
          if( n_hi[r1] > 0 ){
            st_nir1 = 0
            if( r1 > 1 ){ st_nir1 = sum(n_hi[1:(r1-1)]) }
            ir = st_nir1+(1:n_hi[r1])
            for( r2 in r1:Rh ){
              if( n_hi[r2] > 0 ){
                st_nir2 = 0
                if( r2 > 1 ){ st_nir2 = sum(n_hi[1:(r2-1)]) }
                ic = st_nir2+(1:n_hi[r2])
                Sigma_hi = Matrix(0,n_hi[r1],n_hi[r2],sparse=FALSE,
                                  doDiag=FALSE)
                Sigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[r1]][[r2]]] =
                  Sigma_h[r1,r2]
                Omega[[sc]][ir,ic] = Sigma_hi
                if( r2 > r1 ){ Omega[[sc]][ic,ir] = t(Sigma_hi) }
              }
            }
          }
        }
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          Xi_hi1 = Xi_hi1[!sapply(Xi_hi1,is.null)]
          Omega[[sc]] = Omega[[sc]] + bdiag(Xi_hi1)
        }
      }
      Omega = bdiag(Omega[!sapply(Omega,is.null)])
      n_Omega = nrow(Omega)
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
        for( rr in 1:Rh ){
          if( pc$h[[hh]]$pvc_2[rr] > 0 &&
              !is.null(pc$h[[hh]]$Z2[[gg]][[rr]]) ){
            ic = pc$h[[hh]]$Omega_ic[[gg]][[rr]]
            Omega[ic,ic] = Omega[ic,ic] + pc$h[[hh]]$Z2[[gg]][[rr]] %*%
                           kronecker(Diagonal(pc$h[[hh]]$nplev[gg,rr]),
                                     Sigma_hr2[[rr]]) %*%
                           t(pc$h[[hh]]$Z2[[gg]][[rr]])
          }
        }
      }
      if( any(is.infinite(Omega)) ){ return(-Inf) }
      cOmegaCatch = pc$tryCatch.W.E(chol(Omega))
      if( is(cOmegaCatch$value,"Matrix") ){ cOmega = cOmegaCatch$value
      } else { return(-Inf) }
      if( kappa(cOmega) > 1000 ){ return(-Inf) }

      # calculate components of log-likelihood function
      ll = ll - sum(log(diag(cOmega)))
      ll = ll - t(resid) %*% chol2inv(cOmega) %*% resid/2
      ll = ll - n_Omega*log(2*pi)/2
    }
  }
  return(as.numeric(ll))
}
