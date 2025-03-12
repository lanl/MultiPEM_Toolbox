########################################################################
#                                                                      #
# This file contains code for collecting fits from forward models and  #
# fitted components of mixed effect covariance matrices.               #
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

predict = function(x, pc)
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

  # create prediction environment
  pred_out = new.env(hash=TRUE)
  pred_out$h = vector("list",pc$H)

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

    # collect data and fits
    pred_out$h[[hh]]$observed = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$fitted = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$Sigma_eps = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$Omega = vector("list",pc$h[[hh]]$nsource_groups)
    pred_out$h[[hh]]$IOmega = vector("list",pc$h[[hh]]$nsource_groups)
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      pred_out$h[[hh]]$Xi_1 = vector("list",pc$h[[hh]]$nsource)
      pred_out$h[[hh]]$Sigma_vc1 = vector("list",pc$h[[hh]]$nsource)
      pred_out$h[[hh]]$mu_b1 = vector("list",Rh)
      pred_out$h[[hh]]$sd_b1 = vector("list",Rh)
      snames = vector("list",Rh)
      pred_out$h[[hh]]$source_effects =
        vector("list",pc$h[[hh]]$nsource)
    }
    if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
      pred_out$h[[hh]]$Xi_2 = vector("list",pc$h[[hh]]$nsource_groups)
      pred_out$h[[hh]]$Sigma_vc2 =
        vector("list",pc$h[[hh]]$nsource_groups)
      pred_out$h[[hh]]$mu_b2 = vector("list",Rh)
      pred_out$h[[hh]]$sd_b2 = vector("list",Rh)
      panames = vector("list",Rh)
      pred_out$h[[hh]]$path_effects =
        vector("list",pc$h[[hh]]$nsource)
    }

    # iterate over sources in phenomenology "hh"
    for( gg in 1:pc$h[[hh]]$nsource_groups ){
      sc = 0
      tsc = length(pc$h[[hh]]$Source_Groups[[gg]])
      resid = NULL
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        Sigma_h1 = vector("list",tsc)
        Z1 = vector("list",tsc)
      }
      Omega = vector("list",tsc)
      Lcov = vector("list",tsc)
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

        # extract forward model parameters for emplacement condition "tt"
        # associated with source "ii"
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
          Z1_s = vector("list",Rh)
          Xi_hi1 = vector("list",Rh)
        }

        # collect observations and associated fits
        pred_out$h[[hh]]$observed[[ii]] = vector("list",Rh)
        pred_out$h[[hh]]$fitted[[ii]] = vector("list",Rh)
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          pred_out$h[[hh]]$Xi_1[[ii]] = vector("list",Rh)
          pred_out$h[[hh]]$source_effects[[ii]] = vector("list",Rh)
        }
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
          pred_out$h[[hh]]$path_effects[[ii]] = vector("list",Rh)
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

            # calculate residual vector
            resid = c(resid,pc$h[[hh]]$Y[[ii]][[rr]] - yhat)

            # collect observed data and fitted values
            pred_out$h[[hh]]$observed[[ii]][[rr]] =
              pc$h[[hh]]$Y[[ii]][[rr]]
            pred_out$h[[hh]]$fitted[[ii]][[rr]] = yhat

            # calculate components of model covariance matrix
            if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
              if( pc$h[[hh]]$pvc_1[rr] > 0 &&
                  !is.null(pc$h[[hh]]$Z1[[ii]][[rr]]) ){
                Z1_s[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]]
                Xi_hi1[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]] %*%
                               Sigma_hr1[[rr]] %*%
                               t(pc$h[[hh]]$Z1[[ii]][[rr]])
              } else {
                Xi_hi1[[rr]] = Diagonal(n_hi[rr],0)
              }
              pred_out$h[[hh]]$Xi_1[[ii]][[rr]] = Xi_hi1[[rr]]
            }
          }
        }

        # components of posterior source effects covariance matrix
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          Sigma_h1[[sc]] = bdiag(Sigma_hr1[!sapply(Sigma_hr1,is.null)])
          Z1[[sc]] = bdiag(Z1_s[!sapply(Z1_s,is.null)])
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
        pred_out$h[[hh]]$Sigma_eps[[ii]] = Omega[[sc]]
        Lcov[[sc]] = Omega[[sc]]
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          Xi_hi1 = Xi_hi1[!sapply(Xi_hi1,is.null)]
          pred_out$h[[hh]]$Sigma_vc1[[ii]] = bdiag(Xi_hi1)
          Omega[[sc]] = Omega[[sc]] + bdiag(Xi_hi1)
        }
      }
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        Sigma_h1 = bdiag(Sigma_h1[!sapply(Sigma_h1,is.null)])
        Z1 = bdiag(Z1[!sapply(Z1,is.null)])
      }
      Omega = bdiag(Omega[!sapply(Omega,is.null)])
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
        Sigma_h2 = vector("list",Rh)
        pred_out$h[[hh]]$Xi_2[[gg]] = vector("list",Rh)
        n_Omega = nrow(Omega)
        Xi_hg2 = Matrix(0,n_Omega,n_Omega,sparse=FALSE,doDiag=FALSE)
        Z2 = NULL
        for( rr in 1:Rh ){
          if( pc$h[[hh]]$pvc_2[rr] > 0 &&
              !is.null(pc$h[[hh]]$Z2[[gg]][[rr]]) ){
            Sigma_h2[[rr]] =
              kronecker(Diagonal(pc$h[[hh]]$nplev[gg,rr]),
                        Sigma_hr2[[rr]])
            pred_out$h[[hh]]$Xi_2[[gg]][[rr]] =
                             pc$h[[hh]]$Z2[[gg]][[rr]] %*%
                             Sigma_h2[[rr]] %*%
                             t(pc$h[[hh]]$Z2[[gg]][[rr]])
            ic = pc$h[[hh]]$Omega_ic[[gg]][[rr]]
            Xi_hg2[ic,ic] = pred_out$h[[hh]]$Xi_2[[gg]][[rr]]
            tZ2 = Matrix(0,n_Omega,nrow(Sigma_h2[[rr]]),sparse=FALSE,
                         doDiag=FALSE)
            tZ2[ic,] = pc$h[[hh]]$Z2[[gg]][[rr]]
            Z2 = cbind(Z2,tZ2) 
          }
        }
        pred_out$h[[hh]]$Sigma_vc2[[gg]] = Xi_hg2
        Omega = Omega + Xi_hg2
        Sigma_h2 = bdiag(Sigma_h2[!sapply(Sigma_h2,is.null)])
        B_2 = Sigma_h2 %*% t(Z2)
      }
      pred_out$h[[hh]]$Omega[[gg]] = Omega
      cOmega = chol(Omega)
      pred_out$h[[hh]]$IOmega[[gg]] = chol2inv(cOmega)

      # calculate posterior distributions of source and path effects
      Lcov = bdiag(Lcov[!sapply(Lcov,is.null)])
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        Lcov_1 = Lcov
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
          Lcov_1 = Lcov_1 + pred_out$h[[hh]]$Sigma_vc2[[gg]]
        }
        Catch = pc$tryCatch.W.E(chol(Lcov_1))
        if( is(Catch$value,"Matrix") ){ cLcov_1 = Catch$value
        } else {
          cat("Lcov_1 is not positive definite.\n")
          Lcov_1 = nearPD(Lcov_1)$mat
          cLcov_1 = chol(Lcov_1)
        }
        ILcov_1 = chol2inv(cLcov_1)
        B_1 = Sigma_h1 %*% t(Z1)
        ICov_1 = Lcov_1 + Z1 %*% B_1
        Catch = pc$tryCatch.W.E(chol(ICov_1))
        if( is(Catch$value,"Matrix") ){ cICov_1 = Catch$value
        } else {
          cat("ICov_1 is not positive definite.\n")
          ICov_1 = nearPD(ICov_1)$mat
          cICov_1 = chol(ICov_1)
        }
        Cov_1 = chol2inv(cICov_1)
        Cov_1 = Sigma_h1 - B_1 %*% Cov_1 %*% t(B_1)
        Catch = pc$tryCatch.W.E(chol(Cov_1))
        if( !is(Catch$value,"Matrix") ){
          cat(paste("Source effect posterior covariance matrix ",
              "is not positive definite.\n",sep=""))
          Cov_1 = nearPD(Cov_1)$mat
        }
        mu_b1 = t(Z1) %*% ILcov_1 %*% resid
        mu_b1 = Cov_1 %*% mu_b1
        seffects = Z1 %*% mu_b1
        sd_b1 = sqrt(diag(Cov_1))
        kk = 0
        for( ii in pc$h[[hh]]$Source_Groups[[gg]] ){
          for( rr in 1:Rh ){
            if( pc$h[[hh]]$n[[ii]][rr] > 0 ){
              ib1 = kk+1:nrow(Sigma_hr1[[rr]])
              pred_out$h[[hh]]$mu_b1[[rr]] =
                c(pred_out$h[[hh]]$mu_b1[[rr]],mu_b1[ib1])
              pred_out$h[[hh]]$sd_b1[[rr]] =
                c(pred_out$h[[hh]]$sd_b1[[rr]],sd_b1[ib1])
              kk = kk + nrow(Sigma_hr1[[rr]])
            }
          }
        }
        for( rr in 1:Rh ){
          snames[[rr]] = c(snames[[rr]],
                           rep(pc$h[[hh]]$Source[[gg]][[rr]],
                               each=nrow(Sigma_hr1[[rr]])))
          if( "Omega_ic" %in% pnames ){
            ic = pc$h[[hh]]$Omega_ic[[gg]][[rr]]
          } else {
            st_ng = 0
            if( rr > 1 ){
              st_ng = st_ng + sum(pc$h[[hh]]$ng[[gg]][1:(rr-1)])
            }
            ic = st_ng + (1:pc$h[[hh]]$ng[[gg]][rr])
          }
          tseffects = seffects[ic]
          kk = 0
          for( ii in pc$h[[hh]]$Source_Groups[[gg]] ){
            if( pc$h[[hh]]$n[[ii]][rr] > 0 ){
              iseffects = kk + 1:pc$h[[hh]]$n[[ii]][rr]
              pred_out$h[[hh]]$source_effects[[ii]][[rr]] =
                tseffects[iseffects]
              kk = kk + pc$h[[hh]]$n[[ii]][rr]
            }
          }
        }
      }
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
        Lcov_2 = Lcov
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          Lcov_2 = Lcov_2 + t(B_1) %*% t(Z1)
        }
        Catch = pc$tryCatch.W.E(chol(Lcov_2))
        if( is(Catch$value,"Matrix") ){ cLcov_2 = Catch$value
        } else {
          cat("Lcov_2 is not positive definite.\n")
          Lcov_2 = nearPD(Lcov_2)$mat
          cLcov_2 = chol(Lcov_2)
        }
        ILcov_2 = chol2inv(cLcov_2)
        ICov_2 = Lcov_2 + pred_out$h[[hh]]$Sigma_vc2[[gg]]
        Catch = pc$tryCatch.W.E(chol(ICov_2))
        if( is(Catch$value,"Matrix") ){ cICov_2 = Catch$value
        } else {
          cat("ICov_2 is not positive definite.\n")
          ICov_2 = nearPD(ICov_2)$mat
          cICov_2 = chol(ICov_2)
        }
        Cov_2 = chol2inv(cICov_2)
        Cov_2 = Sigma_h2 - B_2 %*% Cov_2 %*% t(B_2)
        Catch = pc$tryCatch.W.E(chol(Cov_2))
        if( !is(Catch$value,"Matrix") ){
          cat(paste("Path effect posterior covariance matrix ",
              "is not positive definite.\n",sep=""))
          Cov_2 = nearPD(Cov_2)$mat
        }
        mu_b2 = t(Z2) %*% ILcov_2 %*% resid
        mu_b2 = Cov_2 %*% mu_b2
        peffects = Z2 %*% mu_b2
        sd_b2 = sqrt(diag(Cov_2))
        kk = 0
        for( rr in 1:Rh ){
          if( pc$h[[hh]]$ng[[gg]][rr] > 0 ){
            ib2 = kk+1:(pc$h[[hh]]$nplev[gg,rr]*nrow(Sigma_hr2[[rr]]))
            pred_out$h[[hh]]$mu_b2[[rr]] =
              c(pred_out$h[[hh]]$mu_b2[[rr]],mu_b2[ib2])
            pred_out$h[[hh]]$sd_b2[[rr]] =
              c(pred_out$h[[hh]]$sd_b2[[rr]],sd_b2[ib2])
            kk = kk + pc$h[[hh]]$nplev[gg,rr]*nrow(Sigma_hr2[[rr]])
          }
        }
        for( rr in 1:Rh ){
          panames[[rr]] = c(panames[[rr]],
                            rep(pc$h[[hh]]$Path[[gg]][[rr]],
                                each=nrow(Sigma_hr2[[rr]])))
          ic = pc$h[[hh]]$Omega_ic[[gg]][[rr]]
          tpeffects = peffects[ic]
          kk = 0
          for( ii in pc$h[[hh]]$Source_Groups[[gg]] ){
            if( pc$h[[hh]]$n[[ii]][rr] > 0 ){
              ipeffects = kk + 1:pc$h[[hh]]$n[[ii]][rr]
              pred_out$h[[hh]]$path_effects[[ii]][[rr]] =
                tpeffects[ipeffects]
              kk = kk + pc$h[[hh]]$n[[ii]][rr]
            }
          }
        }
      }
    }
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      for( rr in 1:Rh ){
        names(pred_out$h[[hh]]$mu_b1[[rr]]) = snames[[rr]]
        names(pred_out$h[[hh]]$sd_b1[[rr]]) = snames[[rr]]
      }
    }
    if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
      for( rr in 1:Rh ){
        names(pred_out$h[[hh]]$mu_b2[[rr]]) = panames[[rr]]
        names(pred_out$h[[hh]]$sd_b2[[rr]]) = panames[[rr]]
      }
    }
  }
  return(pred_out)
}
