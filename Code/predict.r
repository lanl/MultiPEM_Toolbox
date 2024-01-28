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
  # extract level 1 variance components
  if( pc$pvc_1 > 0 ){
    vc1_all = x[1:pc$pvc_1]
    x = x[-(1:pc$pvc_1)]
  }
  # extract level 2 variance components
  if( pc$pvc_1 > 0 && pc$pvc_2 > 0 ){
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

    # construct level 1 variance component
    # covariance matrices
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      Sigma_hr1 = vector("list",Rh)
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                              pc$h[[hh]]$pvc_2 > 0) ){
        Sigma_hr2 = vector("list",Rh)
      }
      for( rr in 1:Rh ){
        pvc_1 = pc$h[[hh]]$pvc_1[rr]
        if( pvc_1 > 0 ){
          Sigma_hr1[[rr]] = Diagonal(pvc_1,exp(vc1_all[1:pvc_1]))
          vc1_all = vc1_all[-(1:pvc_1)]
          # construct level 2 variance component
          # covariance matrices
          if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                  pc$h[[hh]]$pvc_2 > 0) ){
            pvc_2 = pc$h[[hh]]$pvc_2[rr]
            if( pvc_2 > 0 ){
              Sigma_hr2[[rr]] = Diagonal(pvc_2,exp(vc2_all[1:pvc_2]))
              vc2_all = vc2_all[-(1:pvc_2)]
            }
          }
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

    # collect data and fits
    pred_out$h[[hh]]$observed = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$fitted = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$Sigma_eps = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$Omega = vector("list",pc$h[[hh]]$nsource)
    pred_out$h[[hh]]$IOmega = vector("list",pc$h[[hh]]$nsource)
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      pred_out$h[[hh]]$Xi_1 = vector("list",pc$h[[hh]]$nsource)
      pred_out$h[[hh]]$Sigma_vc1 = vector("list",pc$h[[hh]]$nsource)
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                              pc$h[[hh]]$pvc_2 > 0) ){
        pred_out$h[[hh]]$Xi_2 = vector("list",pc$h[[hh]]$nsource)
        pred_out$h[[hh]]$Sigma_vc2 = vector("list",pc$h[[hh]]$nsource)
      }
    }

    # iterate over sources in phenomenology "hh"
    for(ii in 1:pc$h[[hh]]$nsource){
      # number of responses for source "ii"
      n_hi = pc$h[[hh]]$n[[ii]]
      n_hi_tot = sum(n_hi)

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
        Xi_hi1 = vector("list",Rh)
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          Xi_hi2 = vector("list",Rh)
        }
      }

      # collect observations and associated fits
      pred_out$h[[hh]]$observed[[ii]] = vector("list",Rh)
      pred_out$h[[hh]]$fitted[[ii]] = vector("list",Rh)
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        pred_out$h[[hh]]$Xi_1[[ii]] = vector("list",Rh)
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          pred_out$h[[hh]]$Xi_2[[ii]] = vector("list",Rh)
        }
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

          # prepare call to forward model
          Arg = "(beta_t,pm)"
          if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
            W = w_eiv[pc$h[[hh]]$eiv[[ii]]]
            if( !(exists("theta_names",where=pm,inherits=FALSE)) ){
              pm$theta_names = "W"
            }
            Arg = "(c(beta_t,W),pm)"
          }
          if( "nev" %in% pnames && pc$h[[hh]]$nev[ii] ){
            if( "itheta0" %in% pnames ){
              Arg = "(c(beta_t,theta0[pc$h[[hh]]$itheta0]),pm)"
            } else { Arg = "(c(beta_t,theta0),pm)" }
          }
          fcall = paste("pc$ffm$",pc$h[[hh]]$f[rr],Arg,sep="")

          # collect observed data and fitted values
          pred_out$h[[hh]]$observed[[ii]][[rr]] =
            pc$h[[hh]]$Y[[ii]][[rr]]
          pred_out$h[[hh]]$fitted[[ii]][[rr]] = eval(parse(text=fcall))

          # calculate components of model covariance matrix
          if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
            if( pc$h[[hh]]$pvc_1[rr] > 0 ){
              Xi_hi1[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]] %*%
                             Sigma_hr1[[rr]] %*%
                             t(pc$h[[hh]]$Z1[[ii]][[rr]])
            } else {
              Xi_hi1[[rr]] = Diagonal(n_hi[rr],0)
            }
            pred_out$h[[hh]]$Xi_1[[ii]][[rr]] = Xi_hi1[[rr]]
            if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                    pc$h[[hh]]$pvc_2 > 0) ){
              if( pc$h[[hh]]$pvc_1[rr] > 0 ){
                if( pc$h[[hh]]$pvc_2[rr] > 0 ){
                  Xi_hi2[[rr]] = pc$h[[hh]]$Z2[[ii]][[rr]] %*%
                               kronecker(Diagonal(pc$h[[hh]]$nplev[ii,rr]),
                                         Sigma_hr2[[rr]]) %*%
                               t(pc$h[[hh]]$Z2[[ii]][[rr]])
                } else {
                  Xi_hi2[[rr]] = Diagonal(n_hi[rr],0)
                }
              } else {
                Xi_hi2[[rr]] = Diagonal(n_hi[rr],0)
              }
              pred_out$h[[hh]]$Xi_2[[ii]][[rr]] = Xi_hi2[[rr]]
            }
          }
        }
      }

      # calculate model covariance matrix
      Omega = Matrix(0,n_hi_tot,n_hi_tot)
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
              Sigma_hi = Matrix(0,n_hi[r1],n_hi[r2])
              Sigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[r1]][[r2]]] =
                Sigma_h[r1,r2]
              Omega[ir,ic] = Sigma_hi
              if( r2 > r1 ){ Omega[ic,ir] = t(Sigma_hi) }
            }
          }
        }
      }
      pred_out$h[[hh]]$Sigma_eps[[ii]] = Omega
      if( any(pc$h[[hh]]$pvc_1 > 0) ){
        pred_out$h[[hh]]$Sigma_vc1[[ii]] = bdiag(Xi_hi1)
        Xi_hi1 = Xi_hi1[!sapply(Xi_hi1,is.null)]
        Omega = Omega + bdiag(Xi_hi1)
      }
      if( any(pc$h[[hh]]$pvc_1 > 0 & pc$h[[hh]]$pvc_2 > 0) ){
        pred_out$h[[hh]]$Sigma_vc2[[ii]] = bdiag(Xi_hi2)
        Xi_hi2 = Xi_hi2[!sapply(Xi_hi2,is.null)]
        Omega = Omega + bdiag(Xi_hi2)
      }
      pred_out$h[[hh]]$Omega[[ii]] = Omega
      cOmega = chol(Omega)
      pred_out$h[[hh]]$IOmega[[ii]] = chol2inv(cOmega)
    }
  }
  return(pred_out)
}
