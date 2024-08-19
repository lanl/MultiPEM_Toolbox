########################################################################
#                                                                      #
# This file contains code for calculating the gradient of the log-     #
# likelihood of the calibration parameters, with an option to include  #
# errors-in-variables for calibration data yields with their           #
# gradients. New event data can also be optionally incorporated,       #
# allowing for joint log-likelihood gradients of new event and         #
# calibration parameters based on calibration and new event data.      #
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

gll_cal = function(x, pc)
{
  # use R Matrix package
  require(Matrix)

  # extract inference parameters for new event
  if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
    theta0 = x[1:pc$ntheta0]
    if( exists("itransform",where=pc,inherits=FALSE) ){
      if( pc$itransform ){
        dtheta0 = pc$j_tau(theta0, pc=pc)
        theta0 = pc$tau(theta0, pc=pc)
      }
    }
    if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
      ith0_bds = pc$itheta0_bounds
      dtheta0_b = rep(1,pc$ntheta0)
      if( length(ith0_bds[[1]]) > 0 ){
        dtheta0_b[ith0_bds[[1]]] = pc$dnotExp(theta0[ith0_bds[[1]]])
      }
      if( length(ith0_bds[[2]]) > 0 ){
        dtheta0_b[ith0_bds[[2]]] = -pc$dnotExp(theta0[ith0_bds[[2]]])
      }
      if( length(ith0_bds[[3]]) > 0 ){
        tau = pc$notExp(theta0[ith0_bds[[3]]])
        dtheta0_b[ith0_bds[[3]]] = pc$dnotExp(theta0[ith0_bds[[3]]])*
                                   pc$theta0_range/(1+tau)^2
      }
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

  # log-likelihood gradient due to calibration data and
  # (optionally) new event data

  # initialize gradient vectors
  # calibration inference parameters
  if( pc$ncalp > 0 ){ gr_cp = numeric(pc$ncalp)
  } else { gr_cp = NULL }
  # new event parameters
  if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
    gr_th0 = numeric(pc$ntheta0)
  } else { gr_th0 = NULL }
  # errors-in-variables
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    gr_eiv = numeric(pc$nsource)
  } else { gr_eiv = NULL }
  # common forward model parameters
  gr_beta0 = NULL
  # emplacement condition dependent forward model parameters
  gr_betat = NULL
  # level 1 variance components
  gr_vc1 = NULL
  # level 2 variance components
  gr_vc2 = NULL
  # observational error covariance parameters
  gr_eps = NULL

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
    if( exists("dnotExp",where=pc,inherits=FALSE) ){
      pm$dnotExp = pc$dnotExp
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

    # number of emplacement conditions for phenomenology "hh"
    if( ptbeta > 0 ){ Th = pc$h[[hh]]$Th }

    # initialize gradient vectors forward model parameters
    if( pbeta > 0 ){ g_beta0 = numeric(pbeta) }
    if( ptbeta > 0 ){
      g_betat = vector("list",Th)
      for( tt in 1:Th ){
        if( pc$h[[hh]]$ptbeta[tt] > 0 ){
          g_betat[[tt]] = numeric(pc$h[[hh]]$ptbeta[tt])
        }
      }
    }
    # initialize gradient vectors variance components
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      g_vc1 = NULL
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                              pc$h[[hh]]$pvc_2 > 0) ){
        g_vc2 = NULL
      }
      for( rr in 1:Rh ){
        if( pc$h[[hh]]$pvc_1[rr] > 0 ){
          g_vc1 = c(g_vc1,numeric(pc$h[[hh]]$pvc_1[rr]))
          if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                  pc$h[[hh]]$pvc_2 > 0) ){
            if( pc$h[[hh]]$pvc_2[rr] > 0 ){
              g_vc2 = c(g_vc2,numeric(pc$h[[hh]]$pvc_2[rr]))
            }
          }
        }
      }
    }
    g_eps = numeric(Rh*(Rh+1)/2)

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
    for(ii in 1:pc$h[[hh]]$nsource){
      # number of responses for source "ii"
      n_hi = pc$h[[hh]]$n[[ii]]
      n_hi_tot = sum(n_hi)

      # setup argument for forward model/jacobian call
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

      # named parameters in forward model/jacobian call
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

      # residual vector
      resid = NULL

      # covariance matrices and differentials
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        Xi_hi1 = vector("list",Rh)
        dXi_hi1 = vector("list",Rh)
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          Xi_hi2 = vector("list",Rh)
          dXi_hi2 = vector("list",Rh)
        }
      }

      # Jacobian matrix for new event inference parameters
      if( "nev" %in% pnames && pc$h[[hh]]$nev[ii] ){ Jac_th0 = NULL }

      # Jacobian matrix for calibration inference parameters
      if( pm$cal ){ Jac_c = NULL }

      # gradient vector for errors-in-variables yield
      if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
        g_w = NULL
      }

      # Jacobian matrices for forward model parameters
      if( pbeta > 0 ){
        Jac_0 = Matrix(0,n_hi_tot,pbeta,sparse=FALSE,doDiag=FALSE)
      }
      if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
        Jac_t = Matrix(0,n_hi_tot,pc$h[[hh]]$ptbeta[tt],
                       sparse=FALSE,doDiag=FALSE)
      }

      # iterate over responses "rr"
      for(rr in 1:Rh){
        if( n_hi[rr] > 0 ){
          st_nir = 0
          if( rr > 1 ){ st_nir = sum(n_hi[1:(rr-1)]) }

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
          if( any(is.nan(yhat)) ){ return(NaN) }

          # calculate residual vector
          resid = c(resid,pc$h[[hh]]$Y[[ii]][[rr]] - yhat)

          # calculate components of Jacobian matrix
          gcall = paste("pc$gfm$",pc$h[[hh]]$g[rr],Arg,sep="")
          jac = eval(parse(text=gcall))
          if( "nev" %in% pnames && pc$h[[hh]]$nev[ii] ){
            Jac_th0 = rbind(Jac_th0,jac$jtheta)
          }
          if( pm$cal ){
            Jac_c = rbind(Jac_c,jac$jcalp)
          }
          if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
            g_w = c(g_w,jac$jtheta)
          }
          if( is.list(jac) ){ jac = jac$jbeta }
          if( is.null(betatr) && pc$h[[hh]]$pbeta[rr] > 0 ){
            ir = st_nir+(1:n_hi[rr])
            ic = st_beta+(1:pc$h[[hh]]$pbeta[rr])
            Jac_0[ir,ic] = jac
          }
          if( is.null(betar) && pc$h[[hh]]$pbetat[[tt]][rr] > 0 ){
            ir = st_nir+(1:n_hi[rr])
            ic = st_betatr+(1:pc$h[[hh]]$pbetat[[tt]][rr])
            Jac_t[ir,ic] = jac
          }
          if( !is.null(betar) && !is.null(betatr) ){
            ir = st_nir+(1:n_hi[rr])
            ic = st_beta+(1:pc$h[[hh]]$pbeta[rr])
            Jac_0[ir,ic] = jac[,pc$h[[hh]]$ibetar[[(tt-1)*Rh+rr]]]
            ic = st_betatr+(1:pc$h[[hh]]$pbetat[[tt]][rr])
            Jac_t[ir,ic] = jac[,pc$h[[hh]]$ibetatr[[(tt-1)*Rh+rr]]]
          }

          # calculate components of model covariance matrix
          if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
            if( pc$h[[hh]]$pvc_1[rr] > 0 ){
              Xi_hi1[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]] %*%
                             Sigma_hr1[[rr]] %*%
                             t(pc$h[[hh]]$Z1[[ii]][[rr]])
              dXi_hi1[[rr]] = vector("list",pc$h[[hh]]$pvc_1[rr])
              for( ll in 1:pc$h[[hh]]$pvc_1[rr] ){
                dXi_hi1[[rr]][[ll]] = (Sigma_hr1[[rr]][ll,ll] *
                                     pc$h[[hh]]$Z1[[ii]][[rr]][,ll]) %*%
                                     t(pc$h[[hh]]$Z1[[ii]][[rr]][,ll])
              }
            } else {
              Xi_hi1[[rr]] = Diagonal(n_hi[rr],0)
            }
            if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                    pc$h[[hh]]$pvc_2 > 0) ){
              if( pc$h[[hh]]$pvc_1[rr] > 0 ){
                if( pc$h[[hh]]$pvc_2[rr] > 0 ){
                  Xi_hi2[[rr]] = pc$h[[hh]]$Z2[[ii]][[rr]] %*%
                               kronecker(Diagonal(pc$h[[hh]]$nplev[ii,rr]),
                                         Sigma_hr2[[rr]]) %*%
                               t(pc$h[[hh]]$Z2[[ii]][[rr]])
                  dXi_hi2[[rr]] = vector("list",pc$h[[hh]]$pvc_2[rr])
                  for( ll in 1:pc$h[[hh]]$pvc_2[rr] ){
                    el = Matrix(0,pc$h[[hh]]$pvc_2[rr],1,sparse=FALSE,
                                doDiag=FALSE)
                    el[ll] = sqrt(Sigma_hr2[[rr]][ll,ll])
                    dXi_hi2[[rr]][[ll]] = pc$h[[hh]]$Z2[[ii]][[rr]] %*%
                               kronecker(Diagonal(pc$h[[hh]]$nplev[ii,rr]),
                                         el %*% t(el)) %*%
                               t(pc$h[[hh]]$Z2[[ii]][[rr]])
                  }
                } else {
                  Xi_hi2[[rr]] = Diagonal(n_hi[rr],0)
                }
              } else {
                Xi_hi2[[rr]] = Diagonal(n_hi[rr],0)
              }
            }
          }
        }
      }

      # calculate model covariance matrix
      Omega = Matrix(0,n_hi_tot,n_hi_tot,sparse=FALSE,doDiag=FALSE)
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
              Omega[ir,ic] = Sigma_hi
              if( r2 > r1 ){ Omega[ic,ir] = t(Sigma_hi) }
            }
          }
        }
      }
      if( any(pc$h[[hh]]$pvc_1 > 0) ){
        Xi_hi1 = Xi_hi1[!sapply(Xi_hi1,is.null)]
        Omega = Omega + bdiag(Xi_hi1)
      }
      if( any(pc$h[[hh]]$pvc_1 > 0 & pc$h[[hh]]$pvc_2 > 0) ){
        Xi_hi2 = Xi_hi2[!sapply(Xi_hi2,is.null)]
        Omega = Omega + bdiag(Xi_hi2)
      }
      cOmega = chol(Omega)
      IOmega = chol2inv(cOmega)
      resid_io = IOmega %*% resid

      # calculate components of log-likelihood gradient
      # new event inference parameters
      if( "nev" %in% pnames && pc$h[[hh]]$nev[ii] ){
        g_th0 = t(Jac_th0) %*% resid_io
        if( "itheta0" %in% pnames ){
          if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
            g_th0 = g_th0 * dtheta0_b[pc$h[[hh]]$itheta0]
          }
          if( exists("itransform",where=pc,inherits=FALSE) ){
            if( pc$itransform ){
              g_th0 =
                t(dtheta0[pc$h[[hh]]$itheta0,pc$h[[hh]]$itheta0]) %*%
                g_th0
            }
          }
          gr_th0[pc$h[[hh]]$itheta0] = gr_th0[pc$h[[hh]]$itheta0]+
                                       g_th0
        } else {
          if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
            g_th0 = g_th0 * dtheta0_b
          }
          if( exists("itransform",where=pc,inherits=FALSE) ){
            if( pc$itransform ){ g_th0 = t(dtheta0) %*% g_th0 }
          }
          gr_th0 = gr_th0 + g_th0
        }
      }
      # calibration inference parameters
      if( pm$cal ){
        gr_cp[csub] = gr_cp[csub] + t(Jac_c) %*% resid_io
      }
      # errors-in-variables
      if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
        gr_eiv[pc$h[[hh]]$eiv[[ii]]] = gr_eiv[pc$h[[hh]]$eiv[[ii]]]+
                                       t(g_w) %*% resid_io
      }
      # common forward model parameters
      if( pbeta > 0 ){
        g_beta0 = g_beta0+t(Jac_0) %*% resid_io
      }
      # emplacement condition dependent forward model parameters
      if( ptbeta > 0 && !is.null(g_betat[[tt]]) ){
        g_betat[[tt]] = g_betat[[tt]]+t(Jac_t) %*% resid_io
      }
      # matrix utilized in variance component gradients
      IOmegaAdj = IOmega - resid_io %*% t(resid_io)
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        # zero matrix for variance component gradients
        Zero = Matrix(0,n_hi_tot,n_hi_tot,sparse=FALSE,doDiag=FALSE)
        # starting indices for accessing gradient vector for
        # response "rr"
        st_vc1 = 0
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          st_vc2 = 0
        }
        for( rr in 1:Rh ){
          if( n_hi[rr] > 0 ){
            if( pc$h[[hh]]$pvc_1[rr] > 0 ){
              if( rr > 1 ){ st_vc1 = sum(pc$h[[hh]]$pvc_1[1:(rr-1)]) }
              # starting indices for accessing response "rr"
              # contributions
              st_nir = 0
              if( rr > 1 ){ st_nir = sum(n_hi[1:(rr-1)]) }
              ir = st_nir+(1:n_hi[rr])
              # level 1 variance component gradients
              for( ll in 1:pc$h[[hh]]$pvc_1[rr] ){
                dOmega = Zero
                dOmega[ir,ir] = dXi_hi1[[rr]][[ll]]
                g_vc1[ll+st_vc1] = g_vc1[ll+st_vc1]+
                                   sum(diag(IOmegaAdj %*% dOmega))
              }
              if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                      pc$h[[hh]]$pvc_2 > 0) ){
                if( pc$h[[hh]]$pvc_2[rr] > 0 ){
                  if( rr > 1 ){ st_vc2 =
                                  sum(pc$h[[hh]]$pvc_2[1:(rr-1)]) }
                  # level 2 variance component gradients
                  for( ll in 1:pc$h[[hh]]$pvc_2[rr] ){
                    dOmega = Zero
                    dOmega[ir,ir] = dXi_hi2[[rr]][[ll]]
                    g_vc2[ll+st_vc2] = g_vc2[ll+st_vc2]+
                                       sum(diag(IOmegaAdj %*% dOmega))
                  }
                }
              }
            }
          }
        }
      }
      # observation error parameter gradients
      for( rr in 1:Rh ){
        if( n_hi[rr] > 0 ){
          dOmega = Matrix(0,n_hi_tot,n_hi_tot,sparse=FALSE,doDiag=FALSE)
          st_nir = 0
          if( rr > 1 ){ st_nir = sum(n_hi[1:(rr-1)]) }
          ir = st_nir+(1:n_hi[rr])
          for( r2 in rr:Rh ){
            if( n_hi[r2] > 0 ){
              st_nir2 = 0
              if( r2 > 1 ){ st_nir2 = sum(n_hi[1:(r2-1)]) }
              ic = st_nir2+(1:n_hi[r2])
              dSigma_hi = Matrix(0,n_hi[rr],n_hi[r2],sparse=FALSE,
                                 doDiag=FALSE)
              if( r2 == rr ){
                dSigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[rr]][[rr]]] =
                  2*L_h[rr,rr]^2
                dOmega[ir,ic] = dSigma_hi
              } else {
                dSigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[rr]][[r2]]] =
                  L_h[rr,rr]*L_h[rr,r2]
                dOmega[ir,ic] = dSigma_hi
                dOmega[ic,ir] = t(dSigma_hi)
              }
            }
          }
          g_eps[rr] = g_eps[rr]+sum(diag(IOmegaAdj %*% dOmega))
        }
      }
      if( Rh > 1 ){
        kk = 0
        for( ss in 2:Rh ){
          if( n_hi[ss] > 0 ){
            st_nis = 0
            if( ss > 1 ){ st_nis = sum(n_hi[1:(ss-1)]) }
            is = st_nis+(1:n_hi[ss])
            for( rr in 1:(ss-1) ){
              if( n_hi[rr] > 0 ){
                kk = kk+1
                dOmega = Matrix(0,n_hi_tot,n_hi_tot,sparse=FALSE,
                                doDiag=FALSE)
                for( r1 in rr:(ss-1) ){
                  if( n_hi[r1] > 0 ){
                    st_nir1 = 0
                    if( r1 > 1 ){ st_nir1 = sum(n_hi[1:(r1-1)]) }
                    ir = st_nir1+(1:n_hi[r1])
                    dSigma_hi = Matrix(0,n_hi[r1],n_hi[ss],sparse=FALSE,
                                       doDiag=FALSE)
                    dSigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[r1]][[ss]]]=
                      L_h[rr,r1]
                    dOmega[ir,is] = dSigma_hi
                    dOmega[is,ir] = t(dSigma_hi)
                  }
                }
                for( r2 in ss:Rh ){
                  if( n_hi[r2] > 0 ){
                    st_nir2 = 0
                    if( r2 > 1 ){ st_nir2 = sum(n_hi[1:(r2-1)]) }
                    ir = st_nir2+(1:n_hi[r2])
                    dSigma_hi = Matrix(0,n_hi[ss],n_hi[r2],sparse=FALSE,
                                       doDiag=FALSE)
                    if( r2 == ss ){
                      dSigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[ss]][[ss]]]=
                        2*L_h[rr,ss]
                      dOmega[is,ir] = dSigma_hi
                    } else {
                      dSigma_hi[pc$h[[hh]]$i[[ii]]$cov_pairs[[ss]][[r2]]] =
                        L_h[rr,r2]
                      dOmega[is,ir] = dSigma_hi
                      dOmega[ir,is] = t(dSigma_hi)
                    }
                  }
                }
                g_eps[Rh+kk] = g_eps[Rh+kk]+sum(diag(IOmegaAdj %*%
                                                     dOmega))
              }
            }
          }
        }
      }
    }

    # concatenate gradients in correct order
    gr_th0 = as.numeric(gr_th0)
    gr_cp = as.numeric(gr_cp)
    if( pbeta > 0 ){ gr_beta0 = c(gr_beta0,as.numeric(g_beta0)) }
    if( ptbeta > 0 ){
      for( tt in 1:Th ){
        if( !is.null(g_betat[[tt]]) ){
          gr_betat = c(gr_betat,as.numeric(g_betat[[tt]]))
        }
      }
    }
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      gr_vc1 = c(gr_vc1,-0.5*g_vc1)
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                              pc$h[[hh]]$pvc_2 > 0) ){
        gr_vc2 = c(gr_vc2,-0.5*g_vc2)
      }
    }
    gr_eps = c(gr_eps,-0.5*g_eps)
  }
  return(c(gr_th0,gr_cp,gr_eiv,gr_beta0,gr_betat,gr_vc1,gr_vc2,gr_eps))
}
