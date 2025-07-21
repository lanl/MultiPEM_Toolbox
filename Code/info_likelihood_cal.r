########################################################################
#                                                                      #
# This file contains code for calculating the (inverse) information    #
# matrix of the calibration inference parameters.                      #
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

info_ll_cal = function(opt, pc)
{
  # use R Matrix package
  require(Matrix)

  # check that calibration inference parameters exist for at least
  # one phenomenology
  if( !exists("calp",where=opt,inherits=FALSE) ){
    stop(paste("Calibration inference parameters must be present ",
               "for at least one phenomenology.",sep=""))
  }

  # extract calibration inference parameters
  calpHat = opt$calp

  # extract errors-in-variables yield parameters
  if( exists("w_eiv",where=opt,inherits=FALSE) ){ wHat_eiv = opt$w_eiv }

  # extract common forward model parameters
  if( exists("beta",where=opt,inherits=FALSE) ){ betaHat = opt$beta }

  # extract forward model parameters dependent on
  # emplacement condition
  if( exists("tbeta",where=opt,inherits=FALSE) ){ tbetaHat = opt$tbeta }

  # extract source variance components
  if( exists("vc_1",where=opt,inherits=FALSE) ){ vcHat_1 = opt$vc_1 }

  # extract path variance components
  if( exists("vc_2",where=opt,inherits=FALSE) ){ vcHat_2 = opt$vc_2 }

  # extract observational error covariance parameters
  epsHat = opt$eps

  # list containing asymptotic covariance matrix (matrices)
  Cmat = list()
  Cmat$acov_cal = 0

  # information matrices for calibration inference parameters
  I_pp = Matrix(0,pc$ncalp,pc$ncalp,sparse=FALSE,doDiag=FALSE)

  # information matrices for errors-in-variables yields
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    I_wc = Matrix(0,pc$nsource,pc$nsource,sparse=FALSE,doDiag=FALSE)
    I_ww = Matrix(0,pc$nsource,pc$nsource,sparse=FALSE,doDiag=FALSE)
    diag(I_ww) = 1/pc$eiv_w_sd^2
    I_pw = Matrix(0,pc$ncalp,pc$nsource,sparse=FALSE,doDiag=FALSE)
  }

  # list of information matrices for common forward model parameters
  if( pc$pbeta > 0 ){ I_0 = vector("list",pc$H) }

  # list of information matrices for multiple emplacement
  # condition forward model parameters
  if( pc$ptbeta > 0 ){
    I_t = vector("list",pc$H)
    if( pc$pbeta > 0 ){ I_0t = vector("list",pc$H) }
  }

  # list of "cross" information matrices between common forward
  # model coefficients and calibration inference parameters
  if( pc$pbeta > 0 ){ I_0p = vector("list",pc$H) }

  # list of "cross" information matrices between emplacement
  # dependent forward model coefficients and calibration inference
  # parameters
  if( pc$ptbeta > 0 ){ I_tp = vector("list",pc$H) }

  # list of "cross" information matrices between common forward
  # model coefficients and errors-in-variables yields
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    if( pc$pbeta > 0 ){ I_0w = vector("list",pc$H) }
  }

  # list of "cross" information matrices between emplacement
  # dependent forward model coefficients and errors-in-variables
  # yields
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    if( pc$ptbeta > 0 ){ I_tw = vector("list",pc$H) }
  }

  # compute covariance matrices by source
  # compute information matrices by phenomenology
  for( hh in 1:pc$H ){
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
      beta = betaHat[1:pbeta]
      betaHat = betaHat[-(1:pbeta)]
    } else { pbeta = 0 }
    if( pc$ptbeta > 0 && "ptbeta" %in% pnames ){
      ptbeta = sum(pc$h[[hh]]$ptbeta)
      tbeta = tbetaHat[1:ptbeta]
      tbetaHat = tbetaHat[-(1:ptbeta)]
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
          Sigma_hr1[[rr]] = Diagonal(pvc_1,exp(vcHat_1[1:pvc_1]))
          vcHat_1 = vcHat_1[-(1:pvc_1)]
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
          Sigma_hr2[[rr]] = Diagonal(pvc_2,exp(vcHat_2[1:pvc_2]))
          vcHat_2 = vcHat_2[-(1:pvc_2)]
        }
      }
    }
    # construct observational error
    # covariance matrices
    L_h = Diagonal(Rh,exp(epsHat[1:Rh]))
    epsHat = epsHat[-(1:Rh)]
    if( Rh > 1 ){
      L_h[upper.tri(L_h)] = epsHat[1:choose(Rh,2)]
      epsHat = epsHat[-(1:choose(Rh,2))]
    }
    Sigma_h = t(L_h) %*% L_h

    # setup for calibration inference parameters
    pm$ncalp = 0
    if( "cal_par_names" %in% pnames ){
      csub = which(pc$cal_par_names %in% pc$h[[hh]]$cal_par_names)
      if( length(csub) > 0 ){
        Cp = calpHat[csub]
        pm$cal = TRUE
        pm$cal_par_names = pc$h[[hh]]$cal_par_names
        pm$ncalp = length(pm$cal_par_names)
      } else { pm$cal = FALSE }
    } else { pm$cal = FALSE }

    # number of emplacement conditions for phenomenology "hh"
    if( ptbeta > 0 ){ Th = pc$h[[hh]]$Th
    } else { Th = 1 }

    # initialize information matrix for common forward model
    # parameters
    if( pbeta > 0 ){
      I_0[[hh]] = Matrix(0,pbeta,pbeta,sparse=FALSE,doDiag=FALSE)
    }

    # create list of information matrices for multiple emplacement
    # condition forward model parameters
    if( ptbeta > 0 ){
      I_t[[hh]] = vector("list",Th)
      if( pbeta > 0 ){ I_0t[[hh]] = vector("list",Th) }
    }

    # initialize "cross" information matrix between common forward
    # model coefficients and calibration inference parameters
    if( pm$ncalp > 0 ){
      if( pbeta > 0 ){
        I_0p[[hh]] = Matrix(0,pbeta,pc$ncalp,sparse=FALSE,doDiag=FALSE)
      }
    }

    # create list of "cross" information matrices between emplacement
    # dependent forward model coefficients and calibration inference
    # parameters
    if( pm$ncalp > 0 ){
      if( ptbeta > 0 ){ I_tp[[hh]] = vector("list",Th) }
    }

    # initialize "cross" information matrix between common forward
    # model coefficients and errors-in-variables yields
    if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
      if( pbeta > 0 ){
        I_0w[[hh]] = Matrix(0,pbeta,pc$nsource,sparse=FALSE,
                            doDiag=FALSE)
      }
    }

    # create list of "cross" information matrices between emplacement
    # dependent forward model coefficients and errors-in-variables
    # yields
    if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
      if( ptbeta > 0 ){ I_tw[[hh]] = vector("list",Th) }
    }

    for( tt in 1:Th ){
      # extract forward model parameters for emplacement condition
      # "tt"
      if( ptbeta > 0 ){
        st_betat = 0
        if( tt > 1 ){ st_betat = sum(pc$h[[hh]]$ptbeta[1:(tt-1)]) }
        if( pc$h[[hh]]$ptbeta[tt] > 0 ){
          betat = tbeta[st_betat+(1:pc$h[[hh]]$ptbeta[tt])]
        }
        i_source = pc$h[[hh]]$sourceht[[tt]]
      } else { i_source = 1:pc$h[[hh]]$nsource }

      # initialize emplacement condition dependent information matrices
      if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
        I_t[[hh]][[tt]] = Matrix(0,pc$h[[hh]]$ptbeta[tt],
                                 pc$h[[hh]]$ptbeta[tt],sparse=FALSE,
                                 doDiag=FALSE)
        if( pbeta > 0 ){
          I_0t[[hh]][[tt]] = Matrix(0,pbeta,pc$h[[hh]]$ptbeta[tt],
                                    sparse=FALSE,doDiag=FALSE)
        }
      }

      # initialize emplacement condition dependent "cross" information
      # matrices between forward model coefficients and calibration
      # inference parameters
      if( pm$ncalp > 0 ){
        if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
          I_tp[[hh]][[tt]] = Matrix(0,pc$h[[hh]]$ptbeta[tt],pc$ncalp,
                                    sparse=FALSE,doDiag=FALSE)
        }
      }

      # initialize emplacement condition dependent "cross" information
      # matrices between forward model coefficients and
      # errors in variables yields
      if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
        if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
          I_tw[[hh]][[tt]] = Matrix(0,pc$h[[hh]]$ptbeta[tt],pc$nsource,
                                    sparse=FALSE,doDiag=FALSE)
        }
      }

      # iterate over sources of type "tt" in phenomenology "hh"
      for( gg in 1:pc$h[[hh]]$nsource_groups ){
        i_source_g = intersect(i_source,pc$h[[hh]]$Source_Groups[[gg]])
        if( length(i_source_g) == 0 ){ next }
        sc = 0
        tsc = length(i_source_g)
        Jac_c = NULL
        g_w = vector("list",tsc)
        Jac_0 = NULL
        Jac_t = NULL
        Omega = vector("list",tsc)
        for( ii in i_source_g ){
          sc = sc + 1
          # number of responses for source "ii"
          n_hi = pc$h[[hh]]$n[[ii]]
          n_hi_tot = sum(n_hi)

          # setup argument for forward model/jacobian call
          Arg = "(c(beta_t"
          if( pm$cal ){ Arg = paste(Arg,",Cp",sep="") }
          if( "eiv" %in% pnames ){
            if( !is.null(pc$h[[hh]]$eiv[[ii]]) ){
              W = wHat_eiv[pc$h[[hh]]$eiv[[ii]]]
              pm$theta_names = "W"
              Arg = paste(Arg,",W",sep="")
            } else {
              if( exists("theta_names",where=pm,inherits=FALSE) ){
                rm(theta_names, envir=pm)
              }
            }
          }
          Arg=paste(Arg,"),pm)",sep="")

          # named parameters in forward model/jacobian call
          if( "theta_names" %in% pnames &&
              !is.null(pc$h[[hh]]$theta_names[[ii]]) ){
            pm$theta_names = pc$h[[hh]]$theta_names[[ii]]
          }

          # covariance matrices
          if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
            Xi_hi1 = vector("list",Rh)
          }

          # Jacobian matrix for calibration inference parameters
          if( pm$cal ){ Jac_c_r = NULL }

          # gradient vector for errors-in-variables yield
          if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
            g_w_r = NULL
          }

          # Jacobian matrices for forward model parameters
          if( pbeta > 0 ){
            Jac_0_r = Matrix(0,n_hi_tot,pbeta,sparse=FALSE,doDiag=FALSE)
          }
          if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
            Jac_t_r = Matrix(0,n_hi_tot,pc$h[[hh]]$ptbeta[tt],
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
                betatr = betat[st_betatr+
                               (1:pc$h[[hh]]$pbetat[[tt]][rr])]
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

              # calculate forward model Jacobians
              gcall = paste("pc$gfm$",pc$h[[hh]]$g[rr],Arg,sep="")
              jac = eval(parse(text=gcall)) 

              # calculate components of Jacobian matrix
              if( pm$cal ){
                if( is.list(jac) ){ tjac = jac$jcalp
                } else { tjac = jac }
                Jac_c_r = rbind(Jac_c_r,tjac)
              }
              if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
                if( is.list(jac) ){ tjac = jac$jtheta
                } else { tjac = jac }
                g_w_r = c(g_w_r,tjac)
              }
              if( is.list(jac) && exists("jbeta",where=jac,
                  inherits=FALSE) ){ jac = jac$jbeta }
              if( is.null(betatr) && pbeta > 0 &&
                  pc$h[[hh]]$pbeta[rr] > 0 ){
                ir = st_nir+(1:n_hi[rr])
                ic = st_beta+(1:pc$h[[hh]]$pbeta[rr])
                Jac_0_r[ir,ic] = jac
              }
              if( is.null(betar) && ptbeta > 0 &&
                  pc$h[[hh]]$pbetat[[tt]][rr] > 0 ){
                ir = st_nir+(1:n_hi[rr])
                ic = st_betatr+(1:pc$h[[hh]]$pbetat[[tt]][rr])
                Jac_t_r[ir,ic] = jac
              }
              if( !is.null(betar) && !is.null(betatr) ){
                ir = st_nir+(1:n_hi[rr])
                ic = st_beta+(1:pc$h[[hh]]$pbeta[rr])
                Jac_0_r[ir,ic] = jac[,pc$h[[hh]]$ibetar[[(tt-1)*Rh+rr]]]
                ic = st_betatr+(1:pc$h[[hh]]$pbetat[[tt]][rr])
                Jac_t_r[ir,ic] =
                  jac[,pc$h[[hh]]$ibetatr[[(tt-1)*Rh+rr]]]
              }
    
              # calculate components of model covariance matrix
              if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
                if( pc$h[[hh]]$pvc_1[rr] > 0 &&
                    !is.null(pc$h[[hh]]$Z1[[ii]][[rr]]) ){
                  Xi_hi1[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]] %*%
                                 Sigma_hr1[[rr]] %*%
                                 t(pc$h[[hh]]$Z1[[ii]][[rr]])
                } else { Xi_hi1[[rr]] = Diagonal(n_hi[rr],0) }
              }
            }
          }
          if( pm$cal ){ Jac_c = rbind(Jac_c, Jac_c_r) }
          if( "eiv" %in% pnames ){
            if( is.null(pc$h[[hh]]$eiv[[ii]]) ){
              for( isc in 1:tsc ){
                g_w[[isc]] = c(g_w[[isc]], numeric(n_hi_tot))
              }
            } else {
              g_w[[sc]] = c(g_w[[sc]], g_w_r)
              for( isc in setdiff(1:tsc,sc) ){
                g_w[[isc]] = c(g_w[[isc]], numeric(n_hi_tot))
              }
            }
          }
          if( pbeta > 0 ){ Jac_0 = rbind(Jac_0, Jac_0_r) }
          if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
            Jac_t = rbind(Jac_t, Jac_t_r)
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
        Catch = pc$tryCatch.W.E(chol(Omega))
        if( is(Catch$value,"Matrix") ){ cOmega = Catch$value
        } else {
          cat("Omega is not positive definite.\n")
          Omega = nearPD(Omega)$mat
          cOmega = chol(Omega)
        }
        IOmega = chol2inv(cOmega)

        # information matrix for common forward model coefficients
        if( pbeta > 0 ){
          I_0[[hh]] = I_0[[hh]] + t(Jac_0) %*% IOmega %*% Jac_0
        }

        # information matrix for emplacement condition dependent
        # forward model coefficients
        if( ptbeta > 0 && !is.null(I_t[[hh]][[tt]]) ){
          I_t[[hh]][[tt]] = I_t[[hh]][[tt]] +
                            t(Jac_t) %*% IOmega %*% Jac_t
        }

        # "cross" information matrix between common and emplacement
        # condition dependent forward model coefficients
        if( pbeta > 0 && ptbeta > 0 && !is.null(I_0t[[hh]][[tt]]) ){
          I_0t[[hh]][[tt]] = I_0t[[hh]][[tt]] +
                             t(Jac_0) %*% IOmega %*% Jac_t
        }

        # compute components of information matrices needed for
        # calibration inference parameters
        if( pm$cal ){
          I_pp[csub,csub] = I_pp[csub,csub] +
            t(Jac_c) %*% IOmega %*% Jac_c
          if( pbeta > 0 ){
            I_0p[[hh]][,csub] = I_0p[[hh]][,csub] +
              t(Jac_0) %*% IOmega %*% Jac_c
            if( ptbeta > 0 && !is.null(I_tp[[hh]][[tt]]) ){
              I_tp[[hh]][[tt]][,csub] = I_tp[[hh]][[tt]][,csub] +
                t(Jac_t) %*% IOmega %*% Jac_c
            }
          } else if( ptbeta > 0 && !is.null(I_tp[[hh]][[tt]]) ){
            I_tp[[hh]][[tt]][,csub] = I_tp[[hh]][[tt]][,csub] +
              t(Jac_t) %*% IOmega %*% Jac_c
          }
        }

        # compute components of information matrices needed for
        # errors-in-variables correction
        sc_r = 0
        for( ii in i_source_g ){
          sc_r = sc_r + 1
          if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
            sc_c = 0
            for( jj in i_source_g ){
              sc_c = sc_c + 1
              if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[jj]]) ){
                I_wc[pc$h[[hh]]$eiv[[ii]],pc$h[[hh]]$eiv[[jj]]] =
                I_wc[pc$h[[hh]]$eiv[[ii]],pc$h[[hh]]$eiv[[jj]]] +
                t(g_w[[sc_r]]) %*% IOmega %*% g_w[[sc_c]]
              }
            }
            if( pbeta > 0 ){
              I_0w[[hh]][,pc$h[[hh]]$eiv[[ii]]] =
                t(Jac_0) %*% IOmega %*% g_w[[sc_r]]
              if( ptbeta > 0 && !is.null(I_tw[[hh]][[tt]]) ){
                I_tw[[hh]][[tt]][,pc$h[[hh]]$eiv[[ii]]] =
                  t(Jac_t) %*% IOmega %*% g_w[[sc_r]]
              }
            } else if( ptbeta > 0 && !is.null(I_tw[[hh]][[tt]]) ){
              I_tw[[hh]][[tt]][,pc$h[[hh]]$eiv[[ii]]] =
                t(Jac_t) %*% IOmega %*% g_w[[sc_r]]
            }
            if( pm$cal ){
              I_pw[csub,pc$h[[hh]]$eiv[[ii]]] =
                I_pw[csub,pc$h[[hh]]$eiv[[ii]]] +
                t(Jac_c) %*% IOmega %*% g_w[[sc_r]]
            }
          }
        }
      }
    }
  }

  # construct phenomenology-specific information matrices
  if( pc$pbeta > 0 || pc$ptbeta > 0 ){
    I_bb = vector("list",pc$H)
    I_bc = vector("list",pc$H)
    if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
      I_bw = vector("list",pc$H)
    }
  }
  for( hh in 1:pc$H ){
    if( pc$pbeta > 0 && any(pc$h[[hh]]$pbeta > 0) ){
      pbeta = sum(pc$h[[hh]]$pbeta)
    } else { pbeta = 0 }
    pnames = names(pc$h[[hh]])
    if( pc$ptbeta > 0 && "ptbeta" %in% pnames ){
      ptbeta = sum(pc$h[[hh]]$ptbeta)
    } else { ptbeta = 0 }
    if( pbeta > 0 ){ I_bb[[hh]] = I_0[[hh]] }
    if( ptbeta > 0 ){
      D_t = bdiag(I_t[[hh]][!sapply(I_t[[hh]],is.null)])
      if( pbeta > 0 ){
        D_0t = do.call(cbind,I_0t[[hh]])
        I_bb[[hh]] = cbind(I_bb[[hh]],D_0t)
        I_bb[[hh]] = rbind(I_bb[[hh]],cbind(t(D_0t),D_t))
      } else { I_bb[[hh]] = D_t }
    }
    if( length(pc$h[[hh]]$cal_par_names) > 0 ){
      if( pbeta > 0 ){ I_bc[[hh]] = I_0p[[hh]] }
      if( ptbeta > 0 ){
        I_bc[[hh]] = rbind(I_bc[[hh]],do.call(rbind,I_tp[[hh]]))
      }
    } else {
      mbeta = pbeta+ptbeta
      if( mbeta > 0 ){
        I_bc[[hh]] = Matrix(0,mbeta,pc$ncalp,sparse=TRUE,
                            doDiag=FALSE)
      }
    }
    if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
      if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
        if( pbeta > 0 ){ I_bw[[hh]] = I_0w[[hh]] }
        if( ptbeta > 0 ){
          I_bw[[hh]] = rbind(I_bw[[hh]],do.call(rbind,I_tw[[hh]]))
        }
      } else {
        mbeta = pbeta+ptbeta
        if( mbeta > 0 ){
          I_bw[[hh]] = Matrix(0,mbeta,pc$nsource,sparse=TRUE,
                              doDiag=FALSE)
        }
      }
    }
  }

  # construct global information matrices
  I_pp = forceSymmetric(I_pp)
  Catch = pc$tryCatch.W.E(chol(I_pp))
  if( is(Catch$value,"Matrix") ){ C_pp = Catch$value
  } else {
    cat("I_c(calp) is not positive definite.\n")
    I_pp = nearPD(I_pp)$mat
    C_pp = chol(I_pp)
  }
  II_pp = chol2inv(C_pp)
  S_pp = I_pp
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    S_ww = I_wc + I_ww
    S_cw = I_pw
  }
  if( pc$pbeta > 0 || pc$ptbeta > 0 ){
    S_bb = bdiag(I_bb[!sapply(I_bb,is.null)])
    S_bb = forceSymmetric(S_bb)
    Catch = pc$tryCatch.W.E(chol(S_bb))
    if( is(Catch$value,"Matrix") ){ C_bb = Catch$value
    } else {
      cat("I(beta) is not positive definite.\n")
      S_bb = nearPD(S_bb)$mat
      C_bb = chol(S_bb)
    }
    IS_bb = chol2inv(C_bb)
    S_bc = do.call(rbind,I_bc)
    T_cb = t(S_bc) %*% IS_bb
    S_pp = S_pp - T_cb %*% S_bc
    if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
      S_bw = do.call(rbind,I_bw)
      S_ww = S_ww - t(S_bw) %*% IS_bb %*% S_bw
      S_cw = S_cw - T_cb %*% S_bw
    }
  }
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    S_ww = forceSymmetric(S_ww)
    Catch = pc$tryCatch.W.E(chol(S_ww))
    if( is(Catch$value,"Matrix") ){ C_ww = Catch$value
    } else {
      if( pc$pbeta > 0 || pc$ptbeta > 0 ){
        cat("Adjusted I(eiv) is not positive definite.\n")
      } else {
        cat("I(eiv) is not positive definite.\n")
      }
      S_ww = nearPD(S_ww)$mat
      C_ww = chol(S_ww)
    }
    IS_ww = chol2inv(C_ww)
    S_pp = S_pp - S_cw %*% IS_ww %*% t(S_cw)
  }
  if( pc$pbeta > 0 || pc$ptbeta > 0 ||
      (exists("eiv",where=pc,inherits=FALSE) && pc$eiv) ){
    S_pp = forceSymmetric(S_pp)
    Catch = pc$tryCatch.W.E(chol(S_pp))
    if( is(Catch$value,"Matrix") ){ C_pp = Catch$value
    } else {
      cat("Adjusted I(calp) is not positive definite.\n")
      S_pp = nearPD(S_pp)$mat
      C_pp = chol(S_pp)
    }
    II_pp = chol2inv(C_pp)
  }
  Cmat$II_calp = II_pp
  Cmat$acov_cal = 1

  return(Cmat)
}
