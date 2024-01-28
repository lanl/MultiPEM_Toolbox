########################################################################
#                                                                      #
# This file contains code for calculating the (inverse) information    #
# matrix of the new event inference parameters.                        #
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

info_ll = function(opt, pc)
{
  # use R Matrix package
  require(Matrix)

  # check that new event data exists for at least one
  # phenomenology
  if( !exists("nev",where=pc,inherits=FALSE) ){
    stop(paste("New event data must be present ",
               "for all phenomenologies.",sep=""))
  } else {
    if( !pc$nev ){
      stop(paste("New event data must be present ",
                 "for all phenomenologies.",sep=""))
    }
  }

  # extract inference parameters for new event
  theta0Hat = opt$theta0
  if( exists("itransform",where=pc,inherits=FALSE) ){
    if( pc$itransform ){ theta0Hat = pc$tau(theta0Hat, pc=pc) }
  }
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    theta0Hat = pc$transform(theta0Hat, pc=pc)
  }

  # extract errors-in-variables yield parameters
  if( exists("w_eiv",where=opt,inherits=FALSE) ){ wHat_eiv = opt$w_eiv }

  # extract common forward model parameters
  if( exists("beta",where=opt,inherits=FALSE) ){ betaHat = opt$beta }

  # extract forward model parameters dependent on
  # emplacement condition
  if( exists("tbeta",where=opt,inherits=FALSE) ){ tbetaHat = opt$tbeta }

  # extract level 1 variance components
  if( exists("vc_1",where=opt,inherits=FALSE) ){ vcHat_1 = opt$vc_1 }

  # extract level 2 variance components
  if( exists("vc_2",where=opt,inherits=FALSE) ){ vcHat_2 = opt$vc_2 }

  # extract observational error covariance parameters
  epsHat = opt$eps

  # list containing asymptotic covariance matrix (matrices)
  Cmat = list()
  Cmat$acov = 0

  # new event inference parameter information matrix (assuming fixed
  # calibration data yields)
  I_nev = Matrix(0,pc$ntheta0,pc$ntheta0,sparse=FALSE)
  if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
    I_nev_0 = Matrix(0,pc$ntheta0,pc$ntheta0,sparse=FALSE)
  }

  # information matrices for errors-in-variables yields
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    I_wc = Matrix(0,pc$nsource,pc$nsource,sparse=FALSE)
    I_ww = Matrix(0,pc$nsource,pc$nsource,sparse=FALSE)
    diag(I_ww) = 1/pc$eiv_w_sd^2
    Q_cw = Matrix(0,pc$nsource,pc$nsource,sparse=FALSE)
    B_cw = Matrix(0,pc$nsource,pc$ntheta0,sparse=FALSE)
  }

  # compute covariance matrices by source
  # compute information matrices by phenomenology
  for(hh in 1:pc$H){
    # named objects for phenomenology "hh"
    pnames = names(pc$h[[hh]])

    # check that new event data exists for each phenomenology
    if( !("nev" %in% pnames) ){
      stop(paste("New event data must be present ",
                 "for all phenomenologies.",sep=""))
    }

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
          Sigma_hr1[[rr]] = Diagonal(pvc_1,exp(vcHat_1[1:pvc_1]))
          vcHat_1 = vcHat_1[-(1:pvc_1)]
          # construct level 2 variance component
          # covariance matrices
          if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                  pc$h[[hh]]$pvc_2 > 0) ){
            pvc_2 = pc$h[[hh]]$pvc_2[rr]
            if( pvc_2 > 0 ){
              Sigma_hr2[[rr]] = Diagonal(pvc_2,exp(vcHat_2[1:pvc_2]))
              vcHat_2 = vcHat_2[-(1:pvc_2)]
            }
          }
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

    # number of emplacement conditions for phenomenology "hh"
    if( ptbeta > 0 ){ Th = pc$h[[hh]]$Th
    } else { Th = 1 }

    # initialize information matrix for common forward model
    # parameters
    if( pbeta > 0 ){ I_0 = Matrix(0,pbeta,pbeta,sparse=FALSE) }

    # create list of information matrices for multiple emplacement
    # condition forward model parameters
    if( ptbeta > 0 ){
      I_t = vector("list",Th)
      I_0t = vector("list",Th)
    }

    # initialize "cross" information matrix between common forward
    # model coefficients and errors-in-variables yields
    if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
      if( pbeta > 0 ){ I_0w = Matrix(0,pbeta,pc$nsource,sparse=FALSE) }
    }

    # create list of "cross" information matrices between emplacement
    # dependent forward model coefficients and errors-in-variables
    # yields
    if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
      if( ptbeta > 0 ){ I_tw = vector("list",Th) }
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
        I_t[[tt]] = Matrix(0,pc$h[[hh]]$ptbeta[tt],
                           pc$h[[hh]]$ptbeta[tt],sparse=FALSE)
        I_0t[[tt]] = Matrix(0,pbeta,pc$h[[hh]]$ptbeta[tt],
                            sparse=FALSE)
      }

      # initialize emplacement condition dependent "cross" information
      # matrices between forward model coefficients and
      # errors in variables yields
      if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
        if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
          I_tw[[tt]] = Matrix(0,pc$h[[hh]]$ptbeta[tt],pc$nsource,
                              sparse=FALSE)
        }
      }

      for( ii in i_source ){
        # number of responses for source "ii"
        n_hi = pc$h[[hh]]$n[[ii]]
        n_hi_tot = sum(n_hi)

        # setup argument for forward model/jacobian call
        Arg = "(beta_t,pm)"
        if( "eiv" %in% pnames ){
          if( !is.null(pc$h[[hh]]$eiv[[ii]]) ){
            W = wHat_eiv[pc$h[[hh]]$eiv[[ii]]]
            pm$theta_names = "W"
            Arg = "(c(beta_t,W),pm)"
          } else {
            if( exists("theta_names",where=pm,inherits=FALSE) ){
              rm(theta_names, envir=pm)
            }
          }
        }
        if( pc$h[[hh]]$nev[ii] ){
          if( "itheta0" %in% pnames ){
            Arg = "(c(beta_t,theta0Hat[pc$h[[hh]]$itheta0]),pm)"
          } else { Arg = "(c(beta_t,theta0Hat),pm)" }
        }

        # named parameters in forward model/jacobian call
        if( "theta_names" %in% pnames &&
            !is.null(pc$h[[hh]]$theta_names[[ii]]) ){
          pm$theta_names = pc$h[[hh]]$theta_names[[ii]]
        }

        # covariance matrices
        if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
          Xi_hi1 = vector("list",Rh)
          if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                  pc$h[[hh]]$pvc_2 > 0) ){
            Xi_hi2 = vector("list",Rh)
          }
        }

        # Jacobian matrix for new event inference parameters
        if( pc$h[[hh]]$nev[ii] ){ Jac_th0 = NULL }

        # gradient vector for errors-in-variables yield
        if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
          g_w = NULL
        }

        # Jacobian matrices for forward model parameters
        if( pbeta > 0 ){
          Jac_0 = Matrix(0,n_hi_tot,pbeta,sparse=FALSE)
        }
        if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
          Jac_t = Matrix(0,n_hi_tot,pc$h[[hh]]$ptbeta[tt],
                         sparse=FALSE)
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

            # calculate forward model Jacobians
            gcall = paste("pc$gfm$",pc$h[[hh]]$g[rr],Arg,sep="")
            jac = eval(parse(text=gcall)) 

            # calculate components of Jacobian matrix
            if( pc$h[[hh]]$nev[ii] ){
              Jac_th0 = rbind(Jac_th0,jac$jtheta)
              jac = jac$jbeta
            }
            if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
              g_w = c(g_w,jac$jtheta)
              jac = jac$jbeta
            }
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
              } else { Xi_hi1[[rr]] = Diagonal(n_hi[rr],0) }
              if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                      pc$h[[hh]]$pvc_2 > 0) ){
                if( pc$h[[hh]]$pvc_1[rr] > 0 ){
                  if( pc$h[[hh]]$pvc_2[rr] > 0 ){
                    Xi_hi2[[rr]] = pc$h[[hh]]$Z2[[ii]][[rr]] %*%
                            kronecker(Diagonal(pc$h[[hh]]$nplev[ii,rr]),
                                      Sigma_hr2[[rr]]) %*%
                            t(pc$h[[hh]]$Z2[[ii]][[rr]])
                  } else { Xi_hi2[[rr]] = Diagonal(n_hi[rr],0) }
                } else { Xi_hi2[[rr]] = Diagonal(n_hi[rr],0) }
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
        if( any(pc$h[[hh]]$pvc_1 > 0) ){
          Xi_hi1 = Xi_hi1[!sapply(Xi_hi1,is.null)]
          Omega = Omega + bdiag(Xi_hi1)
        }
        if( any(pc$h[[hh]]$pvc_1 > 0 & pc$h[[hh]]$pvc_2 > 0) ){
          Xi_hi2 = Xi_hi2[!sapply(Xi_hi2,is.null)]
          Omega = Omega + bdiag(Xi_hi2)
        }
        Catch = pc$tryCatch.W.E(chol(Omega))
        if( is(Catch$value,"Matrix") ){ cOmega = Catch$value
        } else {
          cat("Omega is not positive definite.\n")
          return(Cmat)
        }
        IOmega = chol2inv(cOmega)

        # information matrix for common forward model coefficients
        if( pbeta > 0 ){
          if( !pc$h[[hh]]$nev[ii] ){
            I_0 = I_0 + t(Jac_0) %*% IOmega %*% Jac_0
          }
        }

        # information matrix for emplacement condition dependent
        # forward model coefficients
        if( ptbeta > 0 && !is.null(I_t[[tt]]) ){
          if( !pc$h[[hh]]$nev[ii] ){
            I_t[[tt]] = I_t[[tt]] + t(Jac_t) %*% IOmega %*% Jac_t
          }
        }

        # "cross" information matrix between common and emplacement
        # condition dependent forward model coefficients
        if( pbeta > 0 && ptbeta > 0 && !is.null(I_0t[[tt]]) ){
          if( !pc$h[[hh]]$nev[ii] ){
            I_0t[[tt]] = I_0t[[tt]] + t(Jac_0) %*% IOmega %*% Jac_t
          }
        }

        # extract details pertaining to new event inference
        if( pc$h[[hh]]$nev[ii] ){
          if( pbeta > 0 ){
            Jac_0_nev = Jac_0
            Jac_h0t = Jac_0_nev
            if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ){
              t_h0 = tt
              Jac_t_nev = Jac_t
              Jac_h0t = cbind(Jac_h0t,Jac_t_nev)
            }
          } else if( ptbeta > 0 && pc$h[[hh]]$ptbeta[tt] > 0 ) {
            t_h0 = tt
            Jac_t_nev = Jac_t
            Jac_h0t = Jac_t_nev
          } else { t_h0 = tt }
          Omega_h0 = Omega
        }

        # compute components of information matrices needed for
        # errors-in-variables correction
        if( "eiv" %in% pnames && !is.null(pc$h[[hh]]$eiv[[ii]]) ){
          I_wc[pc$h[[hh]]$eiv[[ii]],pc$h[[hh]]$eiv[[ii]]] =
            I_wc[pc$h[[hh]]$eiv[[ii]],pc$h[[hh]]$eiv[[ii]]] +
            t(g_w) %*% IOmega %*% g_w
          if( pbeta > 0 ){
            I_0w[,pc$h[[hh]]$eiv[[ii]]] = t(Jac_0) %*% IOmega %*% g_w
            if( ptbeta > 0 && !is.null(I_tw[[tt]]) ){
              I_tw[[tt]][,pc$h[[hh]]$eiv[[ii]]] =
                                          t(Jac_t) %*% IOmega %*% g_w
            }
          } else if( ptbeta > 0 && !is.null(I_tw[[tt]]) ){
            I_tw[[tt]][,pc$h[[hh]]$eiv[[ii]]] =
                                          t(Jac_t) %*% IOmega %*% g_w
          }
        }
      }
    }

    # compute information matrix for new event inference parameters
    # with no errors-in-variables
    if( pbeta > 0 ){ 
      S_0 = I_0
      if( ptbeta > 0 ){
        II_t = C_0t = vector("list",Th)
        for( tt in 1:Th ){
          if( !is.null(I_t[[tt]]) ){
            I_t[[tt]] = forceSymmetric(I_t[[tt]])
            Catch = pc$tryCatch.W.E(chol(I_t[[tt]]))
            if( is(Catch$value,"Matrix") ){ C_t = Catch$value
            } else {
              cat(paste("t = ",tt,"; I(beta_t) is not positive ",
                        "definite.\n",sep=""))
              return(Cmat)
            }
            II_t[[tt]] = chol2inv(C_t)
            C_0t[[tt]] = I_0t[[tt]] %*% II_t[[tt]] 
            S_0 = S_0 - C_0t[[tt]] %*% t(I_0t[[tt]])
          }
        }
        S_0 = forceSymmetric(S_0)
        Catch = pc$tryCatch.W.E(chol(S_0))
        if( is(Catch$value,"Matrix") ){ C_0 = Catch$value
        } else {
          cat("Adjusted I(beta_0) is not positive definite.\n")
          return(Cmat)
        }
        IS_0 = chol2inv(C_0)
        if( !is.null(I_t[[t_h0]]) ){
          C_c_0t = -IS_0 %*% C_0t[[t_h0]]
          C_c = cbind(IS_0, C_c_0t)
          C_c = rbind(C_c, cbind(t(C_c_0t), II_t[[t_h0]] -
                                 t(C_0t[[t_h0]]) %*% C_c_0t))
        } else { C_c = IS_0 }
      } else {
        S_0 = forceSymmetric(S_0)
        Catch = pc$tryCatch.W.E(chol(S_0))
        if( is(Catch$value,"Matrix") ){ C_0 = Catch$value
        } else {
          cat("I(beta_0) is not positive definite.\n")
          return(Cmat)
        }
        C_c = chol2inv(C_0);
      }
    } else if( ptbeta > 0 && !is.null(I_t[[t_h0]]) ){
      I_t[[t_h0]] = forceSymmetric(I_t[[t_h0]])
      Catch = pc$tryCatch.W.E(chol(I_t[[t_h0]]))
      if( is(Catch$value,"Matrix") ){ C_t = Catch$value
      } else {
        cat(paste("t = ",t_h0,"; I(beta_t) is not positive ",
                   "definite.\n",sep=""))
        return(Cmat)
      }
      C_c = chol2inv(C_t)
    } else { C_c == NULL }
    iInvOmega = FALSE
    if( !is.null(C_c) ){
      P_c = Jac_h0t %*% C_c %*% t(Jac_h0t)
      Sigma_h0c = Omega_h0 + P_c
    } else {
      iInvOmega = TRUE
      Sigma_h0c = Omega_h0
    }
    Sigma_h0c = forceSymmetric(Sigma_h0c)
    Catch = pc$tryCatch.W.E(chol(Sigma_h0c))
    if( is(Catch$value,"Matrix") ){ CSigma_h0c = Catch$value
    } else {
      cat("Adjusted Omega is not positive definite.\n")
      return(Cmat)
    }
    ISigma_h0c = chol2inv(CSigma_h0c)
    if( iInvOmega ){ IOmega_h0 = ISigma_h0c
    } else {
      rInvOmega = FALSE
      if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
        Catch = pc$tryCatch.W.E(chol(Omega_h0))
        if( is(Catch$value,"Matrix") ){ COmega_h0 = Catch$value
        } else {
          cat("Omega_0 is not positive definite.\n")
          return(Cmat)
        }
        IOmega_h0 = chol2inv(COmega_h0)
        rInvOmega = TRUE
      }
    }
    if( "itheta0" %in% pnames ){
      I_nev[pc$h[[hh]]$itheta0,pc$h[[hh]]$itheta0] =
      I_nev[pc$h[[hh]]$itheta0,pc$h[[hh]]$itheta0] +
      t(Jac_th0) %*% ISigma_h0c %*% Jac_th0
      if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
        I_nev_0[pc$h[[hh]]$itheta0,pc$h[[hh]]$itheta0] =
        I_nev_0[pc$h[[hh]]$itheta0,pc$h[[hh]]$itheta0] +
        t(Jac_th0) %*% IOmega_h0 %*% Jac_th0
      }
    } else {
      I_nev = I_nev + t(Jac_th0) %*% ISigma_h0c %*% Jac_th0
      if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
        I_nev_0 = I_nev_0 + t(Jac_th0) %*% IOmega_h0 %*% Jac_th0
      }
    }

    # calculate quantities needed for correction of information matrix
    # for new event inference parameters due to errors-in-variables
    # yields
    if( "eiv" %in% pnames && !is.null(unlist(pc$h[[hh]]$eiv)) ){
      if( !iInvOmega ){
        if( !rInvOmega ){
          Catch = pc$tryCatch.W.E(chol(Omega_h0))
          if( is(Catch$value,"Matrix") ){ COmega_h0 = Catch$value
          } else {
            cat("Omega_0 is not positive definite.\n")
            return(Cmat)
          }
          IOmega_h0 = chol2inv(COmega_h0)
        }
      }
      if( pbeta > 0 ){
        S_0 = I_0 + t(Jac_0_nev) %*% IOmega_h0 %*% Jac_0_nev
        I_00_nev = t(Jac_0_nev) %*% IOmega_h0 %*% Jac_th0
        if( ptbeta > 0 ){
          II_t = C_0t = vector("list",Th)
          for( tt in 1:Th ){
            if( !is.null(I_t[[tt]]) ){
              S_t = I_t[[tt]]
              S_0t = I_0t[[tt]]
              if( tt == t_h0 ){
                S_t = S_t + t(Jac_t_nev) %*% IOmega_h0 %*% Jac_t_nev
                S_0t = S_0t + t(Jac_0_nev) %*% IOmega_h0 %*% Jac_t_nev
              }
              S_t = forceSymmetric(S_t)
              Catch = pc$tryCatch.W.E(chol(S_t))
              if( is(Catch$value,"Matrix") ){ C_t = Catch$value
              } else {
                cat(paste("t = ",tt,"; Adjusted full I(beta_t) is ",
                          "not positive definite.\n",sep=""))
                return(Cmat)
              }
              II_t[[tt]] = chol2inv(C_t)
              C_0t[[tt]] = S_0t %*% II_t[[tt]]
              S_0 = S_0 - C_0t[[tt]] %*% t(S_0t)
            }
          }
          S_0 = forceSymmetric(S_0)
          Catch = pc$tryCatch.W.E(chol(S_0))
          if( is(Catch$value,"Matrix") ){ C_0 = Catch$value
          } else {
            cat("Adjusted full I(beta_0) is not positive definite.\n")
            return(Cmat)
          }
          IS_0 = chol2inv(C_0)
          C_cw = vector("list",Th+1)
          C_cw[[1]] = I_0w
          for( tt in 1:Th ){
            if( !is.null(I_t[[tt]]) ){
              C_cw[[1]] = C_cw[[1]] - C_0t[[tt]] %*% I_tw[[tt]]
            }
          }
          C_cw[[1]] = IS_0 %*% C_cw[[1]]
          for( tt in 1:Th ){
            if( !is.null(I_t[[tt]]) ){
              C_cw[[tt+1]] = II_t[[tt]] %*% I_tw[[tt]] -
                             t(C_0t[[tt]]) %*% C_cw[[1]]
            }
          }
          Q_hcw = t(I_0w) %*% C_cw[[1]]
          for( tt in 1:Th ){
            if( !is.null(I_t[[tt]]) ){
              Q_hcw = Q_hcw + t(I_tw[[tt]]) %*% C_cw[[tt+1]]
            }
          }
          B_hcw = t(C_cw[[1]]) %*% I_00_nev
          if( !is.null(I_t[[t_h0]]) ){
            I_t0_nev = t(Jac_t_nev) %*% IOmega_h0 %*% Jac_th0
            B_hcw = B_hcw + t(C_cw[[t_h0+1]]) %*% I_t0_nev
          }
        } else {
          S_0 = forceSymmetric(S_0)
          Catch = pc$tryCatch.W.E(chol(S_0))
          if( is(Catch$value,"Matrix") ){ C_0 = Catch$value
          } else {
            cat("Full I(beta_0) is not positive definite.\n")
            return(Cmat)
          }
          IS_0 = chol2inv(C_0)
          Q_hcw = t(I_0w) %*% IS_0 %*% I_0w
          B_hcw = t(I_0w) %*% IS_0 %*% I_00_nev
        }
      } else if( ptbeta > 0 ){
        C_cw = vector("list",Th)
        Q_hcw = Matrix(0,pc$nsource,pc$nsource,sparse=FALSE)
        for( tt in 1:Th ){
          if( !is.null(I_t[[tt]]) ){
            S_t = I_t[[tt]]
            if( tt == t_h0 ){
              S_t = S_t + t(Jac_t_nev) %*% IOmega_h0 %*% Jac_t_nev
            }
            S_t = forceSymmetric(S_t)
            Catch = pc$tryCatch.W.E(chol(S_t))
            if( is(Catch$value,"Matrix") ){ C_t = Catch$value
            } else {
              cat(paste("t = ",tt,"; Full I(beta_t) is not ",
                        "positive definite.\n",sep=""))
              return(Cmat)
            }
            II_t = chol2inv(C_t)
            C_cw[[tt]] = II_t %*% I_tw[[tt]]
            Q_hcw = Q_hcw + t(I_tw[[tt]]) %*% C_cw[[tt]]
          }
        }
        B_hcw = 0
        if( !is.null(I_t[[t_h0]]) ){
          I_t0_nev = t(Jac_t_nev) %*% IOmega_h0 %*% Jac_th0
          B_hcw = B_hcw + t(C_cw[[t_h0]]) %*% I_t0_nev
        }
      }
      Q_cw = Q_cw + Q_hcw 
      if( "itheta0" %in% pnames ){
        B_cw[,pc$h[[hh]]$itheta0] = B_cw[,pc$h[[hh]]$itheta0] +
                                    B_hcw
      } else {
        B_cw = B_cw + B_hcw
      }
    }
  }

  # calculate corrected information matrix for new event inference
  # parameters due to errors-in-variables yields
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){ 
    I_eiv = I_wc + I_ww - Q_cw
    I_eiv = forceSymmetric(I_eiv)
    Catch = pc$tryCatch.W.E(chol(I_eiv))
    if( is(Catch$value,"Matrix") ){ CI_eiv = Catch$value
    } else {
      cat("Adjusted I(w) is not positive definite.\n")
      return(Cmat)
    }
    II_eiv = chol2inv(CI_eiv)
    I_nev = I_nev - t(B_cw) %*% II_eiv %*% B_cw
  } 

  # calculate asymptotic covariance matrix for new event inference
  # parameters
  I_nev = forceSymmetric(I_nev)
  Catch = pc$tryCatch.W.E(chol(I_nev))
  if( is(Catch$value,"Matrix") ){ CI_nev = Catch$value
  } else {
    cat("Adjusted I(theta_0) is not positive definite.\n")
    return(Cmat)
  }
  II_nev = chol2inv(CI_nev)
  Cmat$II_nev = II_nev
  if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
    I_nev_0 = forceSymmetric(I_nev_0)
    Catch = pc$tryCatch.W.E(chol(I_nev_0))
    if( is(Catch$value,"Matrix") ){ CI_nev_0 = Catch$value
    } else {
      cat("Adjusted I_0(theta_0) is not positive definite.\n")
      return(Cmat)
    } 
    II_nev_0 = chol2inv(CI_nev_0)
    Cmat$II_nev_0 = II_nev_0
  }

  theta0Hat = opt$theta0
  I_nev_it = I_nev
  if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
    I_nev_0_it = I_nev_0
  }
  iit = FALSE
  if( exists("itransform",where=pc,inherits=FALSE) ){
    if( pc$itransform ){
      iit = TRUE
      Tf = pc$j_tau(theta0Hat, pc=pc)
      theta0Hat = pc$tau(theta0Hat, pc=pc)
    }
  }
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    iit = TRUE
    dtheta0Hat = rep(1,pc$ntheta0)
    if( length(pc$itheta0_bounds[[1]]) > 0 ){
      dtheta0Hat[pc$itheta0_bounds[[1]]] =
        pc$dnotExp(theta0Hat[pc$itheta0_bounds[[1]]])
    }
    if( length(pc$itheta0_bounds[[2]]) > 0 ){
      dtheta0Hat[pc$itheta0_bounds[[2]]] =
        -pc$dnotExp(theta0Hat[pc$itheta0_bounds[[2]]])
    }
    if( length(pc$itheta0_bounds[[3]]) > 0 ){
      tau = pc$notExp(theta0Hat[pc$itheta0_bounds[[3]]])
      dtheta0Hat[pc$itheta0_bounds[[3]]] =
        pc$dnotExp(theta0Hat[pc$itheta0_bounds[[3]]])*
        pc$theta0_range/(1+tau)^2
    }
    Tf_b = diag(pc$ntheta0)
    diag(Tf_b) = dtheta0Hat
    I_nev_it = t(Tf_b) %*% I_nev_it %*% Tf_b
    if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
      I_nev_0_it = t(Tf_b) %*% I_nev_0_it %*% Tf_b
    }
  }
  if( exists("itransform",where=pc,inherits=FALSE) ){
    if( pc$itransform ){
      I_nev_it = t(Tf) %*% I_nev_it %*% Tf
      if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
        I_nev_0_it = t(Tf) %*% I_nev_0_it %*% Tf
      }
    }
  }
  if( iit ){
    I_nev_it = forceSymmetric(I_nev_it)
    Catch = pc$tryCatch.W.E(chol(I_nev_it))
    if( is(Catch$value,"Matrix") ){ CI_nev_it = Catch$value
    } else {
      cat(paste("Adjusted I(theta_0) (untransformed theta_0) is ",
                "not positive definite.\n",sep=""))
      return(Cmat)
    }
    II_nev_it = chol2inv(CI_nev_it)
    Cmat$II_nev_it = II_nev_it
    if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
      I_nev_0_it = forceSymmetric(I_nev_0_it)
      Catch = pc$tryCatch.W.E(chol(I_nev_0_it))
      if( is(Catch$value,"Matrix") ){ CI_nev_0_it = Catch$value
      } else {
        cat(paste("Adjusted I_0(theta_0) (untransformed theta_0) is ",
                  "not positive definite.\n",sep=""))
        return(Cmat)
      }
      II_nev_0_it = chol2inv(CI_nev_0_it)
      Cmat$II_nev_0_it = II_nev_0_it
    }
  }
  Cmat$acov = 1

  return(Cmat)
}
