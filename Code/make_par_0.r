########################################################################
#                                                                      #
# This file contains code for organizing the forward model parameter   #
# values and constructing all required fixed error model covariance    #
# matrices for use in reduced log-likelihood and log-prior             #
# calculations.                                                        #
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

pc_0 = function(xfin, pc)
{
  # use R Matrix package
  require(Matrix)

  # remove errors-in-variables inference parameters
  if( pc$nsource > 0 ){ xfin = xfin[-(1:pc$nsource)] }

  # collect common coefficients by response
  # for each phenomenology
  if( pc$pbeta > 0 ){
    for( hh in 1:pc$H ){
      pbeta = sum(pc$h[[hh]]$pbeta)
      if( pbeta > 0 ){
        pc$h[[hh]]$beta = xfin[1:pbeta]
        xfin = xfin[-(1:pbeta)]
      }
    }
  }

  # collect emplacement condition dependent coefficients by response
  # for each phenomenology
  if( pc$ptbeta > 0 ){
    for( hh in 1:pc$H ){
      if( "ptbeta" %in% names(pc$h[[hh]]) ){
        ptbeta = sum(pc$h[[hh]]$ptbeta)
      } else { ptbeta = 0 }
      if( ptbeta > 0 ){
        pc$h[[hh]]$tbeta = xfin[1:ptbeta]
        xfin = xfin[-(1:ptbeta)]
      }
    }
  }

  # extract level 1 variance components
  if( pc$pvc_1 > 0 ){
    vc1_all = xfin[1:pc$pvc_1]    
    xfin = xfin[-(1:pc$pvc_1)]
  }

  # extract level 2 variance components
  if( pc$pvc_1 > 0 && pc$pvc_2 > 0 ){
    vc2_all = xfin[1:pc$pvc_2]
    xfin = xfin[-(1:pc$pvc_2)]
  }

  # initialize forward model quantities used in log-likelihood
  for( hh in 1:pc$H ){
    nsource = pc$h[[hh]]$nsource
    Rh = pc$h[[hh]]$Rh
    # number of responses for source "0"
    pc$h[[hh]]$n0 = pc$h[[hh]]$n[[nsource]]
    # names of new event inference parameters
    pc$h[[hh]]$theta0_names = pc$h[[hh]]$theta_names[[nsource]]
    # covariate matrices and responses for new event data
    pc$h[[hh]]$X0 = vector("list",Rh)
    pc$h[[hh]]$Y0 = vector("list",Rh)
    for( rr in 1:Rh ){
      pc$h[[hh]]$X0[[rr]] = pc$h[[hh]]$X[[nsource]][[rr]]
      pc$h[[hh]]$Y0[[rr]] = pc$h[[hh]]$Y[[nsource]][[rr]]
    }
  }

  # initialize error model quantities used in log-likelihood
  for( hh in 1:pc$H ){
    nsource = pc$h[[hh]]$nsource
    Rh = pc$h[[hh]]$Rh
    # number of responses for source "0"
    n_h0 = pc$h[[hh]]$n0
    n_h0_tot = sum(n_h0)
    # variance component covariance matrices
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      Xi_h01 = vector("list",Rh)
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                              pc$h[[hh]]$pvc_2 > 0) ){
        Xi_h02 = vector("list",Rh)
      }
    }
    for( rr in 1:Rh ){
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        pvc_1 = pc$h[[hh]]$pvc_1[rr]
        if( pvc_1 > 0 ){
          # initialize level 1 covariate matrices
          Z1_0 = pc$h[[hh]]$Z1[[nsource]][[rr]]
          # initialize level 1 bias term covariance matrices
          Sigma_hr1 = Diagonal(pvc_1,exp(vc1_all[1:pvc_1]))
          vc1_all = vc1_all[-(1:pvc_1)]
          # construct level 1 variance component covariance matrices
          Xi_h01[[rr]] = Z1_0 %*% Sigma_hr1 %*% t(Z1_0)
        } else { Xi_h01[[rr]] = Diagonal(n_h0[rr],0) }
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          if( pvc_1 > 0 ){
            pvc_2 = pc$h[[hh]]$pvc_2[rr]
            if( pvc_2 > 0 ){
              # initialize level 2 covariate matrices
              Z2_0 = pc$h[[hh]]$Z2[[nsource]][[rr]]
              # initialize level 2 path levels
              nplev_0 = pc$h[[hh]]$nplev[nsource,rr]
              # initialize level 2 bias term covariance matrices
              Sigma_hr2 = Diagonal(pvc_2,exp(vc2_all[1:pvc_2]))
              vc2_all = vc2_all[-(1:pvc_2)]
              # construct level 2 variance component covariance
              # matrices
              Xi_h02[[rr]] = Z2_0 %*%
                             kronecker(Diagonal(nplev_0),Sigma_hr2) %*%
                             t(Z2_0)
            } else { Xi_h02[[rr]] = Diagonal(n_h0[rr],0) }
          } else { Xi_h02[[rr]] = Diagonal(n_h0[rr],0) }
        }
      }
    }
    # construct observational error covariance matrices
    L_h = Diagonal(Rh,exp(xfin[1:Rh]))
    xfin = xfin[-(1:Rh)]
    if( Rh > 1 ){
      L_h[upper.tri(L_h)] = xfin[1:choose(Rh,2)]
      xfin = xfin[-(1:choose(Rh,2))]
    }
    Sigma_h = t(L_h) %*% L_h
    # calculate model covariance matrix
    Omega = Matrix(0,n_h0_tot,n_h0_tot)
    for( r1 in 1:Rh ){
      if( n_h0[r1] > 0 ){
        st_n0r1 = 0
        if( r1 > 1 ){ st_n0r1 = sum(n_h0[1:(r1-1)]) }
        ir = st_n0r1+(1:n_h0[r1])
        for( r2 in r1:Rh ){
          if( n_h0[r2] > 0 ){
            st_n0r2 = 0
            if( r2 > 1 ){ st_n0r2 = sum(n_h0[1:(r2-1)]) }
            ic = st_n0r2+(1:n_h0[r2])
            Sigma_h0 = Matrix(0,n_h0[r1],n_h0[r2])
            Sigma_h0[pc$h[[hh]]$i[[nsource]]$cov_pairs[[r1]][[r2]]] =
              Sigma_h[r1,r2]
            Omega[ir,ic] = Sigma_h0
            if( r2 > r1 ){ Omega[ic,ir] = t(Sigma_h0) }
          }
        }
      }
    }
    if( any(pc$h[[hh]]$pvc_1 > 0) ){
      Xi_h01 = Xi_h01[!sapply(Xi_h01,is.null)]
      Omega = Omega + bdiag(Xi_h01)
    }
    if( any(pc$h[[hh]]$pvc_1 > 0 & pc$h[[hh]]$pvc_2 > 0) ){
      Xi_h02 = Xi_h02[!sapply(Xi_h02,is.null)]
      Omega = Omega + bdiag(Xi_h02)
    }
    if( any(is.infinite(Omega)) ){
      stop("All elements of Omega must be finite.")
    }
    cOmegaCatch = pc$tryCatch.W.E(chol(Omega))
    if( is(cOmegaCatch$value,"Matrix") ){ cOmega = cOmegaCatch$value
    } else {
      stop("Cholesky factor of Omega is not positive definite.")
    }
    if( kappa(cOmega) > 1000 ){
     stop("Cholesky factor of Omega is ill-conditioned.")
    }
    pc$h[[hh]]$logdet_cOmega = sum(log(diag(cOmega)))
    pc$h[[hh]]$IOmega = chol2inv(cOmega)
  }
  return(pc)
}
