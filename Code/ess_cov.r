########################################################################
#                                                                      #
# This file contains code for calculating the effective sample size    #
# of the calibration (and optionally, new event) data set(s), based on #
# computing the magnitude of the estimated model correlation matrix    #
# evaluated at the MLE.                                                #
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

ess_cov = function(opt, pc)
{
  # use R Matrix package
  require(Matrix)

  # extract level 1 variance components
  if( exists("vc_1",where=opt,inherits=FALSE) ){ vcHat_1 = opt$vc_1 }

  # extract level 2 variance components
  if( exists("vc_2",where=opt,inherits=FALSE) ){ vcHat_2 = opt$vc_2 }

  # extract observational error covariance parameters
  epsHat = opt$eps

  # effective sample size due to calibration data and (optionally)
  # new event data
  ess = 0

  # compute covariance matrices by source
  for(hh in 1:pc$H){
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

    # iterate over sources in phenomenology "hh"
    for(ii in 1:pc$h[[hh]]$nsource){
      # number of responses for source "ii"
      n_hi = pc$h[[hh]]$n[[ii]]
      n_hi_tot = sum(n_hi)

      # covariance matrices
      if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
        Xi_hi1 = vector("list",Rh)
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          Xi_hi2 = vector("list",Rh)
        }
      }

      # iterate over responses "rr"
      for(rr in 1:Rh){
        if( n_hi[rr] > 0 ){
          # calculate components of model covariance matrix
          if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
            if( pc$h[[hh]]$pvc_1[rr] > 0 ){
              Xi_hi1[[rr]] = pc$h[[hh]]$Z1[[ii]][[rr]] %*%
                             Sigma_hr1[[rr]] %*%
                             t(pc$h[[hh]]$Z1[[ii]][[rr]])
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
      Isd = Diagonal(x=1/sqrt(diag(Omega)))
      CorrMat = Isd %*% Omega %*% Isd
      cCorrMat = chol(CorrMat)
      ICorrMat = chol2inv(cCorrMat)

      # calculate effective sample size
      ess = ess + sum(ICorrMat)
    }
  }
  return(as.numeric(ess))
}
