########################################################################
#                                                                      #
# This file contains code for the empirical/physics source and source- #
# to-sensor path attenuation forward models for each observed response #
# within each phenomenology.                                           #
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

#
# Seismic forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information
#           $pbeta : number of empirical parameters
#           $cal : indicator of calibration inference parameters
#           $cal_par_names : names of calibration inference parameters
#           $ncalp : number of calibration inference parameters
#           $theta_names : names of unknown covariates
#           $iresp : TRUE for displacement, FALSE for velocity
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $X : matrix of known covariates
#

f_s <- function(zeta, params)
{
  beta <- zeta[1:params$pbeta]
  zeta <- zeta[-(1:params$pbeta)]
  X <- params$X
  if (is.null(X)) { nev <- 1
  } else { nev <- nrow(X) }
  if (params$cal) {
    calp <- matrix(rep(zeta[1:params$ncalp],each=nev),
                       ncol=params$ncalp)
    colnames(calp) <- params$cal_par_names
    ic2n_cp <- which(colnames(calp) %in% "C2N")
    if ("C2N" %in% colnames(X) && length(ic2n_cp) > 0){
      ic2n_x <- which(colnames(X) %in% "C2N")
      if( ncol(X) > 1 ){ X <- X[,-ic2n_x,drop=FALSE]
      } else { X <- NULL }
      if( ncol(calp) > 1 ){
        X <- cbind(calp[,-ic2n_cp,drop=FALSE], X)
      }
    } else { X <- cbind(calp, X) }
    zeta <- zeta[-(1:params$ncalp)]
  }
  ntheta = length(zeta)
  if (ntheta > 0) {
    theta <- matrix(rep(zeta,each=nev),ncol=ntheta)
    colnames(theta) <- params$theta_names
    X <- cbind(theta, X)
  }
  
  y_s <- beta[1]
  w_s <- X[,"W"]*params$yield_scaling
  if ("C2N" %in% colnames(X)){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  t1 <- w_s+c_s
  y_s <- y_s+beta[2]*(X[,"lRange"]-t1)
  if (params$iresp) { y_s <- y_s+t1 }
  names(y_s) <- NULL

  return(y_s)
}

#
# Acoustic forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information 
#           $pbeta : number of empirical parameters
#           $cal : indicator of calibration inference parameters
#           $cal_par_names : names of calibration inference parameters
#           $ncalp : number of calibration inference parameters
#           $theta_names : names of unknown covariates
#           $iresp : TRUE for impulse, FALSE for duration
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $pressure_scaling : pressure scaling factor (nominally 1/3)
#           $temp_scaling : temperature scaling factor (nominally 1/2)
#           $X : matrix of known covariates
#

f_a <- function(zeta, params)
{
  beta <- zeta[1:params$pbeta]
  zeta <- zeta[-(1:params$pbeta)]
  X <- params$X
  if (is.null(X)) { nev <- 1
  } else { nev <- nrow(X) }
  if (params$cal) {
    calp <- matrix(rep(zeta[1:params$ncalp],each=nev),
                       ncol=params$ncalp)
    colnames(calp) <- params$cal_par_names
    ic2n_cp <- which(colnames(calp) %in% "C2N")
    if ("C2N" %in% colnames(X) && length(ic2n_cp) > 0){
      ic2n_x <- which(colnames(X) %in% "C2N")
      if( ncol(X) > 1 ){ X <- X[,-ic2n_x,drop=FALSE]
      } else { X <- NULL }
      if( ncol(calp) > 1 ){
        X <- cbind(calp[,-ic2n_cp,drop=FALSE], X)
      }
    } else { X <- cbind(calp, X) }
    zeta <- zeta[-(1:params$ncalp)]
  }
  ntheta = length(zeta)
  if (ntheta > 0) {
    theta <- matrix(rep(zeta,each=nev),ncol=ntheta)
    colnames(theta) <- params$theta_names
    X <- cbind(theta, X)
  }

  w_s <- X[,"W"]*params$yield_scaling
  if ("C2N" %in% colnames(X)){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  t_s <- X[,"logTempSc"]*params$temp_scaling
  p_s <- X[,"logPressureSc"]*params$pressure_scaling
  t1 <- w_s+c_s
  t2 <- t1-p_s
  y_s <- beta[1]+t1+p_s-t_s
  if (params$iresp) { y_s <- y_s+p_s }
  y_s <- y_s+beta[2]*(X[,"lRange"]-t2) 
  names(y_s) <- NULL

  return(y_s)
}
