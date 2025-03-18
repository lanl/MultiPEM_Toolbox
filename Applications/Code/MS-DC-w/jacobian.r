########################################################################
#                                                                      #
# This file contains code for computing gradients of the               #
# empirical/physics source and source-to-sensor path attenuation       #
# forward models for each observed response within each phenomenology. #
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

g_s <- function(zeta, params)
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

  jbeta_s <- matrix(rep(1,nev),ncol=1)
  w_s <- X[,"W"]*params$yield_scaling
  ic2n <- "C2N" %in% colnames(X)
  if (ic2n){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  t1 <- w_s+c_s
  jbeta_s <- cbind(jbeta_s,X[,"lRange"]-t1)
  rownames(jbeta_s) <- colnames(jbeta_s) <- NULL

  if (params$cal){
    if (ic2n){
      jcalp_s <- matrix(rep(-beta[2]*params$yield_scaling,nev),ncol=1)
      if (params$iresp) { jcalp_s <- jcalp_s+params$yield_scaling }
    } else {
      jcalp_s <- matrix(0,nev,1)
    }
    rownames(jcalp_s) <- colnames(jcalp_s) <- NULL
  }

  if (ntheta > 0){
    jtheta_s <- matrix(rep(-beta[2]*params$yield_scaling,nev),ncol=1)
    if (params$iresp) { jtheta_s <- jtheta_s+params$yield_scaling }
    rownames(jtheta_s) <- colnames(jtheta_s) <- NULL
  }

  if (params$cal){
    if (ntheta > 0){
      return(list(jbeta=jbeta_s,jcalp=jcalp_s,jtheta=jtheta_s))
    } else {
      return(list(jbeta=jbeta_s,jcalp=jcalp_s))
    }
  } else {
    if (ntheta > 0){
      return(list(jbeta=jbeta_s,jtheta=jtheta_s))
    } else { return(jbeta_s) }
  }
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
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $pressure_scaling : pressure scaling factor (nominally 1/3)
#           $X : matrix of known covariates
#

g_a <- function(zeta, params)
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

  jbeta_s <- matrix(rep(1,nev),ncol=1)
  w_s <- X[,"W"]*params$yield_scaling
  ic2n <- "C2N" %in% colnames(X)
  if (ic2n){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  p_s <- X[,"logPressureSc"]*params$pressure_scaling
  t1 <- w_s+c_s-p_s
  jbeta_s <- cbind(jbeta_s,X[,"lRange"]-t1)
  rownames(jbeta_s) <- colnames(jbeta_s) <- NULL

  if (params$cal){
    if (ic2n){
      jcalp_s <- matrix(rep((1-beta[2])*params$yield_scaling,nev),ncol=1)
    } else { 
      jcalp_s <- matrix(0,nev,1)
    }
    rownames(jcalp_s) <- colnames(jcalp_s) <- NULL
  }

  if (ntheta > 0){
    jtheta_s <- matrix(rep((1-beta[2])*params$yield_scaling,nev),ncol=1)
    rownames(jtheta_s) <- colnames(jtheta_s) <- NULL
  }

  if (params$cal){
    if (ntheta > 0){
      return(list(jbeta=jbeta_s,jcalp=jcalp_s,jtheta=jtheta_s))
    } else {
      return(list(jbeta=jbeta_s,jcalp=jcalp_s))
    }
  } else {
    if (ntheta > 0){
      return(list(jbeta=jbeta_s,jtheta=jtheta_s))
    } else { return(jbeta_s) }
  }
}
