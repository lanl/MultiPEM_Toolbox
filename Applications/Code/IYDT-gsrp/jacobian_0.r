########################################################################
#                                                                      #
# This file contains code for computing gradients of the               #
# empirical/physics source and source-to-sensor path attenuation       #
# forward models for each observed response within each phenomenology. #
# Only new event inference parameters are input to these functions     #
# (forward model parameters are assumed fixed).                        #
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
#           $beta : empirical parameters
#           $theta_names : names of unknown covariates
#           $iresp : TRUE for displacement, FALSE for velocity
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $X : matrix of known covariates
#

g0_s <- function(zeta, params)
{
  beta <- params$beta
  beta[3] <- -params$notExp(beta[3])
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  ntheta <- length(zeta)
  theta <- matrix(rep(zeta,each=nev),ncol=ntheta)
  colnames(theta) <- params$theta_names
  X <- cbind(theta,params$X)
  
  w_s <- X[,"W"]*params$yield_scaling
  if ("C2N" %in% colnames(X)){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  t1 <- exp(-w_s-c_s)
  h_s <- X[,"HOB"]*t1
  t2 <- beta[4]*h_s+beta[5]
  t3 <- beta[3]/(1+exp(t2))/(1+exp(-t2))
  jtheta_s <- matrix(-(beta[2]+beta[4]*h_s*t3)*params$yield_scaling,
                     ncol=1)
  if (params$iresp) { jtheta_s <- jtheta_s+params$yield_scaling }
  if (ntheta > 1){ jtheta_s <- cbind(jtheta_s,beta[4]*t3*t1) }
  rownames(jtheta_s) <- colnames(jtheta_s) <- NULL
  return(jtheta_s)
}

#
# Acoustic forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information 
#           $beta : empirical parameters
#           $theta_names : names of unknown covariates
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $pressure_scaling : pressure scaling factor (nominally 1/3)
#           $X : matrix of known covariates
#

g0_a <- function(zeta, params)
{
  beta <- params$beta
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  ntheta <- length(zeta)
  theta <- matrix(rep(zeta,each=nev),ncol=ntheta)
  colnames(theta) <- params$theta_names
  X <- cbind(theta,params$X)

  w_s <- X[,"W"]*params$yield_scaling
  if ("C2N" %in% colnames(X)){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  p_s <- X[,"logPressureSc"]*params$pressure_scaling
  t1 <- exp(p_s-w_s-c_s)
  h_s <- X[,"HOB"]*t1
  t2 <- 1/(1+exp(beta[3]*h_s))
  jtheta_s <- matrix((1-beta[2]-beta[3]*h_s*t2)*params$yield_scaling,
                     ncol=1)
  if (ntheta > 1){ jtheta_s <- cbind(jtheta_s,beta[3]*t2*t1) }
  rownames(jtheta_s) <- colnames(jtheta_s) <- NULL
  return(jtheta_s)
}

#
# Optical forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information
#           $beta : empirical parameters
#           $theta_names : names of unknown covariates
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $X : matrix of known covariates
#

g0_o <- function(zeta, params)
{
  beta <- params$beta
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  ntheta <- length(zeta)
  theta <- matrix(rep(zeta,each=nev),ncol=ntheta)
  colnames(theta) <- params$theta_names
  X <- cbind(theta,params$X)

  w_s <- X[,"W"]*params$yield_scaling
  t1 <- exp(-w_s)
  h_s <- X[,"HOB"]*t1
  t2 <- params$notExp(beta[4])
  t3 <- h_s/t2
  t4 <- t3^2
  t5 <- exp(-t4)
  t6 <- params$notExp(beta[3])
  t7 <- 1+t6*t5
  jtheta_s <- as.matrix(beta[2]+2*t6*t4*t5*params$yield_scaling/t7)
  if (ntheta > 1){
    jtheta_s <- cbind(jtheta_s,-2*t6*t1*t3*t5/t2/t7)
  }
  rownames(jtheta_s) <- colnames(jtheta_s) <- NULL
  return(jtheta_s)
}

#
# Crater forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information
#           $beta : empirical parameters
#           $theta_names : names of unknown covariates
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $X : matrix of known covariates
#

g0_c <- function(zeta, params)
{
  beta <- params$beta
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  ntheta <- length(zeta)
  theta <- matrix(rep(zeta,each=nev),ncol=ntheta)
  colnames(theta) <- params$theta_names
  X <- cbind(theta, params$X)

  if ("yield_scaling" %in% names(params)){
    jtheta_s <- matrix(rep(params$yield_scaling,nev),ncol=1)
  } else {
    jtheta_s <- matrix(rep(beta[2],nev),ncol=1)
  }
  rownames(jtheta_s) <- colnames(jtheta_s) <- NULL
  return(jtheta_s)
}
