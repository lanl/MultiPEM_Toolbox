########################################################################
#                                                                      #
# This file contains code for the empirical/physics source and source- #
# to-sensor path attenuation forward models for each observed response #
# within each phenomenology. Only new event inference parameters are   #
# input to these functions (forward model parameters are assumed       #
# fixed).                                                              #
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

f0_s <- function(zeta, params)
{
  beta <- params$beta
  beta[3] <- -params$notExp(beta[3])
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  theta <- matrix(rep(zeta,each=nev),ncol=length(zeta))
  colnames(theta) <- params$theta_names
  X <- cbind(theta, params$X)
  
  y_s <- beta[1]
  w_s <- X[,"W"]*params$yield_scaling
  if ("C2N" %in% colnames(X)){
    c_s <- X[,"C2N"]*params$yield_scaling
  } else {
    c_s <- 0
  }
  t1 <- w_s+c_s
  y_s <- y_s+beta[2]*(X[,"lRange"]-t1)
  h_s <- X[,"HOB"]/exp(t1)
  y_s <- y_s+beta[3]/(1+exp(-beta[4]*h_s-beta[5]))
  if (params$iresp) { y_s <- y_s+t1 }
  names(y_s) <- NULL

  return(y_s)
}

#
# Acoustic forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information 
#           $beta : empirical parameters
#           $theta_names : names of unknown covariates
#           $iresp : TRUE for impulse, FALSE for duration
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $pressure_scaling : pressure scaling factor (nominally 1/3)
#           $temp_scaling : temperature scaling factor (nominally 1/2)
#           $X : matrix of known covariates
#

f0_a <- function(zeta, params)
{
  beta <- params$beta
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  theta <- matrix(rep(zeta,each=nev),ncol=length(zeta))
  colnames(theta) <- params$theta_names
  X <- cbind(theta, params$X)

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
  h_s <- X[,"HOB"]*exp(-t2)
  t3 <- beta[3]*h_s
  y_s <- y_s+t3
  it3 <- which(t3 <= 0)
  y_s[it3] <- y_s[it3]-log(1+exp(t3[it3]))
  it3 <- which(t3 > 0)
  y_s[it3] <- y_s[it3]-t3[it3]-log(1+exp(-t3[it3]))
  names(y_s) <- NULL

  return(y_s)
}

#
# Optical forward model
#  zeta   : vector containing all unknown parameters
#  params : list containing all auxiliary information
#           $beta : empirical parameters
#           $theta_names : names of unknown covariates
#           $yield_scaling : yield scaling factor (nominally 1/3)
#           $pressure_scaling : pressure scaling factor (nominally 1/3)
#           $temp_scaling : temperature scaling factor (nominally 1/2)
#           $X : matrix of known covariates
#

f0_o <- function(zeta, params)
{
  beta <- params$beta
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  theta <- matrix(rep(zeta,each=nev),ncol=length(zeta))
  colnames(theta) <- params$theta_names
  X <- cbind(theta, params$X)

  w_s <- X[,"W"]*params$yield_scaling
  t_s <- X[,"logTempSc"]*params$temp_scaling
  p_s <- X[,"logPressureSc"]*params$pressure_scaling
  y_s <- beta[1]+w_s+p_s-t_s
  h_s <- X[,"HOB"]*exp(p_s-w_s)
  y_s <- y_s+beta[2]*exp(-abs(h_s))
  t1 <- 100*(h_s+0.1)
  y_s <- y_s+t1
  it1 <- which(t1 <= 0)
  y_s[it1] <- y_s[it1]-log(1+exp(t1[it1]))
  it1 <- which(t1 > 0)
  y_s[it1] <- y_s[it1]-t1[it1]-log(1+exp(-t1[it1]))
  names(y_s) <- NULL

  return(y_s)
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

f0_c <- function(zeta, params)
{
  beta <- params$beta
  if (is.null(params$X)) { nev <- 1
  } else { nev <- nrow(params$X) }
  theta <- matrix(rep(zeta,each=nev),ncol=length(zeta))
  colnames(theta) <- params$theta_names
  X <- cbind(theta, params$X)

  y_s <- beta[1]
  if ("yield_scaling" %in% names(params)){
    w_s <- X[,"W"]*params$yield_scaling
  } else {
    w_s <- X[,"W"]*beta[2]
  }
  y_s <- y_s + w_s
  names(y_s) <- NULL

  return(y_s)
}

notExp <- function(x)
{
  f <- x
  ind <- x > 1
  f[ind] <- exp(1) * (x[ind]^2 + 1)/2
  ind <- (x <= 1) & (x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  x[ind] <- -x[ind]
  f[ind] <- exp(1) * (x[ind]^2 + 1)/2
  f[ind] <- 1/f[ind]
  f
}
