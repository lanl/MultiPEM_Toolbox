########################################################################
#                                                                      #
# This file contains code for calculating the log-prior density of the #
# calibration parameters, with an option to include new event          #
# inference parameters and errors-in-variables for calibration data    #
# yields.                                                              #
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

lprior = function(x, pc)
{
  # use R Matrix package
  require(Matrix)

  # log-prior density
  lp = 0

  # extract scale parameters for variance component priors
  pnames0 = names(pc)
  if( !("A" %in% pnames0) ){
    # total number of parameters
    npars = length(x)
    lA = x[(npars-pc$p_A+1):npars]
    A = exp(lA)
    x = x[-((npars-pc$p_A+1):npars)]
    iA = 0
  } else {
    A = pc$A
    lA = log(A)
    iA = 1
  }
  # extract FGSN parameters for errors-in-variables yield
  # priors
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    # number of remaining parameters
    npars = length(x)
    fgsn_all = x[(npars-pc$K-1):npars]
    lam2 = exp(fgsn_all[2])
    x = x[-((npars-pc$K-1):npars)]
  }
  # extract inference parameters for new event
  if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
    theta0 = x[1:pc$ntheta0]
    if( exists("itransform",where=pc,inherits=FALSE) ){
      if( pc$itransform ){
        lp = lp + pc$log_absdet_j_tau(theta0, pc=pc)
        theta0 = pc$tau(theta0, pc=pc)
      }
    }
    if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
      ith0_bds = pc$itheta0_bounds
      if( length(ith0_bds[[1]]) > 0 ){
        jt = sum(log(pc$dnotExp(theta0[ith0_bds[[1]]])))
        lp = lp + jt
      }
      if( length(ith0_bds[[2]]) > 0 ){
        jt = sum(log(pc$dnotExp(theta0[ith0_bds[[2]]])))
        lp = lp + jt
      }
      if( length(ith0_bds[[3]]) > 0 ){
        tau = pc$notExp(theta0[ith0_bds[[3]]])
        jt = pc$sum_theta0_logrange +
             sum(log(pc$dnotExp(theta0[ith0_bds[[3]]]))) -
             2*sum(log(1+tau))
        lp = lp + jt
      }
      theta0 = pc$transform(theta0, pc=pc)
    }
    if( "lp_theta0" %in% names(pc) ){
      Arg = "(theta0,pc)"
      lp_theta0_call = paste("pc$flp$",pc$lp_theta0$f,Arg,sep="")
      # evaluate log-prior for new event inference parameters
      lp = lp + eval(parse(text=lp_theta0_call))
    }
    x = x[-(1:pc$ntheta0)]
  }
  if( pc$ncalp > 0 ){
    calp = x[1:pc$ncalp]
    if( "lp_calp" %in% names(pc) ){
      Arg = "(calp,pc)"
      lp_calp_call = paste("pc$flp$",pc$lp_calp$f,Arg,sep="")
      # evaluate log-prior for calibration inference parameters
      lp = lp + eval(parse(text=lp_calp_call))
    }
    x = x[-(1:pc$ncalp)]
  }
  # extract errors-in-variables yield parameters
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    w_eiv = x[1:pc$nsource]
    # FGSN prior
    lp = lp + pc$dfgsn(w_eiv,fgsn_all[1],fgsn_all[2],
                       sqrt(lam2),fgsn_all[-(1:2)])
    # alpha
    lp = lp + 0
    # lambda-squared
    lp = lp + pc$dinvgamma(fgsn_all[2],lam2,0.5,0.5)
    # jacobian
    lp = lp + fgsn_all[2]
    # omega
    lp = lp + sum(dnorm(fgsn_all[-(1:2)],0,10,log=TRUE))
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
  # extract source variance components
  if( pc$pvc_1 > 0 ){
    vc1_all = x[1:pc$pvc_1]
    x = x[-(1:pc$pvc_1)]
  }
  # extract path variance components
  if( pc$pvc_2 > 0 ){
    vc2_all = x[1:pc$pvc_2]
    x = x[-(1:pc$pvc_2)]
  }

  for(hh in 1:pc$H){
    # named objects for phenomenology "hh"
    pnames = names(pc$h[[hh]])

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

    # extract forward model parameters for emplacement condition "tt"
    if( ptbeta > 0 ){
      Th = pc$h[[hh]]$Th
      betat = vector("list",Th)
      for( tt in 1:Th ){
        st_betat = 0
        if( tt > 1 ){ st_betat = sum(pc$h[[hh]]$ptbeta[1:(tt-1)]) }
        if( pc$h[[hh]]$ptbeta[tt] > 0 ){
          betat[[tt]] = tbeta[st_betat+(1:pc$h[[hh]]$ptbeta[tt])]
        }
      }
    }

    # extract variance component parameters for
    # phenomenology "hh"
    if( pc$pvc_1 > 0 && any(pc$h[[hh]]$pvc_1 > 0) ){
      # source
      pvc_1 = sum(pc$h[[hh]]$pvc_1)
      lvc_1 = vc1_all[1:pvc_1]
      vc_1 = exp(lvc_1)
      vc1_all = vc1_all[-(1:pvc_1)]
    }
    if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_2 > 0) ){
      # path
      pvc_2 = sum(pc$h[[hh]]$pvc_2)
      lvc_2 = vc2_all[1:pvc_2]
      vc_2 = exp(lvc_2)
      vc2_all = vc2_all[-(1:pvc_2)]
    }

    # number of responses for phenomenology "hh"
    Rh = pc$h[[hh]]$Rh

    # construct observational error
    # covariance matrices
    ell_h = x[1:Rh]
    L_h = Diagonal(Rh,exp(ell_h))
    x = x[-(1:Rh)]
    if( Rh > 1 ){
      L_h[upper.tri(L_h)] = x[1:choose(Rh,2)]
      x = x[-(1:choose(Rh,2))]
    }
    Sigma_h = t(L_h) %*% L_h
    voe_h = diag(Sigma_h)
    lvoe_h = log(voe_h)
    slvoe_h = sum(lvoe_h)
    if( is.infinite(slvoe_h) || is.nan(slvoe_h) ){ return(-Inf) }
    sdoe_h = sqrt(voe_h)
    ISd_h = Diagonal(Rh,1/sdoe_h)
    C_h = ISd_h %*% Sigma_h %*% ISd_h
    ldetC_h = 2*sum(ell_h)-slvoe_h

    # iterate over responses "rr"
    for(rr in 1:Rh){
      # extract forward model parameters for response "rr"
      if( pbeta > 0 && pc$h[[hh]]$pbeta[rr] > 0 ){
        if( "lp_beta" %in% pnames ){
          st_beta = 0
          if( rr > 1 ){ 
            st_beta = sum(pc$h[[hh]]$pbeta[1:(rr-1)])
          }
          betar = beta[st_beta+(1:pc$h[[hh]]$pbeta[rr])]
          Arg = "(betar,pc)"
          lp_beta_call = paste("pc$flp$",pc$h[[hh]]$lp_beta$f[rr],
                               Arg,sep="")
          # evaluate log-prior for common model parameters
          lp = lp + eval(parse(text=lp_beta_call))
        }
      }
      if( ptbeta > 0 ){
        if( "lp_betat" %in% pnames ){
          for( tt in 1:Th ){
            if( pc$h[[hh]]$pbetat[[tt]][rr] > 0 &&
                !is.null(pc$h[[hh]]$lp_betat[[tt]]) ){
              st_betatr = 0
              if( rr > 1 ){
                st_betatr = sum(pc$h[[hh]]$pbetat[[tt]][1:(rr-1)])
              }
              betatr = betat[[tt]][st_betatr+
                                   (1:pc$h[[hh]]$pbetat[[tt]][rr])]
              Arg = "(betatr,pc)"
              lp_betat_call = paste("pc$flp$",
                                    pc$h[[hh]]$lp_betat[[tt]]$f[rr],
                                    Arg,sep="")
              # evaluate log-prior for emplacement condition
              # dependent parameters
              lp = lp + eval(parse(text=lp_betat_call))
            }
          }
        }
      }

      # evaluate log-prior for variance components
      if( pc$pvc_1 > 0 && pc$h[[hh]]$pvc_1[rr] > 0 ){
        st_vc1 = 0
        if( rr > 1 ){
          st_vc1 = sum(pc$h[[hh]]$pvc_1[1:(rr-1)])
        }
        lvc_r = lvc_1[st_vc1+(1:pc$h[[hh]]$pvc_1[rr])]
        vc_r = vc_1[st_vc1+(1:pc$h[[hh]]$pvc_1[rr])]
        if( !("A" %in% pnames0) ){ iA = iA+1 }
        # source 
        lp = lp + pc$lphc(lvc_r,vc_r,lA[iA],A[iA])
        # jacobian
        lp = lp + sum(lvc_r)
      }
      if( pc$pvc_2 > 0 && pc$h[[hh]]$pvc_2[rr] > 0 ){
        st_vc2 = 0 
        if( rr > 1 ){
          st_vc2 = sum(pc$h[[hh]]$pvc_2[1:(rr-1)])
        }
        lvc_r = lvc_2[st_vc2+(1:pc$h[[hh]]$pvc_2[rr])]
        vc_r = vc_2[st_vc2+(1:pc$h[[hh]]$pvc_2[rr])]
        if( !("A" %in% pnames0) ){
          if( pc$pvc_1 == 0 || pc$h[[hh]]$pvc_1[rr] == 0 ){
            iA = iA+1
          }
        }
        # path 
        lp = lp + pc$lphc(lvc_r,vc_r,lA[iA],A[iA])
        # jacobian
        lp = lp + sum(lvc_r)
      }
      if( (pc$pvc_1 > 0 && pc$h[[hh]]$pvc_1[rr] > 0) ||
          (pc$pvc_2 > 0 && pc$h[[hh]]$pvc_2[rr] > 0) ){
        if( !("A" %in% pnames0) ){
          # log prior for scale parameter A
          lp = lp + 0
          # jacobian
          lp = lp + lA[iA]
        }
      }
    }

    # evaluate log-prior for observational error parameters
    # observational error variances
    lp = lp - slvoe_h
    # observational error correlations
    # Lewandowski-Kurowicka-Joe (LKJ) prior
    eta = pc$lp_corr$eta
    lp = lp + (eta - 1)*ldetC_h
    lp = lp + sum((2*eta-2+Rh-1:(Rh-1))*(Rh-1:(Rh-1)))*log(2)+
         sum((Rh-1:(Rh-1))*lbeta(eta+(Rh-1:(Rh-1)-1)/2,
             eta+(Rh-1:(Rh-1)-1)/2))
    # jacobian
    ncpar = Rh*(Rh+1)/2
    Jac_c = Matrix(0,ncpar,ncpar,sparse=FALSE,doDiag=FALSE)
    qq = 0
    k_ij = pc$k_ij
    for( r2 in 1:Rh ){
      for( r1 in 1:r2 ){
        qq = qq+1
        if( r1 == r2 ){
          Jac_c[qq,k_ij(1,r2):k_ij(r2,r2)] = 2*L_h[1:r2,r2]
          Jac_c[qq,k_ij(r2,r2)] = Jac_c[qq,k_ij(r2,r2)]*L_h[r2,r2]
        } else {
          Jac_c[qq,k_ij(1,r1):k_ij(r1,r1)] =
          L_h[1:r1,r2]/sdoe_h[r1]/sdoe_h[r2]-
          C_h[r1,r2]*L_h[1:r1,r1]/voe_h[r1]
          Jac_c[qq,k_ij(r1,r1)] = Jac_c[qq,k_ij(r1,r1)]*L_h[r1,r1]
          Jac_c[qq,k_ij(1,r2):k_ij(r2,r2)] =
          c(L_h[1:r1,r1]/sdoe_h[r1]/sdoe_h[r2]-
            C_h[r1,r2]*L_h[1:r1,r2]/voe_h[r2],
            -C_h[r1,r2]*L_h[(r1+1):r2,r2]/voe_h[r2])
          Jac_c[qq,k_ij(r2,r2)] = Jac_c[qq,k_ij(r2,r2)]*L_h[r2,r2]
        }
      }
    }
    lp = lp + determinant(Jac_c,logarithm=TRUE)$modulus
  }
  return(as.numeric(lp))
}

dinvgamma = function(lv,v,shape,rate)
{
  sum(shape*log(rate)-lgamma(shape)-(shape+1)*lv-rate/v)
}

dfgsn = function(w,alpha,llambda2,lambda,omega)
{
  M = length(omega)
  t1 = (w - alpha)/lambda
  t2 = NULL
  for (ii in 1:M){ t2 = cbind(t2, t1^(2*ii-1)) }
  t3 = t2 %*% omega
  return(sum(dnorm(t1, log=TRUE) + pnorm(t3, log=TRUE) +
         log(2) - llambda2/2))
}

lphc <- function(lv,v,lA,A)
{
  sum(lA - log(A^2 + v) - lv/2 - log(pi))
}

k_ij = function(i,j)
{
  i+choose(j,2)
}
