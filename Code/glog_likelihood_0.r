########################################################################
#                                                                      #
# This file contains code for calculating the gradient of the log-     #
# likelihood of the new event parameters. The forward model            #
# coefficients and error model variance components, and (optionally)   #
# errors-in-variables yields, are fixed at values determined by the    #
# calibration data.                                                    #
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

gll_0 = function(x, pc)
{
  # transform inference parameters if necessary
  if( exists("itransform",where=pc,inherits=FALSE) ){
    if( pc$itransform ){
      dx = pc$j_tau(x, pc=pc)
      x = pc$tau(x, pc=pc)
    }
  }
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    ith0_bds = pc$itheta0_bounds
    dx_b = rep(1,pc$ntheta0)
    if( length(ith0_bds[[1]]) > 0 ){
      dx_b[ith0_bds[[1]]] = pc$dnotExp(x[ith0_bds[[1]]])
    }
    if( length(ith0_bds[[2]]) > 0 ){
      dx_b[ith0_bds[[2]]] = -pc$dnotExp(x[ith0_bds[[2]]])
    }
    if( length(ith0_bds[[3]]) > 0 ){
      tau = pc$notExp(x[ith0_bds[[3]]])
      dx_b[ith0_bds[[3]]] = pc$dnotExp(x[ith0_bds[[3]]])*
                            pc$theta0_range/(1+tau)^2
    }
    x = pc$transform(x, pc=pc)
  }

  # log-likelihood gradient due to new event data

  # initialize gradient vectors
  # new event parameters
  gr_th0 = numeric(pc$ntheta0)

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
      beta = pc$h[[hh]]$beta
    } else { pbeta = 0 }
    if( pc$ptbeta > 0 && "ptbeta" %in% pnames ){
      ptbeta = sum(pc$h[[hh]]$ptbeta)
      tbeta = pc$h[[hh]]$tbeta
    } else { ptbeta = 0 }

    # number of responses for phenomenology "hh"
    Rh = pc$h[[hh]]$Rh

    # number of responses for source "0"
    n_h0 = pc$h[[hh]]$n0
    n_h0_tot = sum(n_h0)

    # setup argument for forward model/jacobian call
    if( "itheta0" %in% pnames ){
      Arg = "(x[pc$h[[hh]]$itheta0],pm)"
    } else { Arg = "(x,pm)" }

    # named parameters in forward model/jacobian call
    if( "theta0_names" %in% pnames ){
      pm$theta_names = pc$h[[hh]]$theta0_names
    }

    # extract forward model parameters for emplacement condition "tt"
    # associated with source "0"
    if( ptbeta > 0 ){
      for( rr in 1:Rh ){
        if( n_h0[rr] > 0 ){
          tt = as.numeric(as.character(pc$h[[hh]]$X0[[rr]]$Type[1]))
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

    # Jacobian matrix for new event inference parameters
    Jac_th0 = NULL

    # iterate over responses "rr"
    for(rr in 1:Rh){
      if( n_h0[rr] > 0 ){
        # covariate matrix for source "0"
        pm$X = pc$h[[hh]]$X0[[rr]]

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
        pm$beta = beta_t
        if( "iResponse" %in% pnames ){
          pm$iresp = pc$h[[hh]]$iResponse[rr]
        }

        # calculate forward model
        fcall = paste("pc$ffm$",pc$h[[hh]]$f0[rr],Arg,sep="")
        yhat = eval(parse(text=fcall))
        if( any(is.nan(yhat)) ){ return(NaN) }

        # calculate residual vector
        resid = c(resid,pc$h[[hh]]$Y0[[rr]] - yhat)

        # calculate components of Jacobian matrix
        gcall = paste("pc$gfm$",pc$h[[hh]]$g0[rr],Arg,sep="")
        jac = eval(parse(text=gcall)) 
        Jac_th0 = rbind(Jac_th0,jac)
      }
    }

    # calculate components of log-likelihood gradient
    # new event inference parameters
    g_th0 = t(Jac_th0) %*% pc$h[[hh]]$IOmega %*% resid
    if( "itheta0" %in% pnames ){
      if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
        g_th0 = g_th0 * dx_b[pc$h[[hh]]$itheta0]
      }
      if( exists("itransform",where=pc,inherits=FALSE) ){
        if( pc$itransform ){
          g_th0 =
            t(dx[pc$h[[hh]]$itheta0,pc$h[[hh]]$itheta0]) %*%
            g_th0
        }
      }
      gr_th0[pc$h[[hh]]$itheta0] = gr_th0[pc$h[[hh]]$itheta0] +
                                   g_th0
    } else {
      if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
        g_th0 = g_th0 * dx_b
      }
      if( exists("itransform",where=pc,inherits=FALSE) ){
        if( pc$itransform ){ g_th0 = t(dx) %*% g_th0 }
      }
      gr_th0 = gr_th0 + g_th0
    }
  }
  return(as.vector(gr_th0))
}
