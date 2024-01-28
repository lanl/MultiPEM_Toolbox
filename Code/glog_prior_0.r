########################################################################
#                                                                      #
# This file contains code for calculating the gradient of the          #
# log-prior density of the new event parameters.                       #
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

glprior_0 = function(x, pc)
{
  # extract inference parameters for new event
  gr_x = numeric(pc$ntheta0)
  if( exists("itransform",where=pc,inherits=FALSE) ){
    if( pc$itransform ){
      gr_x = gr_x + pc$dlog_absdet_j_tau(x, pc=pc)
      if( exists("itheta0_bounds",where=pc,inherits=FALSE) ||
          ("lp_theta0" %in% names(pc)) ){
        dx = pc$j_tau(x, pc=pc)
      }
      x = pc$tau(x, pc=pc)
    }
  }
  if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
    ith0_bds = pc$itheta0_bounds
    g_x = rep(0,pc$ntheta0)
    if( length(ith0_bds[[1]]) > 0 ){
      djt = pc$d2notExp(x[ith0_bds[[1]]])/
            pc$dnotExp(x[ith0_bds[[1]]])
      g_x[ith0_bds[[1]]] = djt
    }
    if( length(ith0_bds[[2]]) > 0 ){
      djt = pc$d2notExp(x[ith0_bds[[2]]])/
            pc$dnotExp(x[ith0_bds[[2]]])
      g_x[ith0_bds[[2]]] = djt
    }
    if( length(ith0_bds[[3]]) > 0 ){
      tau = pc$notExp(x[ith0_bds[[3]]])
      dtau = pc$dnotExp(x[ith0_bds[[3]]])
      djt = pc$d2notExp(x[ith0_bds[[3]]])/
            dtau - 2*dtau/(1+tau)
      g_x[ith0_bds[[3]]] = djt
    }
    if( exists("itransform",where=pc,inherits=FALSE) ){
      if( pc$itransform ){ g_x = t(dx) %*% g_x }
    }
    gr_x = gr_x + g_x
    if( "lp_theta0" %in% names(pc) ){
      dx_b = rep(1,pc$ntheta0)
      if( length(ith0_bds[[1]]) > 0 ){
        dx_b[ith0_bds[[1]]] = pc$dnotExp(x[ith0_bds[[1]]])
      }
      if( length(ith0_bds[[2]]) > 0 ){
        dx_b[ith0_bds[[2]]] = -pc$dnotExp(x[ith0_bds[[2]]])
      }
      if( length(ith0_bds[[3]]) > 0 ){
        dx_b[ith0_bds[[3]]] = pc$dnotExp(x[ith0_bds[[3]]])*
                                   pc$theta0_range/(1+tau)^2
      }
    }
    x = pc$transform(x, pc=pc)
  }
  if( "lp_theta0" %in% names(pc) ){
    Arg = "(x,pc)"
    glp_theta0_call = paste("pc$glp$",pc$lp_theta0$g,Arg,sep="")
    # evaluate gradient for new event inference parameters
    g_x = eval(parse(text=glp_theta0_call))
    if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
      g_x = g_x * dx_b
    }
    if( exists("itransform",where=pc,inherits=FALSE) ){
       if( pc$itransform ){ g_x = t(dx) %*% g_x }
    }
    gr_x = gr_x + g_x
  }
  return(as.vector(gr_x))
}
