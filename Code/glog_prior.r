########################################################################
#                                                                      #
# This file contains code for calculating the gradient of the          #
# log-prior density of the calibration parameters, with an option to   #
# include new event inference parameters and errors-in-variables for   #
# calibration data yields.                                             #
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

glprior = function(x, pc)
{
  # use R Matrix package
  require(Matrix)

  # extract scale parameters for variance component priors
  pnames0 = names(pc)
  if( !("A" %in% pnames0) ){
    # total number of parameters
    npars = length(x)
    A = exp(x[(npars-pc$p_A+1):npars])
    x = x[-((npars-pc$p_A+1):npars)]
    iA = 0
  } else {
    A = pc$A
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
    gr_th0 = numeric(pc$ntheta0)
    if( exists("itransform",where=pc,inherits=FALSE) ){
      if( pc$itransform ){
        gr_th0 = gr_th0 + pc$dlog_absdet_j_tau(theta0, pc=pc)
        if( exists("itheta0_bounds",where=pc,inherits=FALSE) ||
            ("lp_theta0" %in% names(pc)) ){
          dtheta0 = pc$j_tau(theta0, pc=pc)
        }
        theta0 = pc$tau(theta0, pc=pc)
      }
    }
    if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
      ith0_bds = pc$itheta0_bounds
      g_th0 = rep(0,pc$ntheta0)
      if( length(ith0_bds[[1]]) > 0 ){
        djt = pc$d2notExp(theta0[ith0_bds[[1]]])/
              pc$dnotExp(theta0[ith0_bds[[1]]])
        g_th0[ith0_bds[[1]]] = djt
      }
      if( length(ith0_bds[[2]]) > 0 ){
        djt = pc$d2notExp(theta0[ith0_bds[[2]]])/
              pc$dnotExp(theta0[ith0_bds[[2]]])
        g_th0[ith0_bds[[2]]] = djt
      }
      if( length(ith0_bds[[3]]) > 0 ){
        tau = pc$notExp(theta0[ith0_bds[[3]]])
        dtau = pc$dnotExp(theta0[ith0_bds[[3]]])
        djt = pc$d2notExp(theta0[ith0_bds[[3]]])/
              dtau - 2*dtau/(1+tau)
        g_th0[ith0_bds[[3]]] = djt
      }
      if( exists("itransform",where=pc,inherits=FALSE) ){
        if( pc$itransform ){ g_th0 = t(dtheta0) %*% g_th0 }
      }
      gr_th0 = gr_th0 + g_th0
      if( "lp_theta0" %in% names(pc) ){
        dtheta0_b = rep(1,pc$ntheta0)
        if( length(ith0_bds[[1]]) > 0 ){
          dtheta0_b[ith0_bds[[1]]] = pc$dnotExp(theta0[ith0_bds[[1]]])
        }
        if( length(ith0_bds[[2]]) > 0 ){
          dtheta0_b[ith0_bds[[2]]] = -pc$dnotExp(theta0[ith0_bds[[2]]])
        }
        if( length(ith0_bds[[3]]) > 0 ){
          dtheta0_b[ith0_bds[[3]]] = pc$dnotExp(theta0[ith0_bds[[3]]])*
                                     pc$theta0_range/(1+tau)^2
        }
      }
      theta0 = pc$transform(theta0, pc=pc)
    }
    if( "lp_theta0" %in% names(pc) ){
      Arg = "(theta0,pc)"
      glp_theta0_call = paste("pc$glp$",pc$lp_theta0$g,Arg,sep="")
      # evaluate gradient for new event inference parameters
      g_th0 = eval(parse(text=glp_theta0_call))
      if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
        g_th0 = g_th0 * dtheta0_b
      }
      if( exists("itransform",where=pc,inherits=FALSE) ){
        if( pc$itransform ){ g_th0 = t(dtheta0) %*% g_th0 }
      }
      gr_th0 = gr_th0 + g_th0
    }
    x = x[-(1:pc$ntheta0)]
  } else { gr_th0 = NULL }
  # extract errors-in-variables yield parameters
  if( exists("eiv",where=pc,inherits=FALSE) && pc$eiv ){
    w_eiv = x[1:pc$nsource]
    # FGSN parameters
    tfgsn = pc$gdfgsn(w_eiv,fgsn_all[1],sqrt(lam2),fgsn_all[-(1:2)])
    gr_fgsn = tfgsn$hp
    # include lam2 prior and jacobian
    gr_fgsn[2] = gr_fgsn[2] + 0.5/lam2 - (0.5+1) + 1
    # include omega prior
    gr_fgsn[-(1:2)] = gr_fgsn[-(1:2)] - (fgsn_all[-(1:2)] - 0)/(10^2)
    # errors-in-variables yields
    gr_eiv = tfgsn$w
    x = x[-(1:pc$nsource)]
  } else { gr_fgsn = NULL; gr_eiv = NULL; }
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
  # extract level 1 variance components
  if( pc$pvc_1 > 0 ){
    vc1_all = x[1:pc$pvc_1]
    x = x[-(1:pc$pvc_1)]
  }
  # extract level 2 variance components
  if( pc$pvc_1 > 0 && pc$pvc_2 > 0 ){
    vc2_all = x[1:pc$pvc_2]
    x = x[-(1:pc$pvc_2)]
  }

  # log-prior gradient due to calibration data and
  # (optionally) new event data

  # initialize gradient vectors
  # common forward model parameters
  gr_beta0 = NULL
  # emplacement condition dependent forward model parameters
  gr_betat = NULL
  # level 1 variance components
  gr_vc1 = NULL
  # level 2 variance components
  gr_vc2 = NULL
  # observational error covariance parameters
  gr_eps = NULL
  # variance component prior scale parameters
  gr_A = NULL

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
      # level 1
      pvc_1 = sum(pc$h[[hh]]$pvc_1)
      vc_1 = exp(vc1_all[1:pvc_1])
      vc1_all = vc1_all[-(1:pvc_1)]
      if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                              pc$h[[hh]]$pvc_2 > 0) ){
        # level 2
        pvc_2 = sum(pc$h[[hh]]$pvc_2)
        vc_2 = exp(vc2_all[1:pvc_2])
        vc2_all = vc2_all[-(1:pvc_2)]
      }
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
    sdoe_h = sqrt(voe_h)
    ISd_h = Diagonal(Rh,1/sdoe_h)
    C_h = ISd_h %*% Sigma_h %*% ISd_h
    ldetC_h = 2*sum(ell_h)-slvoe_h

    # initialize gradient vectors forward model parameters
    if( ptbeta > 0 ){
      g_betat = vector("list",Th)
    }

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
          glp_beta_call = paste("pc$glp$",pc$h[[hh]]$lp_beta$g[rr],
                                Arg,sep="")
          # evaluate gradient for common model parameters
          gr_beta0 = c(gr_beta0,eval(parse(text=glp_beta_call)))
        } else { gr_beta0 = c(gr_beta0,numeric(pc$h[[hh]]$pbeta[rr])) }
      }
      if( ptbeta > 0 ){
        for( tt in 1:Th ){
          if( pc$h[[hh]]$pbetat[[tt]][rr] > 0 ){
            if( "lp_betat" %in% pnames &&
                !is.null(pc$h[[hh]]$lp_betat[[tt]]) ){
              st_betatr = 0
              if( rr > 1 ){
                st_betatr = sum(pc$h[[hh]]$pbetat[[tt]][1:(rr-1)])
              }
              betatr = betat[[tt]][st_betatr+
                                   (1:pc$h[[hh]]$pbetat[[tt]][rr])]
              Arg = "(betatr,pc)"
              glp_betat_call = paste("pc$glp$",
                                     pc$h[[hh]]$lp_betat[[tt]]$g[rr],
                                     Arg,sep="")
              # evaluate gradient for emplacement condition
              # dependent parameters
              g_betat[[tt]] = c(g_betat[[tt]],
                                eval(parse(text=glp_betat_call)))
            } else {
              g_betat[[tt]] = c(g_betat[[tt]],
                                numeric(pc$h[[hh]]$pbetat[[tt]][rr]))
            }
          }
        }
      }

      # evaluate gradient for variance components
      if( pc$pvc_1 > 0 && pc$h[[hh]]$pvc_1[rr] > 0 ){
        st_vc1 = 0
        if( rr > 1 ){
          st_vc1 = sum(pc$h[[hh]]$pvc_1[1:(rr-1)])
        }
        vc_r = vc_1[st_vc1+(1:pc$h[[hh]]$pvc_1[rr])]
        if( !("A" %in% pnames0) ){ iA = iA+1 }
        # level 1 including jacobian
        gr_vc1 = c(gr_vc1,1/2-vc_r/(A[iA]^2+vc_r))
        # scale parameter
        if( !("A" %in% pnames0) ){
          g_A = pc$h[[hh]]$pvc_1[rr]-2*A[iA]^2*sum(1/(A[iA]^2+vc_r))
        }
        if( pc$pvc_2 > 0 && any(pc$h[[hh]]$pvc_1 > 0 &
                                pc$h[[hh]]$pvc_2 > 0) ){
          if( pc$h[[hh]]$pvc_2[rr] > 0 ){
            st_vc2 = 0 
            if( rr > 1 ){
              st_vc2 = sum(pc$h[[hh]]$pvc_2[1:(rr-1)])
            }
            vc_r = vc_2[st_vc2+(1:pc$h[[hh]]$pvc_2[rr])]
            # level 2 including jacobian
            gr_vc2 = c(gr_vc2,1/2-vc_r/(A[iA]^2+vc_r))
            # scale parameter
            if( !("A" %in% pnames0) ){
              g_A = g_A+pc$h[[hh]]$pvc_2[rr]-
                    2*A[iA]^2*sum(1/(A[iA]^2+vc_r))
            }
          }
        }
        if( !("A" %in% pnames0) ){
          # jacobian
          g_A = g_A + 1
          gr_A = c(gr_A,g_A)
        }
      }
    }

    # evaluate gradient for observational error parameters
    # jacobian
    ncpar = Rh*(Rh+1)/2
    Jac_c = Matrix(0,ncpar,ncpar)
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
    IJac_c = solve(Jac_c)
    gr_eps_mat = matrix(0,Rh,Rh)
    for( q in 1:Rh ){
      for( p in 1:q ){
        dJac_c = Matrix(0,ncpar,ncpar)
        if( q > 1 ){
          for( r1 in 1:(q-1) ){
            if( r1 >= p ){
              dJac_c[k_ij(r1,q),k_ij(1,r1):k_ij(r1,r1)]=
              C_h[r1,q]*L_h[p,q]*L_h[1:r1,r1]/voe_h[r1]/voe_h[q]-
              L_h[p,r1]*L_h[1:r1,r1]/sdoe_h[r1]^3/sdoe_h[q]-
              L_h[p,q]*L_h[1:r1,q]/sdoe_h[q]^3/sdoe_h[r1]
              dJac_c[k_ij(r1,q),k_ij(p,r1)]=
              dJac_c[k_ij(r1,q),k_ij(p,r1)]+1/sdoe_h[r1]/sdoe_h[q]
              dJac_c[k_ij(r1,q),k_ij(r1,r1)]=
              dJac_c[k_ij(r1,q),k_ij(r1,r1)]*L_h[r1,r1]
              dJac_c[k_ij(r1,q),k_ij(1,q):k_ij(r1,q)]=
              3*C_h[r1,q]*L_h[p,q]*L_h[1:r1,q]/voe_h[q]^2-
              (L_h[p,q]*L_h[1:r1,r1]+L_h[p,r1]*L_h[1:r1,q])/
              sdoe_h[q]^3/sdoe_h[r1]
              dJac_c[k_ij(r1,q),k_ij(p,q)]=
              dJac_c[k_ij(r1,q),k_ij(p,q)]-C_h[r1,q]/voe_h[q]
              dJac_c[k_ij(r1,q),k_ij(r1+1,q):k_ij(q,q)]=
              3*C_h[r1,q]*L_h[p,q]*L_h[(r1+1):q,q]/voe_h[q]^2-
              L_h[p,r1]*L_h[(r1+1):q,q]/sdoe_h[q]^3/sdoe_h[r1]
              dJac_c[k_ij(r1,q),k_ij(q,q)]=
              dJac_c[k_ij(r1,q),k_ij(q,q)]*L_h[q,q]
            } else {
              dJac_c[k_ij(r1,q),k_ij(1,r1):k_ij(r1,r1)]=
              C_h[r1,q]*L_h[p,q]*L_h[1:r1,r1]/voe_h[r1]/voe_h[q]-
              L_h[p,q]*L_h[1:r1,q]/sdoe_h[q]^3/sdoe_h[r1]
              dJac_c[k_ij(r1,q),k_ij(1,q):k_ij(r1,q)]=
              3*C_h[r1,q]*L_h[p,q]*L_h[1:r1,q]/voe_h[q]^2-
              L_h[p,q]*L_h[1:r1,r1]/sdoe_h[r1]/sdoe_h[q]^3                
              dJac_c[k_ij(r1,q),k_ij((r1+1),q):k_ij(q,q)]=
              3*C_h[r1,q]*L_h[p,q]*L_h[(r1+1):q,q]/voe_h[q]^2
              dJac_c[k_ij(r1,q),k_ij(p,q)]=
              dJac_c[k_ij(r1,q),k_ij(p,q)]-C_h[r1,q]/voe_h[q]
              if( p == q ){
                dJac_c[k_ij(r1,q),k_ij(1,r1):k_ij(r1-1,r1)]=
                dJac_c[k_ij(r1,q),k_ij(1,r1):k_ij(r1-1,r1)]*L_h[q,q]
                dJac_c[k_ij(r1,q),k_ij(r1,r1)]=
                dJac_c[k_ij(r1,q),k_ij(r1,r1)]*L_h[r1,r1]*L_h[q,q]
                dJac_c[k_ij(r1,q),k_ij(1,q):k_ij(q-1,q)]=
                dJac_c[k_ij(r1,q),k_ij(1,q):k_ij(q-1,q)]*L_h[q,q]
                dJac_c[k_ij(r1,q),k_ij(q,q)]=
                (dJac_c[k_ij(r1,q),k_ij(q,q)]-C_h[r1,q]/voe_h[q])*
                L_h[q,q]^2
              } else {
                dJac_c[k_ij(r1,q),k_ij(r1,r1)]=
                dJac_c[k_ij(r1,q),k_ij(r1,r1)]*L_h[r1,r1]
                dJac_c[k_ij(r1,q),k_ij(q,q)]=
                dJac_c[k_ij(r1,q),k_ij(q,q)]*L_h[q,q]
              }
            }
          }
        }
        if( q < Rh ){
          for( r2 in (q+1):Rh ){
            dJac_c[k_ij(q,r2),k_ij(1,q):k_ij(q,q)]=
            3*C_h[q,r2]*L_h[p,q]*L_h[1:q,q]/voe_h[q]^2-
            (L_h[p,r2]*L_h[1:q,q]+L_h[p,q]*L_h[1:q,r2])/
            sdoe_h[q]^3/sdoe_h[r2]
            dJac_c[k_ij(q,r2),k_ij(p,q)]=
            dJac_c[k_ij(q,r2),k_ij(p,q)]-C_h[q,r2]/voe_h[q]
            dJac_c[k_ij(q,r2),k_ij(1,r2):k_ij(q,r2)]=
            C_h[q,r2]*L_h[p,q]*L_h[1:q,r2]/voe_h[q]/voe_h[r2]-
            L_h[p,r2]*L_h[1:q,r2]/sdoe_h[r2]^3/sdoe_h[q]-
            L_h[p,q]*L_h[1:q,q]/sdoe_h[q]^3/sdoe_h[r2]
            dJac_c[k_ij(q,r2),k_ij(p,r2)]=
            dJac_c[k_ij(q,r2),k_ij(p,r2)]+1/sdoe_h[q]/sdoe_h[r2]
            dJac_c[k_ij(q,r2),k_ij(q+1,r2):k_ij(r2,r2)]=
            C_h[q,r2]*L_h[p,q]*L_h[(q+1):r2,r2]/voe_h[q]/voe_h[r2]-
            L_h[p,r2]*L_h[(q+1):r2,r2]/sdoe_h[r2]^3/sdoe_h[q]
            if( p == q ){
              dJac_c[k_ij(q,r2),k_ij(1,q):k_ij(q-1,q)]=
              dJac_c[k_ij(q,r2),k_ij(1,q):k_ij(q-1,q)]*L_h[q,q]
              dJac_c[k_ij(q,r2),k_ij(q,q)]=
              (dJac_c[k_ij(q,r2),k_ij(q,q)]-C_h[q,r2]/voe_h[q])*
              L_h[q,q]^2+L_h[q,r2]*L_h[q,q]/sdoe_h[q]/sdoe_h[r2]
              dJac_c[k_ij(q,r2),k_ij(1,r2):k_ij(r2-1,r2)]=
              dJac_c[k_ij(q,r2),k_ij(1,r2):k_ij(r2-1,r2)]*L_h[q,q]
              dJac_c[k_ij(q,r2),k_ij(r2,r2)]=
              dJac_c[k_ij(q,r2),k_ij(r2,r2)]*L_h[r2,r2]*L_h[q,q]
            } else {
              dJac_c[k_ij(q,r2),k_ij(q,q)]=
              dJac_c[k_ij(q,r2),k_ij(q,q)]*L_h[q,q]
              dJac_c[k_ij(q,r2),k_ij(r2,r2)]=
              dJac_c[k_ij(q,r2),k_ij(r2,r2)]*L_h[r2,r2]
            }
          }
        }
        dJac_c[k_ij(q,q),k_ij(p,q)]=2
        if( p == q ){
          dJac_c[k_ij(q,q),k_ij(q,q)]=
          (dJac_c[k_ij(q,q),k_ij(q,q)]+2)*L_h[q,q]^2
        }
        # gradient of observational error variances and correlations
        # log-prior
        g_eps = -2*pc$lp_corr$eta*L_h[p,q]/voe_h[q]
        if( p == q ){ g_eps = g_eps*L_h[q,q]+2*(pc$lp_corr$eta-1) }
        # add in gradient of log-Jacobian
        gr_eps_mat[p,q] = g_eps+sum(diag(IJac_c %*% dJac_c))
      }
    }
    gr_eps = c(gr_eps,diag(gr_eps_mat),gr_eps_mat[upper.tri(gr_eps_mat)]) 
    if( ptbeta > 0 ){
      for( tt in 1:Th ){
        if( !is.null(g_betat[[tt]]) ){
          gr_betat = c(gr_betat,g_betat[[tt]])
        }
      }
    }
  }
  return(c(gr_th0,gr_eiv,gr_beta0,gr_betat,gr_vc1,gr_vc2,gr_eps,gr_fgsn,
           gr_A))
}

gdfgsn = function(w,alpha,lambda,omega)
{
  M = length(omega)
  t1 = (w - alpha)/lambda
  t2 = NULL
  for (ii in 1:M){ t2 = cbind(t2, t1^(2*ii-1)) }
  t3 = t2 %*% omega
  t4 = dnorm(t3,log=TRUE)-pnorm(t3,log.p=TRUE)
  t4 = exp(t4)
  t5 = NULL
  for (ii in 1:M){ t5 = cbind(t5, (2*ii-1)*(t1^(2*ii-2))) }
  t6 = t5 %*% omega
  t7 = NULL
  for (ii in 1:M){ t7 = cbind(t7, t4*t2[,ii]) }
  return(list(hp = c(sum((t1 - t4*t6)/lambda),
                     sum((t1^2 - 1 - t1*t4*t6)/2),
                     apply(t7,2,sum)),
               w = (t4*t6 - t1)/lambda))
}

k_ij = function(i,j)
{
  i+choose(j,2)
}
