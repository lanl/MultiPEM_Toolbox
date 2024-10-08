########################################################################
#                                                                      #
# This file contains code for maximizing the log-likelihood of the     #
# new event parameters based on new event data, where the forward and  #
# error model parameters are fixed at values pre-estimated from the    #
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

calc_mle_0 = function(p_cal,gdir,adir,f0,nst=10,ncor=1,ci_lev=0.95,
                      igrad=TRUE,bfgs=TRUE,igrck=TRUE,t_cal=NULL,
                      g0=NULL,fopt_in=NULL,Xst=NULL,tst=NULL,
                      fopt_out=NULL,pl="multicore")
{
  #
  # FUNCTION INPUTS
  #

  # p_cal: environment storing all objects needed in characterization
  #        calculations
  # gdir: directory for general subroutines
  # adir: directory for application subroutines
  # f0: names of forward models for each response by phenomenology
  # nst: number of starting values for MLE optimization
  # ncor: number of cores for MLE optimization
  # ci_lev: confidence interval levels for calibration and new event
  #         parameter inference
  # igrad: forward model gradients provided (TRUE/FALSE)
  # bfgs: MLE optimization uses BFGS methods (TRUE/FALSE)
  # igrck: Likelihood function gradient check (TRUE/FALSE)
  # t_cal: object used if bounds supplied to MLE optimization
  # g0: names of forward model gradients for each response by
  #     phenomenology
  # fopt_in: location of input R data file with starting values for
  #          MLE optimization
  # Xst: matrix of starting values for MLE optimization if not
  #      generated by this function
  # tst: vector of starting values for new event parameters in
  #      MLE optimization
  # fopt_out: location to write output R data file with results of
  #           MLE optimization
  # pl: strategy for running parallel jobs (see help for plan()
  #     function in future package)

  #
  # END FUNCTION INPUTS
  #

  #
  # SOURCE SUPPORTING R FUNCTIONS
  #

  # Global
  # log-likelihood functions
  source(paste(gdir,"/log_likelihood_0.r",sep=""),local=TRUE)
  if( igrad ){
    # gradient of log-likelihood functions
    source(paste(gdir,"/glog_likelihood_0.r",sep=""),local=TRUE)
    # R function to calculate information matrix for new event
    # parameters
    source(paste(gdir,"/info_likelihood_0.r",sep=""),local=TRUE)
  } else {
    # R function to calculate observed information matrix
    # for new event parameters
    source(paste(gdir,"/observed_information_0.r",sep=""),local=TRUE)
  }

  # Application
  # forward models
  source(paste(adir,"/forward_0.r",sep=""),local=TRUE)
  for( hh in 1:p_cal$H ){
    ufn = unique(f0[[hh]])
    lfn = length(ufn)
    for( qq in 1:lfn ){
      eval(parse(text=paste("p_cal$ffm$",ufn[qq]," = ",ufn[qq],sep="")))
    }
  }
  if( igrad ){
    # gradient of forward models
    source(paste(adir,"/jacobian_0.r",sep=""),local=TRUE)
    if( is.null(g0) ){
      stop("Gradient function names must be provided.")
    }
    for( hh in 1:p_cal$H ){
      ugn = unique(g0[[hh]])
      lgn = length(ugn)
      for( qq in 1:lgn ){
        eval(parse(text=paste("p_cal$gfm$",ugn[qq]," = ",ugn[qq],
                              sep="")))
      }
    }
  }
  # print summary statistics
  source(paste(adir,"/print_sumstats_0.r",sep=""),local=TRUE)

  #
  # END SOURCE SUPPORTING R FUNCTIONS
  #

  #
  # USER SPECIFIED FIELDS
  #

  # Names of forward models and gradients
  for( hh in 1:p_cal$H ){
    p_cal$h[[hh]]$f0 = f0[[hh]]
    if( igrad ){
      if( is.null(g0) ){
        stop("Gradient function names must be provided.")
      }
      p_cal$h[[hh]]$g0 = g0[[hh]]
    }
  }

  #
  # END USER SPECIFIED FIELDS
  #

  #
  # ADDITIONAL QUANTITIES USED IN CALCULATIONS
  #

  # statistical model functions
  if( !igrad ){ gll_0 = NULL }
  p_cal$ll_0 = ll_0; p_cal$gll_0 = gll_0;
  if( !igrad ){ p_cal$obs_info_0 = obs_info_0 }
  # support functions
  p_cal$print_ss_0 = print_ss_0

  # Specifications for MLE
  # Starting value for optimization
  tlb = p_cal$theta0_bounds[,1]
  tub = p_cal$theta0_bounds[,2]
  if( is.null(Xst) ){
    if( !is.null(fopt_in) ){
      opt = vector("list",1)
      opt[[1]] = readRDS(fopt_in)
    }
    Xst = NULL
    for( ii in 1:nst ){
      if( ii > 1 ){
        xst = runif(p_cal$ntheta0, min=-2, max=2)
        itheta0_bounds = vector("list",3)
        if( exists("itheta0_bounds",where=p_cal,inherits=FALSE) ){
          itheta0_bounds=p_cal$itheta0_bounds
        } else if( !is.null(t_cal) ){
          itheta0_bounds=t_cal$itheta0_bounds
        }   
        if( length(itheta0_bounds[[1]]) > 0 ){
          for( jj in itheta0_bounds[[1]] ){
            xst[jj] = runif(1, min=tlb[jj], max=tlb[jj]+4)
          }
        }
        if( length(itheta0_bounds[[2]]) > 0 ){
          for( jj in itheta0_bounds[[2]] ){
            xst[jj] = runif(1, min=tub[jj]-4,max=tub[jj])
          }
        }
        if( length(itheta0_bounds[[3]]) > 0 ){
          for( jj in itheta0_bounds[[3]] ){
            xst[jj] = runif(1, min=tlb[jj],max=tub[jj])
          }
        }
      } else {
        xst = numeric(p_cal$ntheta0)
        if( !is.null(tst) ){ xst = tst }
      }
      if( !p_cal$opt_B ){
        if( exists("itheta0_bounds",where=p_cal,inherits=FALSE) ){
          xst = p_cal$inv_transform(xst,pc=p_cal)
        }
        if( exists("itransform",where=p_cal,inherits=FALSE) ){
          if( p_cal$itransform ){ xst = p_cal$inv_tau(xst,pc=p_cal) }
        }
      }
      if( !is.null(fopt_in) && ii == 1 ){
        if( is.null(tst) ){ xst = opt[[1]]$theta0 }
      }
      Xst = rbind(Xst, xst)
    }
  } else { Xst = Xst }
  p_cal$Xst = Xst

  #
  # END ADDITIONAL QUANTITIES USED IN CALCULATIONS
  #

  #
  # MAXIMUM LIKELIHOOD CALCULATION
  #

  # Maximize log likelihood to find MLE
  p_cal$tlb = tlb; p_cal$tub = tub;
  for( ii in 1:nst ){
    xst = Xst[ii,]
    while( is.infinite(ll_0(xst, pc=p_cal)) ){
      xst = calc_xst(p_cal,t_cal)
    }
    Xst[ii,] = xst
  }
  p_cal$calc_xst = calc_xst
  # parallel optimization using R package "future"
  require(doFuture)
  if( ncor == 1 ){ pl = "sequential"; plan(pl);
  } else {
    if( ncor > nst ){ ncor = nst }
    if( pl != "sequential" ){ plan(pl,workers=ncor)
    } else { plan(pl) }
  }
  cmax = -Inf
  ptm = proc.time()
  fCatch = foreach( qq = 1:nst ) %dofuture% {
             ll_opt(Xst[qq,],bfgs,p_cal,t_cal)
           } %seed% TRUE
  plan(sequential)
  print("Run time:")
  print(proc.time() - ptm)
  cat("\n")
  Mll = NULL
  for( qq in 1:nst ){
    if( is(fCatch[[qq]]$value,"list") ){
      Mle = fCatch[[qq]]$value
      Mll = c(Mll,Mle$value)
    } else { Mll = c(Mll,-Inf) }
  }
  imax = which( Mll == max(Mll) ); imax = imax[1];
  if( Mll[imax] > cmax ){
    cmax = Mll[imax]
    mle = fCatch[[imax]]$value
  }
  qq = 0; irel_tol = 0;
  while( qq <= 5 ){
    cmax = mle$value
    if( bfgs ){
      if( p_cal$opt_B ){
        mle = optim(mle$par, fn=ll_0, gr=gll_0, pc=p_cal,
                    method="L-BFGS-B",
                    lower=tlb, upper=tub,
                    control=list(fnscale = -1,maxit=1000))
      } else {
        mle = optim(mle$par, fn=ll_0, gr=gll_0, pc=p_cal,
                    method="BFGS",
                    control=list(fnscale = -1,maxit=1000))
      }
    } else {
      if( p_cal$opt_B ){
        mle = optim(mle$par, fn=ll_0, pc=p_cal,
                    lower=tlb, upper=tub,
                    control=list(fnscale = -1,maxit=10000))
      } else {
        mle = optim(mle$par, fn=ll_0, pc=p_cal,
                    control=list(fnscale = -1,maxit=10000))
      }
    }
    rel_tol = abs(mle$value - cmax)/abs(cmax)
    if( rel_tol <= 1.e-8 ){ irel_tol = irel_tol+1
    } else { qq = qq+1 }
    if( irel_tol == 2 ){ break }
  }
  # compute Hessian to obtain observed Fisher information matrix
  # in applications where gradients are not available
  if( !igrad ){
    if( p_cal$opt_B ){
      mle = optim(mle$par, fn=ll_0, pc=p_cal,
                  lower=tlb, upper=tub,
                  control=list(fnscale = -1,maxit=10000),
                  hessian=TRUE)
    } else {
      mle = optim(mle$par, fn=ll_0, pc=p_cal,
                  control=list(fnscale = -1,maxit=10000),
                  hessian=TRUE)
    }
  }

  # Print convergence status of log-likelihood optimization
  print("MLE CONVERGENCE STATUS")
  cat("\n")
  print(mle$convergence)
  print(irel_tol)

  # Additional information if convergence obtained
  if( mle$convergence == 0 ){
    p_cal$mle = unname(mle$par)
    opt_mle = p_cal$opt_mle
    opt_mle$theta0 = mle$par
    if( !is.null(fopt_out) ){ saveRDS(opt_mle, file=fopt_out) }
    # Calculate asymptotic covariance matrix
    # for new event inference parameters
    if( igrad ){
      p_cal$rapid = TRUE
      p_cal$Sigma_mle_0 = info_ll_0(opt_mle, p_cal)
    } else {
      if( p_cal$opt_B ){
        if( is.null(t_cal) ){
          stop(paste("List t_cal must be provided for bounded",
                     " optimization.",sep=""))
        }
        p_cal$Sigma_mle_0 = obs_info_0(p_cal, mle, imle=TRUE,
                                       t_cal=t_cal)
      } else {
        p_cal$Sigma_mle_0 = obs_info_0(p_cal, mle, imle=TRUE)
      }
    }
    # Print MLE
    print("MAXIMUM LIKELIHOOD SUMMARY")
    cat("\n")
    # Indicator of prior distribution parameters
    p_cal$iPrior = FALSE
    xmle = mle$par
    if( p_cal$ncalp > 0 ){ xmle = c(xmle,p_cal$mle_calp) }
    p_cal = print_ss_0(xmle, p_cal, ci=ci_lev)
  }

  #
  # END MAXIMUM LIKELIHOOD CALCULATION
  #

  #
  # CHECK LOG-LIKELIHOOD GRADIENT CALCULATIONS
  #

  if( igrck && igrad ){
    # numerical differentiation package
    require(numDeriv)
    print("CHECK LOG-LIKELIHOOD GRADIENTS")
    cat("\n")
    xnom = unname(mle$par)
    Catch = p_cal$tryCatch.W.E(gll_0(xnom,p_cal))
    if( is(Catch$value,"numeric") ){ gll_0_nom = Catch$value
    } else {
      gll_0_nom = NaN
      print("Error computing analytic gradient at MLE.")
    }
    Catch = p_cal$tryCatch.W.E(grad(ll_0, xnom,
                                    method.args=list(r=6),pc=p_cal))
    if( is(Catch$value,"numeric") ){ gll_0_nom_num = Catch$value
    } else {
      gll_0_nom_num = NaN
      print("Error computing numerical gradient at MLE.")
    }
    diff <- gll_0_nom - gll_0_nom_num
    print("Analytic gradient")
    print(gll_0_nom)
    print("Numerical gradient")
    print(gll_0_nom_num)
    print("Difference")
    print(c(min(diff),max(diff)))
    cat("\n")

    nsamp = 5
    qq = 1
    while( qq <= nsamp ){
      xx = xnom + rnorm(p_cal$ntheta0,mean=0,sd=0.1)
      Catch = p_cal$tryCatch.W.E(gll_0(xx,p_cal))
      if( is(Catch$value,"numeric") ){ gll_0_xx = Catch$value
      } else { next }
      Catch = p_cal$tryCatch.W.E(grad(ll_0, xx,
                                      method.args=list(r=6),pc=p_cal))
      if( is(Catch$value,"numeric") ){
        qq = qq+1
        gll_0_xx_num = Catch$value
      } else { next }
      diff <- gll_0_xx - gll_0_xx_num
      print("Analytic gradient")
      print(gll_0_xx)
      print("Numerical gradient")
      print(gll_0_xx_num)
      print("Difference")
      print(c(min(diff),max(diff)))
      cat("\n")
    }
  }

  #
  # END CHECK LOG-LIKELIHOOD GRADIENT CALCULATIONS
  #

  return(p_cal)
}

calc_xst = function(p_cal,t_cal)
{
  xst = runif(p_cal$ntheta0, min=-2, max=2)
  tlb = p_cal$theta0_bounds[,1]
  tub = p_cal$theta0_bounds[,2]
  itheta0_bounds = vector("list",3)
  if( exists("itheta0_bounds",where=p_cal,inherits=FALSE) ){
    itheta0_bounds=p_cal$itheta0_bounds
  } else if( !is.null(t_cal) ){
    itheta0_bounds=t_cal$itheta0_bounds                            
  }   
  if( length(itheta0_bounds[[1]]) > 0 ){
    for( jj in itheta0_bounds[[1]] ){
      xst[jj] = runif(1, min=tlb[jj],max=tlb[jj]+4)                
    }                                                              
  }
  if( length(itheta0_bounds[[2]]) > 0 ){
    for( jj in itheta0_bounds[[2]] ){
      xst[jj] = runif(1, min=tub[jj]-4,max=tub[jj])
    }             
  }
  if( length(itheta0_bounds[[3]]) > 0 ){
    for( jj in itheta0_bounds[[3]] ){
      xst[jj] = runif(1, min=tlb[jj],max=tub[jj])
    }
  }
  if( !p_cal$opt_B ){
    if( exists("itheta0_bounds",where=p_cal,inherits=FALSE) ){
      xst = p_cal$inv_transform(xst,pc=p_cal)
    }
    if( exists("itransform",where=p_cal,inherits=FALSE) ){
      if( p_cal$itransform ){ xst = p_cal$inv_tau(xst,pc=p_cal) }
    }
  }
  return(xst)
}

ll_opt = function(xst,bfgs,pc,tc)
{
  maxiter = 5
  ii = 1
  while( ii <= maxiter ){
    if( bfgs ){
      if( pc$opt_B ){
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$ll_0,
                gr=pc$gll_0, pc=pc, method="L-BFGS-B",
                lower=pc$tlb, upper=pc$tub,
                control=list(fnscale = -1,maxit=1000)))
      } else {
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$ll_0,
                gr=pc$gll_0, pc=pc, method="BFGS",
                control=list(fnscale = -1,maxit=1000)))
      }
    } else {
      if( pc$opt_B ){
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$ll_0,
                pc=pc, lower=pc$tlb, upper=pc$tub,
                control=list(fnscale = -1,maxit=10000)))
      } else {
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$ll_0,
                pc=pc, control=list(fnscale = -1,maxit=10000)))
      }
    }
    if( is(Catch$value,"list") ){
      return(Catch)
    } else {
      xst = pc$calc_xst(pc,tc)
      while( is.infinite(pc$ll_0(xst, pc=pc)) ){
        xst = pc$calc_xst(pc,tc)
      }
      ii = ii+1
      if( ii > maxiter ){ return(Catch) } else { next }
    }
  }
}
