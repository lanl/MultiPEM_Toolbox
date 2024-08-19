########################################################################
#                                                                      #
# This file contains code for maximizing and sampling the log-         #
# posterior of the new event parameters based on new event data, where #
# the forward and error model parameters are fixed at values           #
# pre-estimated from the calibration data.                             #
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

calc_bayes_0 = function(p_cal,gdir,adir,nst=10,nburn=10000,nmcmc=20000,
                        nthin=1,ncor_map=1,ncor_mc=1,igrad=TRUE,
                        igrck_pr=TRUE,igrck_po=TRUE,bfgs=TRUE,
                        itpr=FALSE,fpr_t=NULL,fgpr_t=NULL,imcmc="FME",
                        pl="multicore",ncor_smc=NULL,lb_smc=NULL,
                        ub_smc=NULL,t_cal=NULL)
{
  #
  # FUNCTION INPUTS
  #

  # p_cal: environment storing all objects needed in characterization
  #        calculations
  # gdir: directory for general subroutines
  # adir: directory for application subroutines
  # nst: number of starting values for MAP optimization
  # nburn: number of MCMC burn-in samples (per imputation)
  # nmcmc: number of MCMC production samples (per imputation)
  # nthin: posterior sample thinning rate (per imputation)
  # ncor_map: number of cores for MAP optimization
  # ncor_mc: number of cores for parallel MCMC chains or multiple
  #          imputation
  # igrad: forward model gradients provided (TRUE/FALSE)
  # igrck_pr: Prior distribution gradient check (TRUE/FALSE)
  # igrck_po: Posterior distribution gradient check (TRUE/FALSE)
  # bfgs: MAP optimization uses BFGS methods (TRUE/FALSE)
  # itpr: prior distributions provided for new event
  #       parameters (TRUE/FALSE)
  # fpr_t: location of functions computing log-prior density for
  #        new event parameters
  # fgpr_t: location of functions computing gradients of log-prior
  #         density for new event parameters
  # imcmc: MCMC algorithm (current options:  "RAM", "FME", "NUTS",
  #                                          "SMC")
  # pl: strategy for running parallel jobs (see help for plan()
  #     function in future package)
  # ncor_smc: number of cores for inner parallelism of SMC algorithm
  # lb_smc: lower bounds of new event parameters for SMC sampling
  # ub_smc: upper bounds of new event parameters for SMC sampling
  # t_cal: object used if bounds supplied to MLE optimization

  #
  # END FUNCTION INPUTS
  #

  #   
  # SOURCE SUPPORTING R FUNCTIONS
  #     

  if( !exists("obs_info_0",where=p_cal,inherits=FALSE) ){
    # R function to calculate observed information matrix
    source(paste(gdir,"/observed_information_0.r",sep=""),local=TRUE)
    p_cal$obs_info_0 = obs_info_0
  }

  #   
  # END SOURCE SUPPORTING R FUNCTIONS
  #

  #
  # BAYESIAN ANALYSIS
  #

  if( itpr ){
    # prior distributions for new event parameters
    if( is.null(fpr_t) ){
      source(paste(adir,"/lp_0.r",sep=""),local=TRUE)
    } else { source(fpr_t,local=TRUE) }
    if( "lp_theta0" %in% names(p_cal) ){
      ufn = p_cal$lp_theta0$f
      eval(parse(text=paste("p_cal$flp$",ufn," = ",ufn,sep="")))
    }
    # gradients
    if( igrad ){
      if( is.null(fgpr_t) ){
        source(paste(adir,"/glp_0.r",sep=""),local=TRUE)
      } else { source(fgpr_t,local=TRUE) }
      if( "lp_theta0" %in% names(p_cal) ){
        ugn = p_cal$lp_theta0$g
        eval(parse(text=paste("p_cal$glp$",ugn," = ",ugn,sep="")))
      }
    }
  }

  # starting values for parameters in MAP optimization
  tlb = p_cal$theta0_bounds[,1]
  tub = p_cal$theta0_bounds[,2]
  Xst = matrix(p_cal$mle,nrow=1)
  if( nst != nrow(p_cal$Xst) ){ nst = nrow(p_cal$Xst) }
  if( nst > 1 ){ Xst = rbind(Xst,p_cal$Xst[2:nst,,drop=FALSE]) }

  # Source log-prior and gradient of the log-prior
  source(paste(gdir,"/log_prior_0.r",sep=""),local=TRUE)
  if( igrad ){ source(paste(gdir,"/glog_prior_0.r",sep=""),local=TRUE) }

  # Attach prior functions needed for Bayesian computation
  if( !igrad ){ glprior_0 = NULL }
  p_cal$lprior_0 = lprior_0; p_cal$glprior_0 = glprior_0;

  # Source log-posterior and gradient of the log-posterior
  source(paste(gdir,"/log_posterior_0.r",sep=""),local=TRUE)
  if( igrad ){
    source(paste(gdir,"/glog_posterior_0.r",sep=""),local=TRUE)
  }

  # Attach posterior functions needed for Bayesian computation
  if( !igrad ){ glpost_0 = NULL }
  p_cal$lpost_0 = lpost_0; p_cal$glpost_0 = glpost_0;

  # Maximize log posterior to find MAP estimate
  p_cal$tlb = tlb; p_cal$tub = tub;
  for( ii in 1:nst ){
    xst = Xst[ii,]
    while( is.infinite(lpost_0(xst, pc=p_cal)) ){
      xst = calc_xst(p_cal,t_cal)
    }
    Xst[ii,] = xst
  }
  p_cal$calc_xst = calc_xst
  # parallel optimization using R package "future"
  require(doFuture)
  if( ncor_map == 1 ){ pl = "sequential"; plan(pl);
  } else {
    if( ncor_map > nst ){ ncor_map = nst }
    if( pl != "sequential" ){ plan(pl,workers=ncor_map)
    } else { plan(pl) }
  }
  cmax = -Inf
  ptm = proc.time()
  fCatch = foreach( qq = 1:nst ) %dofuture% {
             lpo_opt(Xst[qq,],bfgs,p_cal,t_cal)
           } %seed% TRUE
  plan(sequential)
  print("Run time:")
  print(proc.time() - ptm)
  cat("\n")
  Mlp = NULL
  for( qq in 1:nst ){
    if( is(fCatch[[qq]]$value,"list") ){
      Map = fCatch[[qq]]$value
      Mlp = c(Mlp,Map$value)
    } else { Mlp = c(Mlp,-Inf) }
  }
  imax = which( Mlp == max(Mlp) ); imax = imax[1];
  if( Mlp[imax] > cmax ){
    cmax = Mlp[imax]
    map = fCatch[[imax]]$value
  }
  qq = 0; irel_tol = 0;
  while( qq <= 5 ){
    cmax = map$value
    if( bfgs ){
      if( p_cal$opt_B ){
        map = optim(map$par, fn=lpost_0, gr=glpost_0, pc=p_cal,
                    method="L-BFGS-B",
                    lower=tlb, upper=tub,
                    control=list(fnscale = -1,maxit=1000))
      } else {
        map = optim(map$par, fn=lpost_0, gr=glpost_0, pc=p_cal,
                    method="BFGS",
                    control=list(fnscale = -1,maxit=1000))
      }
    } else {
      if( p_cal$opt_B ){
        map = optim(map$par, fn=lpost_0, pc=p_cal,
                    lower=tlb, upper=tub,
                    control=list(fnscale = -1,maxit=10000))
      } else {
        map = optim(map$par, fn=lpost_0, pc=p_cal,
                    control=list(fnscale = -1,maxit=10000))
      }
    }
    rel_tol = abs(map$value - cmax)/abs(cmax)
    if( rel_tol <= 1.e-8 ){ irel_tol = irel_tol+1
    } else { qq = qq+1 }
    if( irel_tol == 2 ){ break }
  }
  if( imcmc %in% c("RAM","FME") ){
    # inverse Hessian for starting proposal distribution in
    # adaptive MCMC
    if( bfgs ){
      if( p_cal$opt_B ){
        map = optim(map$par, fn=lpost_0, gr=glpost_0, pc=p_cal,
                    method="L-BFGS-B",
                    lower=tlb, upper=tub,
                    control=list(fnscale = -1,maxit=1000),
                    hessian=TRUE)
      } else {
        map = optim(map$par, fn=lpost_0, gr=glpost_0, pc=p_cal,
                    method="BFGS",
                    control=list(fnscale = -1,maxit=1000),
                    hessian=TRUE)
      }
    } else {
      if( p_cal$opt_B ){
        map = optim(map$par, fn=lpost_0, pc=p_cal,
                    lower=tlb, upper=tub,
                    control=list(fnscale = -1,maxit=10000),
                    hessian=TRUE)
      } else {
        map = optim(map$par, fn=lpost_0, pc=p_cal,
                    control=list(fnscale = -1,maxit=10000),
                    hessian=TRUE)
      }
    }
    if( p_cal$opt_B ){
      if( is.null(t_cal) ){
        stop("List t_cal must be provided for bounded optimization.")
      }
      p_cal = p_cal$obs_info_0(p_cal, map, imle=FALSE, t_cal=t_cal)
    } else { p_cal = p_cal$obs_info_0(p_cal, map, imle=FALSE) }
    lambda = (2.38)^2/p_cal$ntheta0
    p_cal$IHess = lambda*p_cal$IHess
  }

  # Print convergence status of log-posterior optimization
  print("MAP CONVERGENCE STATUS")
  cat("\n")
  print(map$convergence)
  print(irel_tol)

  # Additional information if convergence obtained
  if( map$convergence == 0 ){
    p_cal$map = unname(map$par)
    # Print MAP
    print("MAXIMUM A POSTERIORI SUMMARY")
    cat("\n")
    # Indicator of prior distribution parameters
    p_cal$iPrior = TRUE
    p_cal = p_cal$print_ss_0(map$par, p_cal)
  }

  # Set up unbounded parameters for posterior sampling
  if( p_cal$opt_B ){
    if( is.null(t_cal) ){
      stop("List t_cal must be provided for bounded optimization.")
    }
    p_cal$itheta0_bounds = t_cal$itheta0_bounds
    map$par = p_cal$inv_transform(map$par,pc=p_cal)
    if( exists("itransform",where=p_cal,inherits=FALSE) ){
      if( p_cal$itransform ){
        map$par = p_cal$inv_tau(map$par,pc=p_cal)
      }
    }
  }

  # Check log-prior gradient calculations
  if( igrck_pr && igrad ){
    # numerical differentiation package
    require(numDeriv)
    print("CHECK LOG-PRIOR GRADIENTS")
    cat("\n")
    xnom = unname(map$par)
    Catch = p_cal$tryCatch.W.E(glprior_0(xnom,p_cal))
    if( is(Catch$value,"numeric") ){ glpr_0_nom = Catch$value
    } else {
      glpr_0_nom = NaN
      print("Error computing analytic gradient at MAP.")
    }
    Catch = p_cal$tryCatch.W.E(grad(lprior_0, xnom,
                                    method.args=list(r=6),pc=p_cal))
    if( is(Catch$value,"numeric") ){ glpr_0_nom_num = Catch$value
    } else {
      glpr_0_nom_num = NaN
      print("Error computing numerical gradient at MAP.")
    }
    diff <- glpr_0_nom - glpr_0_nom_num
    print("Analytic gradient")
    print(glpr_0_nom)
    print("Numerical gradient")
    print(glpr_0_nom_num)
    print("Difference")
    print(c(min(diff),max(diff)))
    cat("\n")

    nsamp = 5
    qq = 1
    while( qq <= nsamp ){
      xx = xnom + rnorm(p_cal$ntheta0,mean=0,sd=0.1)
      Catch = p_cal$tryCatch.W.E(glprior_0(xx,p_cal))
      if( is(Catch$value,"numeric") ){ glpr_0_xx = Catch$value
      } else { next }
      Catch = p_cal$tryCatch.W.E(grad(lprior_0, xx,
                                      method.args=list(r=6),pc=p_cal))
      if( is(Catch$value,"numeric") ){
        qq = qq+1
        glpr_0_xx_num = Catch$value
      } else { next }
      diff <- glpr_0_xx - glpr_0_xx_num
      print("Analytic gradient")
      print(glpr_0_xx)
      print("Numerical gradient")
      print(glpr_0_xx_num)
      print("Difference")
      print(c(min(diff),max(diff)))
      cat("\n")
    }
  }

  # Check log-posterior gradient calculations
  if( igrck_po && igrad ){
    # numerical differentiation package
    require(numDeriv)
    print("CHECK LOG-POSTERIOR GRADIENTS")
    cat("\n")
    xnom = unname(map$par)
    Catch = p_cal$tryCatch.W.E(glpost_0(xnom,p_cal))
    if( is(Catch$value,"numeric") ){ glpo_0_nom = Catch$value
    } else {
      glpo_0_nom = NaN
      print("Error computing analytic gradient at MAP.")
    }
    Catch = p_cal$tryCatch.W.E(grad(lpost_0, xnom,
                                    method.args=list(r=6),pc=p_cal))
    if( is(Catch$value,"numeric") ){ glpo_0_nom_num = Catch$value
    } else {
      glpo_0_nom_num = NaN
      print("Error computing numerical gradient at MAP.")
    }
    diff <- glpo_0_nom - glpo_0_nom_num
    print("Analytic gradient")
    print(glpo_0_nom)
    print("Numerical gradient")
    print(glpo_0_nom_num)
    print("Difference")
    print(c(min(diff),max(diff)))
    cat("\n")

    nsamp = 5
    qq = 1
    while( qq <= nsamp ){
      xx = xnom + rnorm(p_cal$ntheta0,mean=0,sd=0.1)
      Catch = p_cal$tryCatch.W.E(glpost_0(xx,p_cal))
      if( is(Catch$value,"numeric") ){ glpo_0_xx = Catch$value
      } else { next }
      Catch = p_cal$tryCatch.W.E(grad(lpost_0, xx,
                                      method.args=list(r=6),pc=p_cal))
      if( is(Catch$value,"numeric") ){
        qq = qq+1
        glpo_0_xx_num = Catch$value
      } else { next }
      diff <- glpo_0_xx - glpo_0_xx_num
      print("Analytic gradient")
      print(glpo_0_xx)
      print("Numerical gradient")
      print(glpo_0_xx_num)
      print("Difference")
      print(c(min(diff),max(diff)))
      cat("\n")
    }
  }

  # number of multiple imputation samples
  # if nimp = 1, MLE used
  nimp = p_cal$nimp

  if (imcmc == "RAM"){
    # Robust Adaptive Metropolis sampling
    # parallel chains using R package "future"
    require(doFuture)
    # using R package "adaptMCMC"
    require(adaptMCMC)
    if( ncor_mc == 1 ){ pl = "sequential"; plan(pl);
    } else {
      if( ncor_mc > nimp && nimp > 1 ){ ncor_mc = nimp }
      if( pl != "sequential" ){ plan(pl,workers=ncor_mc)
      } else { plan(pl) }
    }
    if( nimp == 1 ){
      xi = NULL
      if( ncor_mc > 1 ){
        eps = 0.001
        for (qq in 1:ncor_mc){
          xi = rbind(xi, map$par*(1+runif(p_cal$ntheta0,-eps,eps)))
        }
      } else { xi = rbind(xi, map$par) }
      si = sample(10000,ncor_mc)
      ptm = proc.time()
      fram = lapply(1:ncor_mc, function(qq) { future(MCMC(lpost_0,
                               nburn+round(nmcmc/ncor_mc),init=xi[qq,],
                               scale = as.matrix(p_cal$IHess),
                               acc.rate=0.234,showProgressBar=FALSE,
                               pc=p_cal),seed=si[qq]) })
      ram = lapply(fram,value)
      plan(sequential)
      print("Run time:")
      print(proc.time() - ptm)
      cat("\n")
      # acceptance rate
      print("ACCEPTANCE RATES:")
      cat("\n")
      for( qq in 1:ncor_mc ){
        print(paste("Core ",qq,": ",ram[[qq]]$acceptance.rate,sep=""))
      }
      cat("\n")
      # extract post-burnin posterior samples
      mpi = NULL
      for (qq in 1:ncor_mc){
        tmpi = ram[[qq]]$samples
        mpi = rbind(mpi,as.matrix(tmpi[-(1:nburn),]))
      }
    } else {
      xi = NULL
      eps = 0.001
      for( qq in 1:nimp ){
        xi = rbind(xi, map$par*(1+runif(p_cal$ntheta0,-eps,eps)))
      }
      ptm = proc.time()
      fram = foreach( qq = 1:nimp ) %dofuture% {
               q_cal = p_cal$pc_0(p_cal$mpi[qq,], p_cal)
               MCMC(lpost_0,nburn+nmcmc,init=xi[qq,],
                    scale = as.matrix(p_cal$IHess),
                    acc.rate=0.234,showProgressBar=FALSE,
                    pc=q_cal)
             } %seed% TRUE
      plan(sequential)
      print("Run time:")
      print(proc.time() - ptm)
      cat("\n")
      # acceptance rate
      print("ACCEPTANCE RATES:")
      cat("\n")
      for( qq in 1:nimp ){
        print(paste("Imputation ",qq,": ",
                    fram[[qq]]$acceptance.rate,sep=""))
      }
      cat("\n")
      # extract post-burnin posterior samples
      mpi = NULL
      ipi = seq(1,nmcmc,by=nthin)
      for (qq in 1:nimp){
        tmpi = fram[[qq]]$samples
        tmpi = as.matrix(tmpi[-(1:nburn),])
        # thin posterior samples
        tmpi = tmpi[ipi,,drop=FALSE]
        mpi = rbind(mpi,tmpi)
      }
    }
  }
  if (imcmc == "FME"){
    # Constrained MCMC Sampler
    # parallel chains using R package "future"
    require(doFuture)
    # using R package "FME"
    require(FME)
    if( ncor_mc == 1 ){ pl = "sequential"; plan(pl);
    } else {
      if( ncor_mc > nimp && nimp > 1 ){ ncor_mc = nimp }
      if( pl != "sequential" ){ plan(pl,workers=ncor_mc)
      } else { plan(pl) }
    }
    if( nimp == 1 ){
      xi = NULL
      if( ncor_mc > 1 ){
        eps = 0.001
        for (qq in 1:ncor_mc){
          xi = rbind(xi, map$par*(1+runif(p_cal$ntheta0,-eps,eps)))
        }
      } else { xi = rbind(xi, map$par) }
      si = sample(10000,ncor_mc)
      ll_0_mod = function(x, pc=p_cal){ -2*pc$ll_0(x, pc) }
      lprior_0_mod = function(x, pc=p_cal){ -2*pc$lprior_0(x, pc) }
      ptm = proc.time()
      ffme = lapply(1:ncor_mc, function(qq) { future(modMCMC(
                            ll_0_mod,xi[qq,],pc=p_cal,
                            jump=as.matrix(p_cal$IHess),
                            lower=rep(-Inf,p_cal$ntheta0),
                            upper=rep(Inf,p_cal$ntheta0),
                            prior=lprior_0_mod,
                            niter=nburn+round(nmcmc/ncor_mc),
                            burninlength=nburn,updatecov=100,ntrydr=2,
                            verbose=FALSE),seed=si[qq]) })
      fme = lapply(ffme,value)
      plan(sequential)
      print("Run time:")
      print(proc.time() - ptm)
      cat("\n")
      # acceptance rate
      print("ACCEPTANCE RATES:")
      cat("\n")
      for( qq in 1:ncor_mc ){
        print(paste("Core ",qq,": ",fme[[qq]]$naccepted/
                                    (nburn+round(nmcmc/ncor_mc)),
                    sep=""))
      }
      cat("\n")
      # extract post-burnin posterior samples
      mpi = NULL
      for (qq in 1:ncor_mc){
        tmpi = fme[[qq]]$pars
        mpi = rbind(mpi,tmpi)
      }
      if (is.vector(mpi)) { mpi = matrix(mpi,ncol=1) }
    } else {
      xi = NULL
      eps = 0.001
      for( qq in 1:nimp ){
        xi = rbind(xi, map$par*(1+runif(p_cal$ntheta0,-eps,eps)))
      }
      ll_0_mod = function(x, pc=p_cal){ -2*pc$ll_0(x, pc) }
      lprior_0_mod = function(x, pc=p_cal){ -2*pc$lprior_0(x, pc) }
      ptm = proc.time()
      fme = foreach( qq = 1:nimp ) %dofuture% {
              q_cal = p_cal$pc_0(p_cal$mpi[qq,], p_cal)
              modMCMC(ll_0_mod,xi[qq,],pc=q_cal,
                      jump=as.matrix(p_cal$IHess),
                      lower=rep(-Inf,p_cal$ntheta0),
                      upper=rep(Inf,p_cal$ntheta0),
                      prior=lprior_0_mod,
                      niter=nburn+nmcmc,
                      burninlength=nburn,updatecov=100,ntrydr=2,
                      verbose=FALSE)
            } %seed% TRUE
      plan(sequential)
      print("Run time:")
      print(proc.time() - ptm)
      cat("\n")
      # acceptance rate
      print("ACCEPTANCE RATES:")
      cat("\n")
      for( qq in 1:nimp ){
        print(paste("Imputation ",qq,": ",fme[[qq]]$naccepted/
                    (nburn+nmcmc),sep=""))
      }
      cat("\n")
      # extract post-burnin posterior samples
      mpi = NULL
      ipi = seq(1,nmcmc,by=nthin)
      for (qq in 1:nimp){
        tmpi = as.matrix(fme[[qq]]$pars)
        # thin posterior samples
        tmpi = tmpi[ipi,,drop=FALSE]
        mpi = rbind(mpi,tmpi)
      }
    }
  }
  if (imcmc == "NUTS"){
    # Hamiltonian Monte Carlo No-U-Turn (NUTS) sampling
    # parallel chains using R package "future"
    require(doFuture)
    # using R code from https://github.com/kasparmartens/NUTS
    source(paste(gdir,"/helpers.r",sep=""),local=TRUE)
    source(paste(gdir,"/nuts.r",sep=""),local=TRUE)
    if( ncor_mc == 1 ){ pl = "sequential"; plan(pl);
    } else {
      if( ncor_mc > nimp && nimp > 1 ){ ncor_mc = nimp }
      if( pl != "sequential" ){ plan(pl,workers=ncor_mc)
      } else { plan(pl) }
    }
    if( nimp == 1 ){
      xi = NULL
      if( ncor_mc > 1 ){
        eps = 0.001
        for (qq in 1:ncor_mc){
          xi = rbind(xi, map$par*(1+runif(p_cal$ntheta0,-eps,eps)))
        }
      } else { xi = rbind(xi, map$par) }
      si = sample(10000,ncor_mc)
      if( !igrad ){
        stop("Gradient of log-posterior must be provided.")
      }
      ptm = proc.time()
      fham = lapply(1:ncor_mc, function(qq) { future(NUTS(xi[qq,],
                               f=lpost_0,grad_f=glpost_0,
                               n_iter=nburn+round(nmcmc/ncor_mc),
                               delta=0.8,verbose=FALSE),
                               seed=si[qq]) })
      ham = lapply(fham, value)
      plan(sequential)
      print("Run time:")
      print(proc.time() - ptm)
      cat("\n")
      # extract post-burnin posterior samples
      mpi = NULL
      for (qq in 1:ncor_mc){
        tmpi = ham[[qq]]
        mpi = rbind(mpi, as.matrix(tmpi[-(1:nburn),]))
      }
    } else {
      xi = NULL
      eps = 0.001
      for( qq in 1:nimp ){
        xi = rbind(xi, map$par*(1+runif(p_cal$ntheta0,-eps,eps)))
      }
      if( !igrad ){
        stop("Gradient of log-posterior must be provided.")
      }
      ptm = proc.time()
      fham = foreach( qq = 1:nimp ) %dofuture% {
               q_cal = p_cal$pc_0(p_cal$mpi[qq,], p_cal)
               q_cal$niter = nburn+nmcmc
               q_cal$lpost_0 = lpost_0
               q_cal$glpost_0 = glpost_0
               nuts_imp(xi[qq,],q_cal)
             } %seed% TRUE
      plan(sequential)
      print("Run time:")
      print(proc.time() - ptm)
      cat("\n")
      # extract post-burnin posterior samples
      mpi = NULL
      ipi = seq(1,nmcmc,by=nthin)
      for (qq in 1:nimp){
        tmpi = fham[[qq]]
        tmpi = as.matrix(tmpi[-(1:nburn),])
        # thin posterior samples
        tmpi = tmpi[ipi,,drop=FALSE]
        mpi = rbind(mpi,tmpi)
      }
    }
  }
  if (imcmc == "SMC"){
    # Sequential Monte Carlo (SMC) sampling
    # parallel chains using R package "future"
    require(doFuture)
    # using R code modified from Jason Loeppky
    source(paste(gdir,"/pSMC-ram-mv.r",sep=""),local=TRUE)
    if( !is.null(ncor_smc) ){
      ncor_av = as.numeric(availableCores())
      if( ncor_mc > nimp ){ ncor_mc = nimp }
      if( ncor_smc >= ncor_av ){ ncor_mc = 1; ncor_smc = ncor_av; }
      ncor_tot = ncor_mc*ncor_smc
      if( ncor_tot > ncor_av ){ ncor_mc = ncor_av %/% ncor_smc }
    } else { ncor_smc = 1 }
    if( ncor_mc == 1 ){ pl_1 = "sequential"
    } else { pl_1 = pl }
    if( ncor_smc == 1 ){ pl_2 = "sequential"
    } else { pl_2 = pl }
    plan(list(tweak(pl_1,workers=ncor_mc),tweak(pl_2,workers=ncor_smc)))
    # initial sampling ranges
    if( is.null(lb_smc) ){ slb = rep(-Inf,p_cal$ntheta0)
    } else { slb = lb_smc }
    if( is.null(ub_smc) ){ sub = rep(Inf,p_cal$ntheta0)
    } else { sub = ub_smc }
    islb = is.infinite(slb); isub = is.infinite(sub);
    if( any(islb) || any(isub) ){
      sd_mle = sqrt(diag(p_cal$Sigma_mle$II_nev_it))
      if( any(islb) ){
        jslb = which(islb)
        slb[jslb] = p_cal$mle[jslb] - 10*sd_mle[jslb]
      }
      if( any(isub) ){
        jsub = which(isub)
        sub[jsub] = p_cal$mle[jsub] + 10*sd_mle[jsub]
      }
    }
    ptm = proc.time()
    fsmc = foreach( qq = 1:nimp ) %dofuture% {
             if( nimp == 1 ){ q_cal = p_cal
             } else { q_cal = p_cal$pc_0(p_cal$mpi[qq,], p_cal) }
             SMC(N=nmcmc,M=100,nuseq_T=1,range=cbind(slb,sub),
                 pc=q_cal)
           } %seed% TRUE
    plan(sequential)
    print("Run time:")
    print(proc.time() - ptm)
    cat("\n")
    # acceptance rate
    print("ACCEPTANCE RATES:")
    cat("\n")
    for( qq in 1:nimp ){
      print(paste("Imputation ",qq,": ",fsmc[[qq]]$acr,sep=""))
    }
    cat("\n")
    # extract post-burnin posterior samples
    mpi = NULL
    ipi = seq(1,nmcmc,by=nthin)
    for (qq in 1:nimp){
      tmpi = fsmc[[qq]]$sample
      # thin posterior samples
      tmpi = tmpi[ipi,,drop=FALSE]
      mpi = rbind(mpi,tmpi)
    }
  }

  # Print posterior summary statistics
  print("POSTERIOR SUMMARY")
  cat("\n")
  # Indicator of prior distribution parameters
  p_cal$iPrior = TRUE
  p_cal = p_cal$print_ss_0(mpi,p_cal,
                           levels=c(0.025,0.05,0.5,0.95,0.975))

  #
  # END BAYESIAN ANALYSIS
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

lpo_opt = function(xst,bfgs,pc,tc)
{
  maxiter = 5
  ii = 1
  while( ii <= maxiter ){
    if( bfgs ){
      if( pc$opt_B ){
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$lpost_0,
                gr=pc$glpost_0, pc=pc, method="L-BFGS-B",
                lower=pc$tlb, upper=pc$tub,
                control=list(fnscale = -1,maxit=1000)))
      } else {
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$lpost_0,
                gr=pc$glpost_0, pc=pc, method="BFGS",
                control=list(fnscale = -1,maxit=1000)))
      }
    } else {
      if( pc$opt_B ){
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$lpost_0,
                pc=pc, lower=pc$tlb, upper=pc$tub,
                control=list(fnscale = -1,maxit=10000)))
      } else {
        Catch = pc$tryCatch.W.E(optim(xst, fn=pc$lpost_0,
                pc=pc, control=list(fnscale = -1,maxit=10000)))
      }
    }
    if( is(Catch$value,"list") ){
      return(Catch)
    } else {
      xst = pc$calc_xst(pc,tc)
      while( is.infinite(pc$lpost_0(xst, pc=pc)) ){
        xst = pc$calc_xst(pc,tc)
      }
      ii = ii+1
      if( ii > maxiter ){ return(Catch) } else { next }
    }
  }
}

nuts_imp = function(xi,p_cal)
{
  NUTS(xi,f=p_cal$lpost_0,grad_f=p_cal$glpost_0,n_iter=p_cal$niter,
       delta=0.8,verbose=FALSE)
}
