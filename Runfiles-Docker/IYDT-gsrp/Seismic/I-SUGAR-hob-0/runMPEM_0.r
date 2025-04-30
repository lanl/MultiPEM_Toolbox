########################################################################
#                                                                      #
# This file is the input deck for MultiPEM Toolbox rapid post-         #
# detonation analysis, based on using fixed values of the forward and  #
# error model parameters obtained from calibration data.               #
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
# LOAD R DATA FILE FROM
# CALIBRATION INFERENCE
#

# Specify local directory
loc_dir = "I-SUGAR-hob-0"

# Specify directory for reading .RData
# single-stage analysis
#r_dir = loc_dir
# two-stage analysis
r_dir = "seismic_cal"

load(paste(r_dir,"/.RData",sep=""))

#
# END LOAD R DATA FILE FROM
# CALIBRATION INFERENCE
#

#
# REQUIRED R PACKAGES
#

require(Matrix)

#
# END REQUIRED R PACKAGES
#

#
# PREPROCESSING
#

# Specify directory for general subroutines
gen_dir = "Code"

# Source supporting R function
source(paste(gen_dir,"/prepro_0.r",sep=""))

# Specify directory for application subroutines
app_dir = "Code/IYDT-gsrp"

# Specify root data directory
dat_dir = "Data"

# Specify new event data directories
dat_new = "seismic_new.csv"

# Phenomenologies for this analysis
# 1 - seismic

# Names of new event inference parameters
theta_names = c("W","HOB")

# Number of calibration parameter imputations utilized in
# Markov chain Monte Carlo (MCMC) for new event parameters
nimpute = 1000

# Set flag for bounded optimization
opt_B = FALSE

# Indicate nev parameter transform
itransform = FALSE

# Specify fixed parameters for nev parameter transform
if( itransform ){
  tPars = TRUE

  if( tPars ){
    tPars = vector("list",0)
    tPars$yield_scaling = 1/3
  } else { tPars = NULL }
} else { tPars = NULL }

# Set up parameter constraints
# lower and upper bounds (use -Inf and Inf if unbounded)
lb_theta0 = rep(-Inf,length(theta_names))
lb_theta0[2] = 0
ub_theta0 = rep(Inf,length(theta_names))
ub_theta0[2] = 160

# Set up parameter subsets by phenomenology
tsub = FALSE

if( tsub ) {
  tsub = vector("list",length(Rh))
  #tsub[[1]] =
} else { tsub = NULL }

# Preprocessing for statistical analysis routines
tmp = prepro_0(p_cal,gen_dir,app_dir,dat_dir,dat_new,theta_names,
               nimp=nimpute,bopt=opt_B,itr=itransform,fp_tr=tPars,
               tlb=lb_theta0,tub=ub_theta0,tsub=tsub)
if( opt_B ){
  p_cal = tmp$p_cal
  t_cal = tmp$t_cal
} else {
  p_cal = tmp
  t_cal = NULL
}
rm(tmp)
save.image(file=paste(loc_dir,"/.RData",sep=""))

#
# END PREPROCESSING
#

#
# MAXIMUM LIKELIHOOD CALCULATION
#

# Source supporting R function
source(paste(gen_dir,"/calc_mle_0.r",sep=""))

# Set seed for repeatability of analysis
set.seed(231)

# Names of forward models for each response
# by phenomenology
fm0 = vector("list",length(Rh))
fm0[[1]] = c("f0_s","f0_s")

# Indicate if forward model gradients provided
igrad = TRUE

if( igrad ){
  # Names of forward model gradients for each response
  # by phenomenology
  gfm0 = vector("list",length(Rh))
  gfm0[[1]] = c("g0_s","g0_s")
} else { gfm0 = NULL }

# Specify number of starting values for optimization
nstart = 30

# number of cores to use for optimization
ncores_mle = 30

# Indicate use of BFGS optimization methods
bfgs = TRUE

# Location of R data files with starting values
# for input to MLE optimization
opt_files_in = NULL

# Location of R data file to write the results of
# MLE optimization
opt_files_out = paste(loc_dir,"/opt_nev.RData",sep="")

# Initial start value for theta0
tst = FALSE

if( tst ){
  tst = numeric(p_cal$ntheta0)
  #tst[2] =
} else { tst = NULL }

# Confidence interval levels for new event parameter inference
ci_lev = 0.95

# Indicator of MLE gradient check
mle_grad_ck = TRUE

# Strategy for running parallel jobs (future package)
parallel_plan = "multicore"

# MLE calculations
p_cal = calc_mle_0(p_cal,gen_dir,app_dir,fm0,nst=nstart,ncor=ncores_mle,
                   ci_lev=ci_lev,igrad=igrad,bfgs=bfgs,
                   igrck=mle_grad_ck,t_cal=t_cal,g0=gfm0,
                   fopt_in=opt_files_in,Xst=NULL,tst=tst,
                   fopt_out=opt_files_out,pl=parallel_plan)
save.image(file=paste(loc_dir,"/.RData",sep=""))

#
# END MAXIMUM LIKELIHOOD CALCULATION
#

#
# BAYESIAN ANALYSIS
#

# Specify if Bayesian analysis is to be conducted
iBayes = TRUE

if( iBayes ){
  # Source supporting R function
  source(paste(gen_dir,"/calc_bayes_0.r",sep=""))

  # Indicator of prior distribution for theta0
  iTheta0Prior = FALSE

  if( iTheta0Prior ){
    # location of code for computing log-prior densities and gradients
    prior_files_theta0 = NULL
    if( igrad ){
      gr_prior_files_theta0 = NULL
    } else { gr_prior_files_theta0 = NULL }

    # prior distribution for new event parameters (theta0)
    p_cal$lp_theta0$f = "lp_0"
    if( igrad ){ p_cal$lp_theta0$g = "lq_0" }

    # parameters for log yield parameter prior (Gaussian)
    p_cal$pi_w_mu = (log(10)+log(10000000))/2
    p_cal$pi_w_sd = (log(10000000)-log(10))/6
    # parameters for HOB/DOB parameter prior (Gaussian)
    p_cal$pi_h_mu = 0
    p_cal$pi_h_sd = 160/3
  } else {
    prior_files_theta0 = NULL
    gr_prior_files_theta0 = NULL
  }

  # specify MCMC algorithm
  # options: "RAM", "FME", "NUTS", or "SMC"
  iMCMC = "FME"

  # burn-in
  nburn = 1000

  # production
  nmcmc = 4000

  # posterior sample thinning rate (for multiple imputation)
  nthin = 200

  # number of cores to use for optimization
  ncores_map = 30

  # number of cores to use for generating parallel MCMC chains
  ncores_mc = 30

  # Indicator of prior gradient check
  prior_grad_ck = TRUE

  # Indicator of posterior gradient check
  post_grad_ck = TRUE

  # Options for Sequential Monte Carlo (SMC) sampling
  # (iMCMC = "SMC")
  # number of cores to use for parallelization within SMC algorithm
  ncores_smc = NULL
  # new event parameter ranges for initial SMC sample
  lb_smc = rep(-Inf,length(theta_names))
  ub_smc = rep(Inf,length(theta_names))

  # Bayesian calculations
  p_cal = calc_bayes_0(p_cal,gen_dir,app_dir,nst=nstart,nburn=nburn,
                       nmcmc=nmcmc,nthin=nthin,ncor_map=ncores_map,
                       ncor_mc=ncores_mc,igrad=igrad,
                       igrck_pr=prior_grad_ck,igrck_po=post_grad_ck,
                       bfgs=bfgs,itpr=iTheta0Prior,
                       fpr_t=prior_files_theta0,
                       fgpr_t=gr_prior_files_theta0,imcmc=iMCMC,
                       pl=parallel_plan,ncor_smc=ncores_smc,
                       lb_smc=lb_smc,ub_smc=ub_smc,t_cal=t_cal)
  save.image(file=paste(loc_dir,"/.RData",sep=""))
}

#
# END BAYESIAN ANALYSIS
#
