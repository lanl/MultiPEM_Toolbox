########################################################################
#                                                                      #
# This file is the input deck for MultiPEM Toolbox estimation of       #
# forward and error model parameters based on calibration data.        #
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
gen_dir = "../../../Code"

# Source supporting R function
source(paste(gen_dir,"/prepro_cal.r",sep=""))

# Specify directory for application subroutines
app_dir = "../../Code"

# Specify root data directory
dat_dir = "../../Data"

# Specify calibration data directories
dat_cal = c("seismic_cal.csv",
            "acoustic_cal.csv",
            "optical_cal.csv",
            "crater_cal.csv")

# Phenomenologies for this analysis
# 1 - seismic
# 2 - acoustic
# 3 - optical
# 4 - crater (surface effects)

# Indicate presence of calibration inference parameters
calp = FALSE

if( calp ){
  # Names of calibration inference parameters
  #cal_par_names =
} else { cal_par_names = NULL }

# Specify number of responses for each phenomenology
Rh = c(2,2,2,2)

# Empirical model parameter count: common
# list with elements corresponding to phenomenologies
pbeta = vector("list",length(Rh))
for( hh in 1:length(Rh) ){ pbeta[[hh]] = numeric(Rh[hh]) }
# phenomenology 2
pbeta[[2]] = c(2,2)
# phenomenology 3
pbeta[[3]] = c(4,4)
# phenomenology 4
pbeta[[4]] = c(2,2)

# Specify number of emplacement conditions for each phenomenology
Th = TRUE

if( Th ){ Th = c(3,3,0,0)
} else { Th = NULL }

# Empirical model parameter count: emplacement condition
# list with elements corresponding to phenomenologies
if( !is.null(Th) ){
  pbetat = vector("list",length(Rh))
  for( hh in 1:length(Rh) ){
    if( Th[hh] > 1 ){ pbetat[[hh]] = vector("list",Th[hh]) }
  }
  # phenomenology 1
  for( tt in 1:Th[1] ){
    pbetat[[1]][[tt]] = numeric(Rh[1])
    pbetat[[1]][[tt]] = c(5,5)
  }
  # phenomenology 2
  for( tt in 1:Th[2] ){
    pbetat[[2]][[tt]] = numeric(Rh[2])
    pbetat[[2]][[tt]] = c(1,1)
  }
} else { pbetat = NULL }

# Locations of common parameters in full parameter vector
# list with elements corresponding to phenomenologies
if( !is.null(Th) ){
  ibetar = vector("list",length(Rh))
  for( hh in 1:length(Rh) ){
    if( Th[hh] > 1 ){
      # lists with elements for each response within
      # emplacement condition
      ibetar[[hh]] = vector("list",Th[hh]*Rh[hh])
    }
  }
  # phenomenology 2
  for( tt in 1:Th[2] ){
    for( rr in 1:Rh[2] ){
      ibetar[[2]][[(tt-1)*Rh[2]+rr]] = 1:2
    }
  }
} else { ibetar = NULL }

# Indicate analysis with errors-in-variables (eiv)
eiv = FALSE

# Specifications for errors-in-variables
if( eiv ){
  # Specify phenomenologies utilizing
  # errors-in-variables yields
  ieiv = 3:4

  # Errors-in-variables source lists by
  # phenomenology
  seiv = vector("list",length(Rh))
  for( hh in ieiv ){ seiv[[hh]] = "ALL" }

  # Set standard deviation of eiv Gaussian likelihood
  eiv_w_sd = 0.3/3
} else {
  ieiv = NULL
  seiv = NULL
  eiv_w_sd = NULL
}

# Specify Error Model
# Source variance component parameter count
pvc_1 = TRUE

if( pvc_1 ){
  pvc_1 = vector("list",length(Rh))
  for( hh in 1:length(Rh) ){ pvc_1[[hh]] = numeric(Rh[hh]) }
  # phenomenology 1
  pvc_1[[1]] = c(1,1)
  # phenomenology 2
  pvc_1[[2]] = c(1,1)
} else { pvc_1 = NULL }

# Path variance component parameter count
pvc_2 = TRUE

if( pvc_2 ){
  pvc_2 = vector("list",length(Rh))
  for( hh in 1:length(Rh) ){ pvc_2[[hh]] = numeric(Rh[hh]) }
  # phenomenology 1
  pvc_2[[1]] = c(1,1)
  # phenomenology 2
  pvc_2[[2]] = c(1,1)

  # path error models by phenomenology
  ptype = vector("list",length(Rh))
  # phenomenology 1
  ptype[[1]] = "Crossed"
  # phenomenology 2
  ptype[[2]] = "Crossed"
} else { pvc_2 = NULL; ptype = NULL; }

# Set flag for user-provided code to calculate variance
# component coefficient matrices
calc_Z = FALSE

# Preprocessing for statistical analysis routines
p_cal = prepro_cal(gen_dir,app_dir,dat_dir,dat_cal,Rh,pbeta,
                   izmat=calc_Z,ieiv=ieiv,seiv=seiv,ewsd=eiv_w_sd,
                   Th=Th,pbetat=pbetat,ibetar=ibetar,pvc_1=pvc_1,
                   pvc_2=pvc_2,ptype=ptype,cnames=cal_par_names)
save.image()

#
# END PREPROCESSING
#

#
# MAXIMUM LIKELIHOOD CALCULATION
#

# Source supporting R function
source(paste(gen_dir,"/calc_mle_cal.r",sep=""))

# Set seed for repeatability of analysis
set.seed(601)

# Names of forward models for each response
# by phenomenology
fm = vector("list",length(Rh))
fm[[1]] = c("f_s","f_s")
fm[[2]] = c("f_a","f_a")
fm[[3]] = c("f_o","f_o")
fm[[4]] = c("f_c","f_c")

# Indicate if forward model gradients provided
igrad = TRUE

if( igrad ){
  # Names of forward model gradients for each response
  # by phenomenology
  gfm = vector("list",length(Rh))
  gfm[[1]] = c("g_s","g_s")
  gfm[[2]] = c("g_a","g_a")
  gfm[[3]] = c("g_o","g_o")
  gfm[[4]] = c("g_c","g_c")
} else { gfm = NULL }

# Specifications for forward model calculations
# a) flags for modified forward model calculation by
#    response for each relevant phenomenology
iResponse = TRUE

if( iResponse ){
  iResponse = vector("list",length(Rh))
  iResponse[[1]] = c(TRUE,FALSE)
  iResponse[[2]] = c(TRUE,FALSE)
} else { iResponse = NULL }

# b) fixed quantities required by forward models
fPars = TRUE

if( fPars ){
  fPars = vector("list",length(Rh))
  fPars[[1]]$yield_scaling = 1/3
  fPars[[2]]$yield_scaling = 1/3
  fPars[[2]]$pressure_scaling = 1/3
  fPars[[2]]$temp_scaling = 1/2
  fPars[[3]]$yield_scaling = 1/3
} else { fPars = NULL }

# Specify number of starting values for optimization
nstart = 30

# number of cores to use for optimization
ncores_mle = 30

# Indicate use of BFGS optimization methods
bfgs = TRUE

# Location of R data files with starting values
# for input to MLE optimization
opt_files_in = c("../Opt/opt_1_0.RData",
                 "../Opt/opt_2_0.RData",
                 "../Opt/opt_3_0.RData",
                 "../Opt/opt_4_0.RData")

# Location of R data file to write the results of
# MLE optimization
opt_files_out = "./opt.RData"

if( calp ){
  # Initial start value for calibration inference parameters
  cst = FALSE

  if( cst ){
    cst = numeric(p_cal$ncalp)
    #cst[1] =
  } else { cst = NULL }

  # Confidence interval levels for calibration parameter inference
  ci_lev = 0.95
} else {
  cst = NULL
  ci_lev = NULL
}

# Indicate phenomenology number and type (if needed
# for postprocessing)
Phen = TRUE

if( Phen ){
  Phen = matrix(c(1,3,"Seismic","Optical"),nrow=2)
} else { Phen = NULL }

# Indicator of MLE gradient check
mle_grad_ck = TRUE

# Strategy for running parallel jobs (future package)
parallel_plan = "multicore"

# MLE calculations
p_cal = calc_mle_cal(p_cal,gen_dir,app_dir,fm,nst=nstart,
                     ncor=ncores_mle,ci_lev=ci_lev,igrad=igrad,
                     bfgs=bfgs,igrck=mle_grad_ck,g=gfm,iresp=iResponse,
                     fp_fm=fPars,fopt_in=opt_files_in,Xst=NULL,
                     cst=cst,fopt_out=opt_files_out,phen=Phen,
                     pl=parallel_plan)
save.image()

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
  source(paste(gen_dir,"/calc_bayes_cal.r",sep=""))

  # Indicator of prior distribution for forward model
  # coefficients
  iBetaPrior = TRUE

  if( iBetaPrior ){
    # location of code for computing log-prior densities and gradients
    prior_files_beta = c("../Code/lp_beta_s.r","../Code/lp_beta_o.r")
    if( igrad ){
      gr_prior_files_beta = c("../Code/glp_beta_s.r",
                              "../Code/glp_beta_o.r")
    } else { gr_prior_files_beta = NULL }

    # prior distribution for phenomenology 1
    # forward model coefficients
    p_cal$h[[1]]$lp_betat = vector("list",Th[1])
    for( tt in 1:Th[1] ){
      p_cal$h[[1]]$lp_betat[[tt]]$f = c("lp_s","lp_s")
      if( igrad ){
        p_cal$h[[1]]$lp_betat[[tt]]$g = c("lq_s","lq_s")
      }
    }

    # prior distribution for phenomenology 3
    # forward model coefficients
    p_cal$h[[3]]$lp_beta$f = c("lp_o","lp_o")
    if( igrad ){
      p_cal$h[[3]]$lp_beta$g = c("lq_o","lq_o")
    }
  } else {
    prior_files_beta = NULL
    gr_prior_files_beta = NULL
  }

  # Indicator of prior distribution for calibration parameters
  iCalPrior = FALSE

  if( calp && iCalPrior ){
    # location of code for computing log-prior densities and gradients
    prior_files_calp = NULL
    if( igrad ){
      gr_prior_files_calp = NULL
    } else { gr_prior_files_calp = NULL }

    # prior distribution for calibration parameters (calp)
    p_cal$lp_calp$f = "lp_c"
    if( igrad ){ p_cal$lp_calp$g = "lq_c" }

    # parameters for calibration parameter prior distribution
    #p_cal$pi_c_mu =
    #p_cal$pi_c_sd =
  } else {
    prior_files_calp = NULL
    gr_prior_files_calp = NULL
  }

  # fixed scale parameters for variance component prior
  # comment out if these parameters should vary
  p_cal$A = 20

  # eta parameter in Lewandowski-Kurowicka-Joe (LKJ) prior
  # distribution for correlation parameters
  p_cal$lp_corr$eta = 1

  # FGSN parameters for errors-in-variables yields prior
  # number of components
  p_cal$K = 0
  # total number of FGSN parameters
  p_cal$p_fgsn = 0
  if( eiv ){
    p_cal$K = 2
    p_cal$p_fgsn = p_cal$K + 2
  }

  # specify Markov chain Monte Carlo (MCMC) algorithm
  # options: "RAM", "FME", or "NUTS"
  iMCMC = "RAM"

  # burn-in
  nburn = 10000

  # production
  nmcmc = 20000

  # posterior sample thinning rate
  nthin = 20

  # number of cores to use for optimization
  ncores_map = 30

  # number of cores to use for generating parallel MCMC chains
  ncores_mc = 30

  # Indicator of prior gradient check
  prior_grad_ck = TRUE

  # Indicator of posterior gradient check
  post_grad_ck = TRUE

  # Bayesian calculations
  p_cal = calc_bayes_cal(p_cal,gen_dir,app_dir,nst=nstart,nburn=nburn,
                         nmcmc=nmcmc,nthin=nthin,ncor_map=ncores_map,
                         ncor_mc=ncores_mc,igrad=igrad,
                         igrck_pr=prior_grad_ck,igrck_po=post_grad_ck,
                         bfgs=bfgs,ibpr=iBetaPrior,icpr=iCalPrior,
                         fpr_b=prior_files_beta,
                         fgpr_b=gr_prior_files_beta,
                         fpr_c=prior_files_calp,
                         fgpr_c=gr_prior_files_calp,
                         Xnom=NULL,imcmc=iMCMC,pl=parallel_plan)
  save.image()
}

#
# END BAYESIAN ANALYSIS
#
