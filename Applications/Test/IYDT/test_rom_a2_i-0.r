########################################################################
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

require(numDeriv) # numerical derivatives package
set.seed(3001) # set random number seed

ntest <- 4 # number of test inputs

# Acoustic
ncov <- 3 # number of known covariates
X <- matrix(0,ntest,ncov) # known covariate matrix
X[,1] <- log(runif(ntest,min=1,max=10)) # log range
X[,2] <- log(runif(ntest,min=0.5,max=1.5)) # log scaled temp
X[,3] <- log(runif(ntest,min=0.3,max=1.7)) # log scaled pressure 
colnames(X) <- c("lRange","logTempSc","logPressureSc")
               # known covariate names
psim <- list() # list collecting info for calculations
psim$theta_names <- c("W","HOB")
psim$iresp <- TRUE # indicator of impulse
psim$yield_scaling <- 1/3 # yield scaling coefficient
psim$pressure_scaling <- 1/3 # pressure scaling coefficient
psim$temp_scaling <- 1/2 # temperature scaling coefficient
psim$X <- X # known covariates

nsim <- 5 # number of test betas
pbeta <- 3 # number of statistical coefficients
for(ii in 1:nsim){
  psim$beta <- runif(pbeta,min=-10,max=10)
  zeta <- c(log(runif(1,min=0.5,max=10^6)),
            runif(1,min=-10,max=100)) # unknown variables
  fsim <- f0_a(zeta, psim) # call acoustic forward model
  gsim <- g0_a(zeta, psim) # call acoustic jacobian(s)
  jnum <- jacobian(f0_a, zeta, method.args=list(r=6),
                   params=psim) # numerical jacobian
  print("zeta")
  print(zeta)
  print("calculated forward model")
  print(fsim)
  print("calculated jacobian(s)")
  print(gsim)
  print("numerical jacobian")
  print(jnum)
  print("min/max calc-num jacobian")
  print(range(gsim - jnum))
  cat("\n\n")
}
