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
set.seed(10001) # set random number seed

pbeta <- 1 # number of statistical coefficients
ntest <- 2 # number of test inputs

# Crater
ncov <- 1 # number of known covariates
X <- matrix(0,ntest,ncov) # known covariate matrix
X[,1] <- log(runif(ntest,min=0.5,max=10^6)) # log yield
colnames(X) <- "W" # known covariate names
psim <- list() # list collecting info for calculations
psim$pbeta <- pbeta # number of statistical coefficients
psim$yield_scaling <- 1/3 # yield scaling coefficient
psim$X <- X # known covariates
psim$cal <- FALSE # indicator of global calibration parameters

nsim <- 5 # number of test betas
for(ii in 1:nsim){
  zeta <- runif(pbeta,min=-10,max=10) # unknown variables
  fsim <- f_c(zeta, psim) # call crater forward model
  gsim <- g_c(zeta, psim) # call crater jacobian(s)
  jnum <- jacobian(f_c, zeta, method.args=list(r=6),
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
