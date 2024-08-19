########################################################################
#                                                                      #
# This file contains code for calculating the (inverse) observed       #
# information matrix of the calibration inference parameters.          #
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

obs_info_cal = function(p_cal,opt,imle=TRUE)
{
  #
  # FUNCTION INPUTS
  #

  # p_cal: environment storing all objects needed in characterization
  #        calculations
  # opt: object resulting from call to optim()
  # imle: flag for opt from likelihood maximization (TRUE/FALSE)

  #
  # END FUNCTION INPUTS
  #

  Cmat = list()
  Cmat$acov_cal = 0

  Hess = -opt$hessian
  Hess = forceSymmetric(Hess)
  cHessCatch = p_cal$tryCatch.W.E(chol(Hess))
  delta = 1.e-17
  idelta = FALSE
  while( !is(cHessCatch$value,"Matrix") ){
    idelta = TRUE
    delta = 10*delta
    dHess = Hess
    diag(dHess) = diag(dHess) + delta
    cHessCatch = p_cal$tryCatch.W.E(chol(dHess))
  }
  cHess = cHessCatch$value
  IHess = chol2inv(cHess)
  if( imle ){
    Cmat$II_calp = IHess[1:p_cal$ncalp, 1:p_cal$ncalp]
    Cmat$acov_cal = 1
  } else {
    p_cal$IHess = IHess
  }
  if( idelta ){
    print(paste("Perturbation added to Hessian diagonals: ",delta,
                sep=""))
  }

  if( imle ){ return(Cmat)
  } else { return(p_cal) }
}
