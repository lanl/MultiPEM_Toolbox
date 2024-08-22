########################################################################
#                                                                      #
# This file contains pre-processing code for reading new event data,   #
# and updating the environment that stores all objects required for    #
# characterization calculations with information from estimation of    #
# the forward and error model parameters using calibration data, and   #
# from the new event.                                                  #
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

prepro_0 = function(p_cal,gdir,adir,rdir,ndir,tnames,nimp=1,bopt=FALSE,
                    itr=FALSE,fp_tr=NULL,tlb=NULL,tub=NULL,tsub=NULL)
{
  #
  # FUNCTION INPUTS
  #

  # p_cal: environment storing all objects needed in characterization
  #        calculations
  # gdir: directory for general subroutines
  # adir: directory for application subroutines
  # rdir: root directory for location of data
  # ndir: new event data directories
  # tnames: names of new event inference parameters
  # nimp: number of calibration parameter imputations utilized in
  #       Markov chain Monte Carlo (MCMC) for new event parameters
  # bopt: bounds supplied to MLE optimization (TRUE/FALSE)
  # itr: new event parameter transform (TRUE/FALSE)
  # fp_tr: fixed parameters for new event parameter transform
  # tlb: lower bounds for new event parameters
  # tub: upper bounds for new event parameters
  # tsub: list containing index sets identifying new event
  #       parameters by phenomenology if full set not present

  #
  # END FUNCTION INPUTS
  #

  #
  # READ IN NEW EVENT DATASETS
  #

  # Set up new event data list
  new_data = vector("list",p_cal$H)

  # Indicator of boundary optimization
  p_cal$opt_B = bopt

  # flag for new event
  p_cal$nev = TRUE
  # names of new event inference parameters
  p_cal$theta_names = tnames
  # number of new event inference parameters
  p_cal$ntheta0 = length(tnames)
  # nev parameter transform
  p_cal$itransform = itr
  if( itr ){
    source(paste(adir,"/transform.r",sep=""),local=TRUE)
    p_cal$tau = tau
    p_cal$j_tau = j_tau
    p_cal$log_absdet_j_tau = log_absdet_j_tau
    p_cal$dlog_absdet_j_tau = dlog_absdet_j_tau
    p_cal$inv_tau = inv_tau
    if( !is.null(fp_tr) ){
      for( na in names(fp_tr) ){
        eval(parse(text=paste("p_cal$tpars$",na," = ","fp_tr$",na,
                              sep="")))
      }
    }
  }
  # nev parameter bounds
  if( is.null(tlb) ){ tlb = rep(-Inf,p_cal$ntheta0) }
  if( is.null(tub) ){ tub = rep(Inf,p_cal$ntheta0) }
  p_cal$theta0_bounds = cbind(tlb,tub)
  if( any(is.finite(tlb)) || any(is.finite(tub)) ){
    itheta0_bounds = vector("list",3)
    itheta0_bounds[[1]] = which(is.finite(tlb) & is.infinite(tub))
    itheta0_bounds[[2]] = which(is.finite(tub) & is.infinite(tlb))
    itheta0_bounds[[3]] = which(is.finite(tlb) & is.finite(tub))
    if( bopt ){
      t_cal = new.env(hash=TRUE)
      t_cal$theta0_bounds = p_cal$theta0_bounds
      t_cal$itheta0_bounds = itheta0_bounds
    } else {
      p_cal$itheta0_bounds = itheta0_bounds
    }
    if( length(itheta0_bounds[[3]]) > 0 ){
      p_cal$theta0_range = tub[itheta0_bounds[[3]]] -
                           tlb[itheta0_bounds[[3]]]
      if( bopt ){ t_cal$theta0_range = p_cal$theta0_range }
      p_cal$sum_theta0_logrange = sum(log(p_cal$theta0_range))
    }
    if( bopt ){
      t_cal$notExp = p_cal$notExp
      t_cal$notLog = p_cal$notLog
    }
  }
  # Read new event data
  for( hh in 1:p_cal$H ){
    tmp_new = read.csv(paste(rdir,"/",ndir[hh],sep=""))
    if( "Source" %in% names(tmp_new) ){
      tmp_new$Source = factor(tmp_new$Source)
    }
    if( "Path" %in% names(tmp_new) ){
      tmp_new$Path = factor(tmp_new$Path)
    }
    if( "Type" %in% names(tmp_new) ){
      tmp_new$Type = factor(tmp_new$Type)
    }
    new_data[[hh]] = tmp_new
  }

  # Fill response and covariate matrices
  # count number of sources by phenomenology
  source_levels = p_cal$source_levels
  nsource = p_cal$ncsource
  # include new event source name
  for( hh in 1:p_cal$H ){
    if( "Source" %in% names(new_data[[hh]]) ){
      source_levels$h[[hh]] = c(source_levels$h[[hh]],
                                levels(new_data[[hh]]$Source))
    } else {
      source_levels$h[[hh]] = c(source_levels$h[[hh]],
                                length(source_levels$h[[hh]])+
                                (1:nrow(new_data[[hh]])))
    }
    nsource[hh] = length(source_levels$h[[hh]])  
    if( nsource[hh] != p_cal$ncsource[hh]+1 ){
      stop(paste("Only one new event permitted for ",
                 "Phenomenology ",hh,".",sep=""))
    }
  }
  # count outputs by phenomenology, source, path, response
  for( hh in 1:p_cal$H ){
    # total number of sources
    p_cal$h[[hh]]$nsource = nsource[hh]
    # response matrix
    p_cal$h[[hh]]$Y = c(p_cal$h[[hh]]$Y,vector("list",1))
    # covariate matrix
    p_cal$h[[hh]]$X = c(p_cal$h[[hh]]$X,vector("list",1))
    # sample size (hir)
    p_cal$h[[hh]]$n = c(p_cal$h[[hh]]$n,vector("list",1))
    # covariance pairs
    p_cal$h[[hh]]$i = c(p_cal$h[[hh]]$i,vector("list",1))
    # path count
    p_cal$h[[hh]]$nplev = rbind(p_cal$h[[hh]]$nplev,
                                rep(0,p_cal$h[[hh]]$Rh))
    # sample size (hijr)
    p_cal$nh[[hh]]$i = c(p_cal$nh[[hh]]$i,vector("list",1))
    # include new event source information
    if( "Source" %in% names(new_data[[hh]]) ){
      isource = (new_data[[hh]]$Source ==
                 source_levels$h[[hh]][nsource[hh]])
    } else { isource = 1 }
    p_cal$h[[hh]]$Y[[nsource[hh]]] = vector("list",p_cal$h[[hh]]$Rh)
    p_cal$h[[hh]]$X[[nsource[hh]]] = vector("list",p_cal$h[[hh]]$Rh)
    p_cal$h[[hh]]$n[[nsource[hh]]] = numeric(p_cal$h[[hh]]$Rh)
    Y = new_data[[hh]][isource,1:p_cal$h[[hh]]$Rh,drop=FALSE]
    if( ncol(new_data[[hh]]) > p_cal$h[[hh]]$Rh ){
      X = new_data[[hh]][isource,-(1:p_cal$h[[hh]]$Rh),drop=FALSE]
    } else { X = NULL }
    for( rr in 1:p_cal$h[[hh]]$Rh ){
      iresp_rr = which(!is.na(Y[,rr]))
      nir = length(iresp_rr)
      if( nir > 0 ){
        p_cal$h[[hh]]$Y[[nsource[hh]]][[rr]] = Y[iresp_rr,rr]
        if( !is.null(X) ){
          p_cal$h[[hh]]$X[[nsource[hh]]][[rr]] =
            X[iresp_rr,,drop=FALSE]
        } else {
          p_cal$h[[hh]]$X[[nsource[hh]]][[rr]] = NULL
        }
        p_cal$h[[hh]]$n[[nsource[hh]]][rr] = nir
      } else { p_cal$h[[hh]]$n[[nsource[hh]]][rr] = 0 }
    }
    p_cal$h[[hh]]$i[[nsource[hh]]]$cov_pairs =
      vector("list",p_cal$h[[hh]]$Rh)
    for( r1 in 1:p_cal$h[[hh]]$Rh ){
      iresp_r1 = !is.na(Y[,r1])
      nr1 = sum(iresp_r1)
      if( nr1 > 0 ){ iresp_r1[which(iresp_r1)] = 1:nr1 }
      p_cal$h[[hh]]$i[[nsource[hh]]]$cov_pairs[[r1]] =
        vector("list",p_cal$h[[hh]]$Rh)
      for( r2 in r1:p_cal$h[[hh]]$Rh ){
        iresp_r2 = !is.na(Y[,r2])
        nr2 = sum(iresp_r2)
        if( nr2 > 0 ){ iresp_r2[which(iresp_r2)] = 1:nr2 }
        if( nr1 > 0 && nr2 > 0 ){
          iresp_r1_r2 = (iresp_r1 & iresp_r2)
          n_r1_r2 = sum(iresp_r1_r2)
          if( n_r1_r2 > 0 ){
            iresp_r1 = iresp_r1[which(iresp_r1_r2)]
            iresp_r2 = iresp_r2[which(iresp_r1_r2)]
            p_cal$h[[hh]]$i[[nsource[hh]]]$cov_pairs[[r1]][[r2]] =
              (iresp_r2-1)*nr1+iresp_r1
          }
        }
      }
    }
    p_cal$h[[hh]]$nev = logical(nsource[hh])
    p_cal$h[[hh]]$nev[nsource[hh]] = TRUE
    # number of paths for new event
    p_cal$nh[[hh]]$i[[nsource[hh]]]$r = vector("list",p_cal$h[[hh]]$Rh)
    for( rr in 1:p_cal$h[[hh]]$Rh ){
      if( p_cal$h[[hh]]$n[[nsource[hh]]][rr] > 0 ){
        if( !is.null(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]) &&
            "Path" %in% names(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]) ){
          lpath =
            levels(factor(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]$Path))
          npath = length(lpath)
          p_cal$h[[hh]]$nplev[nsource[hh],rr] = npath
          p_cal$nh[[hh]]$i[[nsource[hh]]]$r[[rr]] = numeric(npath)
          for( ss in 1:npath ){
            ipath = (p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]$Path ==
                     lpath[ss])
            p_cal$nh[[hh]]$i[[nsource[hh]]]$r[[rr]][ss] = sum(ipath)
          }
        }
      }
    }
  }

  #
  # END READ IN NEW EVENT DATASETS
  #

  #
  # USER SPECIFIED FIELDS
  #

  # Specification of new event parameters by phenomenology
  for( hh in 1:p_cal$H ){
    p_cal$h[[hh]]$theta_names = vector("list",nsource[hh])
  }
  # subset theta0 if necessary
  if( !is.null(tsub) ){
    for( hh in 1:p_cal$H ){
      if( !is.null(tsub[[hh]]) ){
        p_cal$h[[hh]]$itheta0 = tsub[[hh]]
        p_cal$h[[hh]]$theta_names[[nsource[hh]]] =
                            p_cal$theta_names[tsub[[hh]]]
      } else {
        p_cal$h[[hh]]$theta_names[[nsource[hh]]] = p_cal$theta_names
      }
    }
  } else {
    for( hh in 1:p_cal$H ){
      p_cal$h[[hh]]$theta_names[[nsource[hh]]] = p_cal$theta_names
    }
  }

  # Variance component coefficient matrices
  if( p_cal$pvc_1 > 0 ){
    if( p_cal$izmat ){
      # place code calc_zmat_0.r in application directory
      source(paste(adir,"/calc_zmat_0.r",sep=""),local=TRUE)
    } else {
      source(paste(gdir,"/calc_zmat_0.r",sep=""),local=TRUE)
    }
    p_cal = calc_zmat_0(p_cal)
  }

  #
  # END USER SPECIFIED FIELDS
  #

  #
  # ADDITIONAL QUANTITIES USED IN CALCULATIONS
  #

  # Collect sources by emplacement condition
  for( hh in 1:p_cal$H ){
    if( "Th" %in% names(p_cal$h[[hh]]) && p_cal$h[[hh]]$Th > 1 ){
      for( tt in 1:p_cal$h[[hh]]$Th ){
        for( rr in 1:p_cal$h[[hh]]$Rh ){
          if( p_cal$h[[hh]]$n[[nsource[hh]]][rr] > 0 ){
            if( !is.null(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]) &&
                as.numeric(as.character(
                p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]$Type[1])) == tt ){
              p_cal$h[[hh]]$sourceht[[tt]] =
                c(p_cal$h[[hh]]$sourceht[[tt]],nsource[hh])
            }
            break
          }
        }
      }
    }
  }

  # Remove known new event inference parameter values from covariate
  # matrix if present
  for( hh in 1:p_cal$H ){
    if( "theta_names" %in% names(p_cal$h[[hh]]) ){
      ltn = length(p_cal$h[[hh]]$theta_names[[nsource[hh]]])
      for( qq in 1:ltn ){
        for( rr in 1:p_cal$h[[hh]]$Rh ){
          if( p_cal$h[[hh]]$n[[nsource[hh]]][rr] > 0 ){
            if( !is.null(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]) ){
              icol =
                which(colnames(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]])
                  %in% p_cal$h[[hh]]$theta_names[[nsource[hh]]][qq])
              if( length(icol) > 0 ){
                if( ncol(p_cal$h[[hh]]$X[[nsource[hh]]][[rr]]) > 1 ){
                  p_cal$h[[hh]]$X[[nsource[hh]]][[rr]] =
                    p_cal$h[[hh]]$X[[nsource[hh]]][[rr]][,-icol,
                                                         drop=FALSE]
                } else {
                  p_cal$h[[hh]]$X[[nsource[hh]]][[rr]] = NULL
                }
              }
            }
          }
        }
      }
    }
  }

  # Set errors-in-variables yield for new event to FALSE
  if( exists("eiv",where=p_cal,inherits=FALSE) && p_cal$eiv ){
    for( hh in p_cal$ieiv ){
      p_cal$h[[hh]]$eiv = c(p_cal$h[[hh]]$eiv, vector("list",1))
    }
  }

  # Update p_cal for fast computations
  source(paste(gdir,"/make_par_0.r",sep=""),local=TRUE)
  p_cal$pc_0 = pc_0
  p_cal = pc_0(p_cal$mle_cal, p_cal)
  if( exists("mpi",where=p_cal,inherits=FALSE) ){
    if( nimp > 1 ){
      mpi = p_cal$mpi
      nmpi = nrow(mpi)
      p_cal$mpi = NULL
      if( nmpi < nimp ){
        print(paste("The amount of posterior samples available is ",
                    "insufficient for the requested number of ",
                    "imputation samples.",sep=""))
        cat("\n")
        print(paste("Reducing number of imputation samples to match ",
                    "number of available posterior samples.",sep=""))
        cat("\n")
        nimp = nmpi
      }
      if( nmpi > 1 ){
        imp_samp = round(seq(1,nmpi,length=nimp))
        p_cal$mpi = mpi[imp_samp,]
      } else {
        print(paste("Insufficient posterior samples for multiple ",
                    "imputation.",sep=""))
        cat("\n")
      }
    } else { p_cal$mpi = NULL }
  }
  p_cal$nimp = nimp

  #
  # END ADDITIONAL QUANTITIES USED IN CALCULATIONS
  #

  if( bopt ){ return(list(p_cal=p_cal,t_cal=t_cal))
  } else { return(p_cal) }
}
