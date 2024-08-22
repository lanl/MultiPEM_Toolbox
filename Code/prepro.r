########################################################################
#                                                                      #
# This file contains pre-processing code for reading calibration and   #
# (optionally) new event data, and setting up the environment that     #
# stores all objects required for characterization calculations.       #
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

prepro = function(gdir,adir,rdir,cdir,Rh,pbeta,bopt=FALSE,nev=FALSE,
                  itr=FALSE,izmat=FALSE,ieiv=NULL,seiv=NULL,ewsd=NULL,
                  Th=NULL,pbetat=NULL,ibetar=NULL,pvc_1=NULL,pvc_2=NULL,
                  tnames=NULL,cnames=NULL,fp_tr=NULL,tlb=NULL,tub=NULL,
                  ndir=NULL,tsub=NULL)
{
  #
  # FUNCTION INPUTS
  #

  # gdir: directory for general subroutines
  # adir: directory for application subroutines
  # rdir: root directory for location of data
  # cdir: directory locations for calibration data
  # Rh: number of responses for each phenomenology
  # pbeta: list containing empirical model common parameter counts
  #        by phenomenology
  # bopt: bounds supplied to MLE optimization (TRUE/FALSE)
  # nev: analysis of new event (TRUE/FALSE)
  # itr: new event parameter transform (TRUE/FALSE)
  # izmat: user-provided code for computing variance component
  #        coefficient matrices (TRUE/FALSE)
  # ieiv: numerical identifier of phenomenologies utilizing
  #       errors-in-variables yields
  # seiv: list containing identifiers of sources containing errors-
  #       in-variables yields by phenomenology ("ALL" - every source)
  # ewsd: standard deviation of errors-in-variables Gaussian
  #       likelihood
  # Th: number of emplacement conditions for each phenomenology
  # pbetat: list containing empirical model emplacement-dependent
  #         parameter counts by phenomenology
  # ibetar: list containing locations of empirical model common
  #         parameters in full parameter vector by phenomenology
  # pvc_1: list containing level 1 variance component parameter
  #        counts by phenomenology
  # pvc_2: list containing level 2 variance component parameter
  #        counts by phenomenology
  # tnames: names of new event inference parameters
  # cnames: names of calibration inference parameters
  # fp_tr: fixed parameters for new event parameter transform
  # tlb: lower bounds for new event parameters
  # tub: upper bounds for new event parameters
  # ndir: new event data directories
  # tsub: list containing index sets identifying new event
  #       parameters by phenomenology if full set not present

  #
  # END FUNCTION INPUTS
  #

  #
  # READ IN ALL DATASETS
  #

  # Determine number of phenomenologies
  H = length(cdir)

  # Set up calibration data list
  cal_data = vector("list",H)
  if( nev ){
    # Set up new event data list
    new_data = vector("list",H)
  }

  # Set up fixed parameter list for function calls
  p_cal = new.env(hash=TRUE)
  p_cal$h = vector("list",H)

  # Indicator of boundary optimization
  p_cal$opt_B = bopt

  # names of calibration inference parameters
  p_cal$cal_par_names = cnames
  # number of calibration inference parameters
  p_cal$ncalp = length(cnames)

  # Collect unique source levels by phenomenology
  source_levels = vector("list",0)
  source_levels$h = vector("list",H)

  # Read calibration data
  p_cal$H = H
  for( hh in 1:H ){
    tmp_cal = read.csv(paste(rdir,"/",cdir[hh],sep=""))
    if( "Source" %in% names(tmp_cal) ){
      tmp_cal$Source = factor(tmp_cal$Source)
    }
    if( "Path" %in% names(tmp_cal) ){
      tmp_cal$Path = factor(tmp_cal$Path)
    }
    if( "Type" %in% names(tmp_cal) ){
      tmp_cal$Type = factor(tmp_cal$Type)
    }
    cal_data[[hh]] = tmp_cal
  }

  # Source transform code
  source(paste(gdir,"/transform.r",sep=""),local=TRUE)
  p_cal$transform = transform
  p_cal$notExp = notExp
  p_cal$dnotExp = dnotExp
  p_cal$d2notExp = d2notExp
  p_cal$inv_transform = inv_transform
  p_cal$notLog = notLog

  # Code for inclusion of new event data in combined analysis
  if( nev ){
    p_cal$nev = TRUE
    # names of new event inference parameters
    if( is.null(tnames) ){
      stop("New event inference parameter names must be provided.")
    }
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
    if( is.null(ndir) ){
      stop("New event data directories must be provided.")
    }
    for( hh in 1:H ){
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
  } else { p_cal$ntheta0 = 0 }

  # Fill response and covariate matrices
  # count number of sources by phenomenology
  nsource = ncsource = numeric(H)
  for( hh in 1:H ){
    # source levels
    if( "Source" %in% names(cal_data[[hh]]) ){
      source_levels$h[[hh]] = levels(cal_data[[hh]]$Source)
    } else {
      source_levels$h[[hh]] = 1:nrow(cal_data[[hh]])
    }
    nsource[hh] = ncsource[hh] = length(source_levels$h[[hh]])
  }
  if( nev ){
    # include new event source name
    for( hh in 1:H ){
      if( "Source" %in% names(new_data[[hh]]) ){
        source_levels$h[[hh]] = c(source_levels$h[[hh]],
                                  levels(new_data[[hh]]$Source))
      } else {
        source_levels$h[[hh]] = c(source_levels$h[[hh]],
                                  length(source_levels$h[[hh]])+
                                  (1:nrow(new_data[[hh]])))
      }
      nsource[hh] = length(source_levels$h[[hh]])  
      if( nsource[hh] != ncsource[hh]+1 ){
        stop(paste("Only one new event permitted for ",
                   "Phenomenology ",hh,".",sep=""))
      }
    }
  }
  # count outputs by phenomenology, source, path, response
  p_cal$nh = vector("list",H)
  for( hh in 1:H ){
    # total number of sources
    p_cal$h[[hh]]$nsource = nsource[hh]
    # response matrix
    p_cal$h[[hh]]$Y = vector("list",nsource[hh])
    # covariate matrix
    p_cal$h[[hh]]$X = vector("list",nsource[hh])
    # sample size (hir)
    p_cal$h[[hh]]$n = vector("list",nsource[hh])
    # covariance pairs
    p_cal$h[[hh]]$i = vector("list",nsource[hh])
    # response count
    p_cal$h[[hh]]$Rh = Rh[hh]
    # path count
    p_cal$h[[hh]]$nplev = matrix(0,nsource[hh],Rh[hh])
    # sample size (hijr)
    p_cal$nh[[hh]]$i = vector("list",nsource[hh])
    for( qq in 1:ncsource[hh] ){
      if( "Source" %in% names(cal_data[[hh]]) ){
        isource = (cal_data[[hh]]$Source == source_levels$h[[hh]][qq])
      } else { isource = qq }
      p_cal$h[[hh]]$Y[[qq]] = vector("list",Rh[hh])
      p_cal$h[[hh]]$X[[qq]] = vector("list",Rh[hh])
      p_cal$h[[hh]]$n[[qq]] = numeric(Rh[hh])
      Y = cal_data[[hh]][isource,1:Rh[hh],drop=FALSE]
      if( ncol(cal_data[[hh]]) > Rh[hh] ){
        X = cal_data[[hh]][isource,-(1:Rh[hh]),drop=FALSE]
      } else { X = NULL }
      for( rr in 1:Rh[hh] ){
        iresp_rr = which(!is.na(Y[,rr]))
        nir = length(iresp_rr)
        if( nir > 0 ){
          p_cal$h[[hh]]$Y[[qq]][[rr]] = Y[iresp_rr,rr]
          if( !is.null(X) ){
            p_cal$h[[hh]]$X[[qq]][[rr]] = X[iresp_rr,,drop=FALSE]
          } else {
            p_cal$h[[hh]]$X[[qq]][[rr]] = NULL
          }
          p_cal$h[[hh]]$n[[qq]][rr] = nir
        } else { p_cal$h[[hh]]$n[[qq]][rr] = 0 }
      }
      p_cal$h[[hh]]$i[[qq]]$cov_pairs = vector("list",Rh[hh])
      for( r1 in 1:Rh[hh] ){
        iresp_r1 = !is.na(Y[,r1])
        nr1 = sum(iresp_r1)
        if( nr1 > 0 ){ iresp_r1[which(iresp_r1)] = 1:nr1 }
        p_cal$h[[hh]]$i[[qq]]$cov_pairs[[r1]] = vector("list",Rh[hh])
        for( r2 in r1:Rh[hh] ){
          iresp_r2 = !is.na(Y[,r2])
          nr2 = sum(iresp_r2)
          if( nr2 > 0 ){ iresp_r2[which(iresp_r2)] = 1:nr2 }
          if( nr1 > 0 && nr2 > 0 ){
            iresp_r1_r2 = (iresp_r1 & iresp_r2)
            n_r1_r2 = sum(iresp_r1_r2)
            if( n_r1_r2 > 0 ){
              iresp_r1 = iresp_r1[which(iresp_r1_r2)]
              iresp_r2 = iresp_r2[which(iresp_r1_r2)]
              p_cal$h[[hh]]$i[[qq]]$cov_pairs[[r1]][[r2]] =
                (iresp_r2-1)*nr1+iresp_r1
            }
          }
        }
      }
      # number of paths per source
      p_cal$nh[[hh]]$i[[qq]]$r = vector("list",Rh[hh])
      for( rr in 1:Rh[hh] ){
        if( p_cal$h[[hh]]$n[[qq]][rr] > 0 ){
          if( !is.null(p_cal$h[[hh]]$X[[qq]][[rr]]) &&
              "Path" %in% names(p_cal$h[[hh]]$X[[qq]][[rr]]) ){
            lpath = levels(factor(p_cal$h[[hh]]$X[[qq]][[rr]]$Path))
            npath = length(lpath)
            p_cal$h[[hh]]$nplev[qq,rr] = npath
            p_cal$nh[[hh]]$i[[qq]]$r[[rr]] = numeric(npath)
            for( ss in 1:npath ){
              ipath = (p_cal$h[[hh]]$X[[qq]][[rr]]$Path == lpath[ss])
              p_cal$nh[[hh]]$i[[qq]]$r[[rr]][ss] = sum(ipath)
            }
          }
        }
      }
    }
    if( nev ){
      # include new event source information
      if( "Source" %in% names(new_data[[hh]]) ){
        isource = (new_data[[hh]]$Source ==
                   source_levels$h[[hh]][nsource[hh]])
      } else { isource = 1 }
      p_cal$h[[hh]]$Y[[nsource[hh]]] = vector("list",Rh[hh])
      p_cal$h[[hh]]$X[[nsource[hh]]] = vector("list",Rh[hh])
      p_cal$h[[hh]]$n[[nsource[hh]]] = numeric(Rh[hh])
      Y = new_data[[hh]][isource,1:Rh[hh],drop=FALSE]
      if( ncol(new_data[[hh]]) > Rh[hh] ){
        X = new_data[[hh]][isource,-(1:Rh[hh]),drop=FALSE]
      } else { X = NULL }
      for( rr in 1:Rh[hh] ){
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
      p_cal$h[[hh]]$i[[nsource[hh]]]$cov_pairs = vector("list",Rh[hh])
      for( r1 in 1:Rh[hh] ){
        iresp_r1 = !is.na(Y[,r1])
        nr1 = sum(iresp_r1)
        if( nr1 > 0 ){ iresp_r1[which(iresp_r1)] = 1:nr1 }
        p_cal$h[[hh]]$i[[nsource[hh]]]$cov_pairs[[r1]] =
          vector("list",Rh[hh])
        for( r2 in r1:Rh[hh] ){
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
      p_cal$nh[[hh]]$i[[nsource[hh]]]$r = vector("list",Rh[hh])
      for( rr in 1:Rh[hh] ){
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
    for( rr in 1:Rh[hh] ){
      iresp_rr = which(!is.na(cal_data[[hh]][,rr]))
      lsource = length(unique(cal_data[[hh]]$Source[iresp_rr]))
      vcFlag = FALSE
      if( lsource == 1 ){
        vcFlag = TRUE
        print(paste("Warning: Insufficient Sources ",
                    "for Variance Component models with ",
                    "Phenomenology ",hh," and Response ",rr,".",
                    sep=""))
      } else {
        tsource = table(cal_data[[hh]]$Source[iresp_rr])
        if( all(tsource <= 1) ){
          vcFlag = TRUE
          print(paste("Warning: Insufficient number of observations ",
                      "per Source for Variance Component models with ",
                      "Phenomenology ",hh," and Response ",rr,".",
                      sep=""))
        }
      }
      if( !vcFlag ){
        if( all(p_cal$h[[hh]]$nplev[1:ncsource[hh],rr] <= 1) ){
                                      print(paste(
                                      "Warning: Insufficient Paths for",
                                      " Level 2 Variance Component",
                                      " models with Phenomenology ",hh,
                                      " and Response ",rr,".",sep=""))
        } else {
          if( "Path" %in% names(cal_data[[hh]]) ){
            tsource = names(tsource)
            npobs = NULL
            for( qq in tsource ){
              tpath = cal_data[[hh]]$Path[iresp_rr]
              psource = which(cal_data[[hh]]$Source[iresp_rr] == qq)
              ppath = table(tpath[psource])
              if( length(ppath) > 1 ){ npobs = c(npobs,max(ppath)) }
            }
            if( all(npobs <= 1) ){
              print(paste("Warning: Insufficient number of ",
                          "observations per Path for Level 2 ",
                          "Variance Component models with ",
                          "Phenomenology ",hh,"and Response ",
                          rr,".",sep=""))
            }
          }
        }
      }
    }
  }

  #
  # END READ IN ALL DATASETS
  #

  #
  # USER SPECIFIED FIELDS
  #

  # Specification of new event (nev) parameters by phenomenology
  if( nev ){
    for( hh in 1:H ){
      p_cal$h[[hh]]$theta_names = vector("list",nsource[hh])
    }
    # subset theta0 if necessary
    if( !is.null(tsub) ){
      for( hh in 1:H ){
        if( !is.null(tsub[[hh]]) ){
          p_cal$h[[hh]]$itheta0 = tsub[[hh]]
          p_cal$h[[hh]]$theta_names[[nsource[hh]]] =
                              p_cal$theta_names[tsub[[hh]]]
        } else {
          p_cal$h[[hh]]$theta_names[[nsource[hh]]] = p_cal$theta_names
        }
      }
    } else {
      for( hh in 1:H ){
        p_cal$h[[hh]]$theta_names[[nsource[hh]]] = p_cal$theta_names
      }
    }
  }

  # Specification of errors-in-variables (eiv) yields by phenomenology
  if( !is.null(ieiv) ){
    if( is.null(seiv) ){
      stop(paste("Sources with errors-in-variables yields must be",
                 " provided for relevant phenomenologies.",sep=""))
    }
    p_cal$eiv = TRUE
    eiv_sources = vector("list",0)
    eiv_sources$h = vector("list",H)
    for( hh in ieiv ){ eiv_sources$h[[hh]] = seiv[[hh]] }
  } else { p_cal$nsource = 0 }

  # Empirical model parameter count: common
  for( hh in 1:H ){ p_cal$h[[hh]]$pbeta = pbeta[[hh]] }

  if( !is.null(Th) ){
    # Empirical model parameter count: emplacement condition
    for( hh in 1:H ){
      if( Th[hh] > 1 ){
        if( is.null(pbetat) ){
          stop(paste("Emplacement condition parameter counts",
                     " must be provided.",sep=""))
        }
        p_cal$h[[hh]]$pbetat = pbetat[[hh]]
      }
    }
  }

  if( !is.null(Th) ){
    # Locations of common and emplacement-dependent parameters in
    # full parameter vector
    if( is.null(ibetar) ){
      ibetar = vector("list",H)
      for( hh in 1:H ){
        if( Th[hh] > 1 ){
          # lists with elements for each response within
          # emplacement condition
          ibetar[[hh]] = vector("list",Th[hh]*Rh[hh])
        }
      }
    }
    for( hh in 1:H ){
      if( Th[hh] > 1 ){
        p_cal$h[[hh]]$ibetar = ibetar[[hh]]
        p_cal$h[[hh]]$ibetatr = vector("list",Th[hh]*Rh[hh])
        for( tt in 1:Th[hh] ){
          if( !is.null(pbetat[[hh]][[tt]]) ){
            for( rr in 1:Rh[hh] ){
              if( pbetat[[hh]][[tt]][rr] > 0 ){
                p_cal$h[[hh]]$ibetatr[[(tt-1)*Rh[hh]+rr]] =
                  setdiff(1:(pbeta[[hh]][rr]+pbetat[[hh]][[tt]][rr]),
                          ibetar[[hh]][[(tt-1)*Rh[hh]+rr]])
              }
            }
          }
        }
      }
    }
  }

  # Level 1 variance component parameter count
  p_cal$pvc_1 = 0
  if( !is.null(pvc_1) ){
    for( hh in 1:H ){ p_cal$h[[hh]]$pvc_1 = pvc_1[[hh]] }
    for( hh in 1:H ){
      p_cal$pvc_1 = p_cal$pvc_1+sum(p_cal$h[[hh]]$pvc_1)
    }
  }

  # Level 2 variance component parameter count
  p_cal$pvc_2 = 0
  if( !is.null(pvc_2) ){
    for( hh in 1:H ){ p_cal$h[[hh]]$pvc_2 = pvc_2[[hh]] }
    for( hh in 1:H ){
      p_cal$pvc_2 = p_cal$pvc_2+sum(p_cal$h[[hh]]$pvc_2)
    }
  }

  # Variance component coefficient matrices
  if( p_cal$pvc_1 > 0 ){
    if( izmat ){
      # place code calc_zmat.r in application directory
      source(paste(adir,"/calc_zmat.r",sep=""),local=TRUE)
    } else {
      source(paste(gdir,"/calc_zmat.r",sep=""),local=TRUE)
    }
    p_cal = calc_zmat(p_cal)
  }

  # subset calibration inference parameters if necessary
  if( !is.null(cnames) ){
    for( hh in 1:p_cal$H ){
      csub = which(p_cal$cal_par_names %in% colnames(cal_data[[hh]]))
      if( length(csub) > 0 ){
        # Specification of calibration parameters by phenomenology
        p_cal$h[[hh]]$cal_par_names = p_cal$cal_par_names[csub]
      }
    }
  }

  #
  # END USER SPECIFIED FIELDS
  #

  #
  # ADDITIONAL QUANTITIES USED IN CALCULATIONS
  #

  if( !is.null(Th) ){
    # Collect sources by emplacement condition
    for( hh in 1:H ){
      if( Th[hh] > 1 ){
        p_cal$h[[hh]]$sourceht = vector("list",Th[hh])
        for( tt in 1:Th[hh] ){
          for( ii in 1:nsource[hh] ){
            for( rr in 1:Rh[hh] ){
              if( p_cal$h[[hh]]$n[[ii]][rr] > 0 ){
                if( !is.null(p_cal$h[[hh]]$X[[ii]][[rr]]) &&
                      as.numeric(as.character(
                      p_cal$h[[hh]]$X[[ii]][[rr]]$Type[1])) == tt ){
                  p_cal$h[[hh]]$sourceht[[tt]] =
                    c(p_cal$h[[hh]]$sourceht[[tt]],ii)
                }
                break
              }
            }
          }
        }
      }
    }
  }

  # Remove known new event inference parameter values from covariate
  # matrix if present (relevant only for "new event" test cases)
  if( nev ){
    for( hh in 1:H ){
      if( "theta_names" %in% names(p_cal$h[[hh]]) ){
        ltn = length(p_cal$h[[hh]]$theta_names[[nsource[hh]]])
        for( qq in 1:ltn ){
          for( rr in 1:Rh[hh] ){
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
  }

  # Setup for errors-in-variables yields
  if( !is.null(ieiv) ){
    # Source lists by phenomenology "hh" for errors-in-variables yields
    # Use "ALL" for every source in phenomenology "hh" data
    eiv_source_list = NULL
    for( hh in ieiv ){
      if( eiv_sources$h[[hh]] == "ALL" ){
        eiv_source_list = c(eiv_source_list,
                            source_levels$h[[hh]][1:ncsource[hh]])
      } else {
        ilist = which(source_levels$h[[hh]][1:ncsource[hh]] %in%
                      eiv_sources$h[[hh]])
        eiv_source_list=c(eiv_source_list,source_levels$h[[hh]][ilist])
      }
    }
    eiv_source_list = unique(eiv_source_list)
    p_cal$eiv_sources = eiv_source_list
    p_cal$nsource = length(eiv_source_list)
    # Set mean of eiv Gaussian likelihood to nominal value given
    # in data
    p_cal$eiv_w = vector("list",p_cal$nsource)
    # Mapping from phenomenology/source to unique source list
    for( hh in ieiv ){
      p_cal$h[[hh]]$eiv = vector("list",nsource[hh])
      for( qq in 1:ncsource[hh] ){
        for( rr in 1:Rh[hh] ){
          if( p_cal$h[[hh]]$n[[qq]][rr] > 0 ){
            if( !is.null(p_cal$h[[hh]]$X[[qq]][[rr]]) ){
              if( "Source" %in% names(p_cal$h[[hh]]$X[[qq]][[rr]]) ){
                source = p_cal$h[[hh]]$X[[qq]][[rr]]$Source[1]
              } else { source = qq }
              if( source %in% eiv_source_list ){
                ilist = which(eiv_source_list %in% source)
                if( is.null(p_cal$h[[hh]]$eiv[[qq]]) ){
                  p_cal$h[[hh]]$eiv[[qq]] = ilist
                }
                if( is.null(p_cal$eiv_w[[ilist]]) ){
                  p_cal$eiv_w[[ilist]] =
                    p_cal$h[[hh]]$X[[qq]][[rr]]$W[1]
                }
              }
              # remove nominal yield from covariate matrix (use eiv
              # yield)
              icol =
                which(colnames(p_cal$h[[hh]]$X[[qq]][[rr]]) %in% "W")
              if( length(icol) > 0 ){
                if( ncol(p_cal$h[[hh]]$X[[qq]][[rr]]) > 1 ){
                  p_cal$h[[hh]]$X[[qq]][[rr]] =
                    p_cal$h[[hh]]$X[[qq]][[rr]][,-icol,drop=FALSE]
                } else {
                  p_cal$h[[hh]]$X[[qq]][[rr]] = NULL
                }
              }
            }
          }
        }
      }
    }
    # Set standard deviation of eiv Gaussian likelihood
    p_cal$eiv_w = unlist(p_cal$eiv_w)
    if( is.null(ewsd) ){
      stop("Errors-in-variables standard deviation must be provided.")
    }
    p_cal$eiv_w_sd = ewsd
  }

  # Total empirical model parameter count: common
  p_cal$pbeta = 0
  for( hh in 1:H ){ p_cal$pbeta = p_cal$pbeta+sum(p_cal$h[[hh]]$pbeta) }

  # Total emplacement-dependent parameter count across responses
  # by emplacement condition
  p_cal$ptbeta = 0
  if( !is.null(Th) ){
    for( hh in 1:H ){
      if( Th[hh] > 1 ){
        p_cal$h[[hh]]$ptbeta = numeric(Th[hh])
        for( tt in 1:Th[hh] ){
          p_cal$h[[hh]]$ptbeta[tt] = sum(p_cal$h[[hh]]$pbetat[[tt]])
        }
        ptbeta = sum(p_cal$h[[hh]]$ptbeta)
        if( ptbeta > 0 ){ p_cal$h[[hh]]$Th = Th[hh] }
        p_cal$ptbeta = p_cal$ptbeta+ptbeta
      }
    }
  }

  # Remove known calibration inference parameter values from covariate
  # matrix if present
  for( hh in 1:p_cal$H ){
    if( "cal_par_names" %in% names(p_cal$h[[hh]]) ){
      ltn = length(p_cal$h[[hh]]$cal_par_names)
      for( ii in 1:ncsource[hh] ){
        for( qq in 1:ltn ){
          for( rr in 1:p_cal$h[[hh]]$Rh ){
            if( p_cal$h[[hh]]$n[[ii]][rr] > 0 ){
              if( !is.null(p_cal$h[[hh]]$X[[ii]][[rr]]) ){
                icol =
                  which(colnames(p_cal$h[[hh]]$X[[ii]][[rr]])
                    %in% p_cal$h[[hh]]$cal_par_names[qq])
                if( length(icol) > 0 ){
                  if( !is.na(p_cal$h[[hh]]$X[[ii]][[rr]][1,icol]) ){
                    if( ncol(p_cal$h[[hh]]$X[[ii]][[rr]]) > 1 ){
                      p_cal$h[[hh]]$X[[ii]][[rr]] =
                        p_cal$h[[hh]]$X[[ii]][[rr]][,-icol,drop=FALSE]
                    } else {
                      p_cal$h[[hh]]$X[[ii]][[rr]] = NULL
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }

  # Total number of model parameters
  p_cal$nmpars = p_cal$ntheta0+p_cal$ncalp+p_cal$nsource+p_cal$pbeta+
                 p_cal$ptbeta+p_cal$pvc_1+p_cal$pvc_2+sum(Rh*(Rh+1))/2

  #
  # END ADDITIONAL QUANTITIES USED IN CALCULATIONS
  #

  if( bopt ){ return(list(p_cal=p_cal,t_cal=t_cal))
  } else { return(p_cal) }
}
