########################################################################
#                                                                      #
# This file contains code for organizing and printing the parameter    #
# estimates provided to the function. If a parameter matrix is         #
# provided along with quantile levels, the parameter quantiles are     #
# computed in addition to the parameter mean. All results are printed. #
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

print_ss = function(xfin, pc, ci=NULL, levels=NULL)
{
  # function to compute variances and correlations from
  # elements of Cholesky factor
  obs_err_cov = function(x,Rh)
  {           
    L_h = diag(exp(x[1:Rh]),nrow=Rh) 
    x = x[-(1:Rh)]
    if( Rh > 1 ){ L_h[upper.tri(L_h)] = x }
    Sigma_h = t(L_h) %*% L_h
    varp = diag(Sigma_h)
    if( Rh > 1 ){
      Sd_h = diag(1/sqrt(varp))
      Corr_h = Sd_h %*% Sigma_h %*% Sd_h
      return(c(varp,Corr_h[upper.tri(Corr_h)]))
    } else { return(varp) }
  }

  # objects in pc
  pnames = names(pc)

  # If xfin is a matrix and percentiles (levels) are provided,
  # extract them prior to collapsing xfin to the column mean
  # vector
  if( !is.vector(xfin) ){
    sfin = xfin
    if( !is.null(levels) ){ nlevels = length(levels) }
  }

  # print new event inference parameters
  if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
    print("NEW EVENT INFERENCE PARAMETERS")
    cat("\n")
    if( !is.vector(xfin) ){
      stheta0 = as.matrix(sfin[,1:pc$ntheta0])
      colnames(stheta0) = pc$theta_names
      if( exists("itransform",where=pc,inherits=FALSE) ){
        if( pc$itransform ){
          for( ii in 1:nrow(stheta0) ){
            stheta0[ii,] = pc$tau(stheta0[ii,],pc=pc)
          }
        }
      }
      if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
        ith0_bds = pc$itheta0_bounds
        if( length(ith0_bds[[1]]) > 0 ){
          for( jj in ith0_bds[[1]] ){
            stheta0[,jj] = pc$theta0_bounds[jj,1] +
                           pc$notExp(stheta0[,jj])
          }
        }
        if( length(ith0_bds[[2]]) > 0 ){
          for( jj in ith0_bds[[2]] ){
            stheta0[,jj] = pc$theta0_bounds[jj,2] -
                           pc$notExp(stheta0[,jj])
          }
        }
        if( length(ith0_bds[[3]]) > 0 ){
          kk = 0
          for( jj in ith0_bds[[3]] ){
            tau = pc$notExp(stheta0[,jj])
            kk = kk + 1
            stheta0[,jj] = pc$theta0_bounds[jj,1] +
                           pc$theta0_range[kk]*tau/(1+tau)
          }
        }
      }
      pc$tmpi_0 = stheta0
      print(paste("POSTERIOR MEAN: ",round(apply(stheta0,2,mean),2),
            sep=""))
      cat("\n")
      print(paste("POSTERIOR SD: ",round(apply(stheta0,2,sd),2),
            sep=""))
      cat("\n")
      if( !is.null(levels) ){
        qfin = apply(stheta0,2,quantile,probs=levels)
        for( qq in 1:nlevels ){
          print(paste("LEVEL ",100*levels[qq],"%: ",
                      round(qfin[qq,],2),sep=""))
          cat("\n")
        }
      }
      print("CORRELATION MATRIX:")
      cat("\n")
      print(round(cor(stheta0),2))
      cat("\n")
      sfin = as.matrix(sfin[,-(1:pc$ntheta0)])
    } else {
      theta0_it = xfin[1:pc$ntheta0]
      names(theta0_it) = pc$theta_names
      theta0 = theta0_it
      iit = FALSE
      if( exists("itransform",where=pc,inherits=FALSE) ){
        if( pc$itransform ){
          iit = TRUE
          theta0 = pc$tau(theta0, pc=pc)
        }
      }
      if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
        iit = TRUE
        theta0 = pc$transform(theta0, pc=pc)
      }
      if( !pc$iPrior ){ pc$tmle_0 = theta0
      } else { pc$tmap_0 = theta0 }
      print("ESTIMATE: ")
      cat("\n")
      print(round(theta0,2))
      cat("\n")
      if( pc$Sigma_mle_0$acov_0 > 0 && !pc$iPrior ){
        if( pc$Sigma_mle_0$acov_0 >= 1 ){
          print("STANDARD DEVIATION FIXED MODEL PARAMETERS: ")
          cat("\n")
          theta0_sd_0 = sqrt(diag(pc$Sigma_mle_0$II_nev_0))
          names(theta0_sd_0) = pc$theta_names
          print(round(theta0_sd_0,2))
          cat("\n")
          print("CORRELATION MATRIX FIXED MODEL PARAMETERS: ")
          cat("\n")
          if( pc$ntheta0 > 1 ){ ISD_nev_0 = diag(1/theta0_sd_0)
          } else { ISD_nev_0 = 1/theta0_sd_0 }
          C_nev_0 = ISD_nev_0 %*% pc$Sigma_mle_0$II_nev_0 %*% ISD_nev_0
          rownames(C_nev_0) = pc$theta_names
          colnames(C_nev_0) = pc$theta_names
          print(round(as.matrix(C_nev_0),2))
          cat("\n")
        }
        if( pc$Sigma_mle_0$acov_0 == 2 ){
          print("STANDARD DEVIATION: ")
          cat("\n")
          theta0_sd = sqrt(diag(pc$Sigma_mle_0$II_nev))
          names(theta0_sd) = pc$theta_names
          print(round(theta0_sd,2))
          cat("\n")
          print("CORRELATION MATRIX: ")
          cat("\n")
          if( pc$ntheta0 > 1 ){ ISD_nev = diag(1/theta0_sd)
          } else { ISD_nev = 1/theta0_sd }
          C_nev = ISD_nev %*% pc$Sigma_mle_0$II_nev %*% ISD_nev
          rownames(C_nev) = pc$theta_names
          colnames(C_nev) = pc$theta_names
          print(round(as.matrix(C_nev),2))
          cat("\n")
        }
        if( !is.null(ci) ){
          for( qq in 1:length(ci) ){
            z_alpha = qnorm((1-ci[qq])/2,lower.tail=FALSE)
            if( pc$Sigma_mle_0$acov_0 >= 1 ){
              if( iit ){
                lb_0 = theta0_it -
                       z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev_0_it))
                ub_0 = theta0_it +
                       z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev_0_it))
                if( exists("itransform",where=pc,inherits=FALSE) ){
                  if( pc$itransform ){
                    bmat = rbind(lb_0,ub_0)
                    if( pc$ntheta0 > 1 ){
                      bcall = "expand.grid("
                      for( rr in 1:(pc$ntheta0-1) ){
                        bcall = paste(bcall,"bmat[,",rr,"],",sep="")
                      }
                      bcall = paste(bcall,"bmat[,pc$ntheta0])",sep="")
                      bmat = eval(parse(text=bcall))
                    }
                    tbmat = NULL
                    for( rr in 1:nrow(bmat) ){
                      tbmat = rbind(tbmat,pc$tau(bmat[rr,],pc=pc))
                    }
                    lb_0 = apply(tbmat,2,min)
                    ub_0 = apply(tbmat,2,max)
                  }
                }
                if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
                  lb_0 = pc$transform(lb_0, pc=pc)
                  ub_0 = pc$transform(ub_0, pc=pc)
                }
              } else {
                lb_0 = theta0 -
                       z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev_0))
                ub_0 = theta0 +
                       z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev_0))
              }
              print(paste(100*ci[qq],"%: ",
                          "CONFIDENCE INTERVAL FIXED MODEL PARAMETERS:",
                          sep=""))
              cat("\n")
              ci_mat_0 = rbind(lb_0,ub_0)
              colnames(ci_mat_0) = pc$theta_names
              print(round(ci_mat_0,2))
              cat("\n")
            }
            if( pc$Sigma_mle_0$acov_0 == 2 ){
              if( iit ){
                lb = theta0_it -
                     z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev_it))
                ub = theta0_it +
                     z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev_it))
                if( exists("itransform",where=pc,inherits=FALSE) ){
                  if( pc$itransform ){
                    bmat = rbind(lb,ub)
                    if( pc$ntheta0 > 1 ){
                      bcall = "expand.grid("
                      for( rr in 1:(pc$ntheta0-1) ){
                        bcall = paste(bcall,"bmat[,",rr,"],",sep="")
                      }
                      bcall = paste(bcall,"bmat[,pc$ntheta0])",sep="")
                      bmat = eval(parse(text=bcall))
                    }
                    tbmat = NULL
                    for( rr in 1:nrow(bmat) ){
                      tbmat = rbind(tbmat,pc$tau(bmat[rr,],pc=pc))
                    }
                    lb = apply(tbmat,2,min)
                    ub = apply(tbmat,2,max)
                  }
                }
                if( exists("itheta0_bounds",where=pc,inherits=FALSE) ){
                  lb = pc$transform(lb, pc=pc)
                  ub = pc$transform(ub, pc=pc)
                }
              } else {
                lb = theta0 - z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev))
                ub = theta0 + z_alpha*sqrt(diag(pc$Sigma_mle_0$II_nev))
              }
              print(paste(100*ci[qq],"%: ","CONFIDENCE INTERVAL:",
                          sep=""))
              cat("\n")
              ci_mat = rbind(lb,ub)
              colnames(ci_mat) = pc$theta_names
              print(round(ci_mat,2))
              cat("\n")
            }
          }
        }
      }
      xfin = xfin[-(1:pc$ntheta0)]
      cat("\n")
    }
  }

  # print calibration inference parameters
  if( pc$ncalp > 0 ){
    if( exists("nev",where=pc,inherits=FALSE) && pc$nev ){
      Sigma_mle = pc$Sigma_mle_0
    } else {
      Sigma_mle = pc$Sigma_mle_cal
    }
    print("CALIBRATION INFERENCE PARAMETERS")
    cat("\n")
    if( !is.vector(xfin) ){
      scalp = as.matrix(sfin[,1:pc$ncalp])
      colnames(scalp) = pc$cal_par_names
      pc$mpi_calp = scalp
      print(paste("POSTERIOR MEAN: ",round(apply(scalp,2,mean),2),
            sep=""))
      cat("\n")
      print(paste("POSTERIOR SD: ",round(apply(scalp,2,sd),2),
            sep=""))
      cat("\n")
      if( !is.null(levels) ){
        qfin = apply(scalp,2,quantile,probs=levels)
        for( qq in 1:nlevels ){
          print(paste("LEVEL ",100*levels[qq],"%: ",
                      round(qfin[qq,],2),sep=""))
          cat("\n")
        }
      }
      print("CORRELATION MATRIX:")
      cat("\n")
      print(round(cor(scalp),2))
      cat("\n")
      sfin = as.matrix(sfin[,-(1:pc$ncalp)])
    } else {
      calp = xfin[1:pc$ncalp]
      names(calp) = pc$cal_par_names
      if( !pc$iPrior ){ pc$mle_calp = calp
      } else { pc$map_calp = calp }
      print("ESTIMATE: ")
      cat("\n")
      print(round(calp,2))
      cat("\n")
      if( Sigma_mle$acov_cal && !pc$iPrior ){
        print("STANDARD DEVIATION: ")
        cat("\n")
        calp_sd = sqrt(diag(Sigma_mle$II_calp))
        names(calp_sd) = pc$cal_par_names
        print(round(calp_sd,2))
        cat("\n")
        print("CORRELATION MATRIX: ")
        cat("\n")
        if( pc$ncalp > 1 ){ ISD_calp = diag(1/calp_sd)
        } else { ISD_calp = 1/calp_sd }
        C_calp = ISD_calp %*% Sigma_mle$II_calp %*% ISD_calp
        rownames(C_calp) = pc$cal_par_names
        colnames(C_calp) = pc$cal_par_names
        print(round(as.matrix(C_calp),2))
        cat("\n")
        if( !is.null(ci) ){
          for( qq in 1:length(ci) ){
            z_alpha = qnorm((1-ci[qq])/2,lower.tail=FALSE)
            lb = calp - z_alpha*sqrt(diag(Sigma_mle$II_calp))
            ub = calp + z_alpha*sqrt(diag(Sigma_mle$II_calp))
            print(paste(100*ci[qq],"%: ","CONFIDENCE INTERVAL:",sep=""))
            cat("\n")
            ci_mat = rbind(lb,ub)
            colnames(ci_mat) = pc$cal_par_names
            print(round(ci_mat,2))
            cat("\n")
          }
        }
      }
      xfin = xfin[-(1:pc$ncalp)]
      cat("\n")
    }
  }

  # print errors-in-variables inference parameters
  if( pc$nsource > 0 ){
    print("ERRORS-IN-VARIABLES YIELDS")
    cat("\n")
    if( !is.null(levels) ){
      sfin_eiv = sfin[,1:pc$nsource]
      colnames(sfin_eiv) = pc$eiv_sources
      print(paste("POSTERIOR MEAN: ",round(apply(as.matrix(
            sfin_eiv),2,mean),2),sep=""))
      cat("\n")
      qfin = apply(as.matrix(sfin_eiv),2,quantile,
                   probs=levels)
      for( qq in 1:nlevels ){
        print(paste("LEVEL ",100*levels[qq],"%: ",
                    round(qfin[qq,],2),sep=""))
        cat("\n")
      }
      sfin = as.matrix(sfin[,-(1:pc$nsource)])
    } else {
      xfin_eiv = xfin[1:pc$nsource]
      names(xfin_eiv) = pc$eiv_sources
      print(round(xfin_eiv,2))
      xfin = xfin[-(1:pc$nsource)]
      cat("\n")
    }
  }

  # print common coefficients by response
  if( pc$pbeta > 0 ){
    print("COMMON COEFFICIENTS")
    cat("\n")
    for( hh in 1:pc$H ){
      pnamesh = names(pc$h[[hh]])
      pbeta = sum(pc$h[[hh]]$pbeta)
      if( pbeta > 0 ){
        if( !is.null(levels) ){
          qbeta = as.matrix(sfin[,1:pbeta])
        } else {
          beta = xfin[1:pbeta]
        }
        for( rr in 1:pc$h[[hh]]$Rh ){
          if( pc$h[[hh]]$pbeta[rr] > 0 ){
            cat("\t")
            print(paste("Phenomenology: ",hh,
                        "; Response: ",rr,sep=""))
            cat("\n")
            if( !is.null(levels) ){
              cat("\t")
              qbetar = qbeta[,1:pc$h[[hh]]$pbeta[rr]]
              if( "Phen" %in% pnamesh &&
                  pc$h[[hh]]$Phen == "Optical" ){
                qbetar[,3] = notExp(qbetar[,3])
                qbetar[,4] = notExp(qbetar[,4])
              }
              print(paste("POSTERIOR MEAN: ",round(apply(as.matrix(
                    qbetar),2,mean),2),sep=""))
              cat("\n")
              qfin = apply(as.matrix(qbetar),
                           2,quantile,probs=levels)
              for( qq in 1:nlevels ){
                cat("\t")
                print(paste("LEVEL ",100*levels[qq],"%: ",
                            round(qfin[qq,],2),sep=""))
                cat("\n")
              }
              qbeta = as.matrix(qbeta[,-(1:pc$h[[hh]]$pbeta[rr])])
            } else {
              cat("\t")
              betar = beta[1:pc$h[[hh]]$pbeta[rr]]
              if( "Phen" %in% pnamesh &&
                  pc$h[[hh]]$Phen == "Optical" ){
                betar[3] = notExp(betar[3])
                betar[4] = notExp(betar[4])
              }
              print(round(betar,2))
              beta = beta[-(1:pc$h[[hh]]$pbeta[rr])]
              cat("\n")
            }
          }
        }
        if( !is.null(levels) ){
          sfin = as.matrix(sfin[,-(1:pbeta)])
        } else {
          xfin = xfin[-(1:pbeta)]
        }
      }
    }
  }

  # print emplacement condition dependent coefficients by response
  if( pc$ptbeta > 0 ){
    print("EMPLACEMENT CONDITION DEPENDENT COEFFICIENTS")
    cat("\n")
    for( hh in 1:pc$H ){
      pnamesh = names(pc$h[[hh]])
      ptbeta = sum(pc$h[[hh]]$ptbeta)
      if( ptbeta > 0 ){
        if( !is.null(levels) ){
          qbetat = as.matrix(sfin[,1:ptbeta])
        } else {
          betat = xfin[1:ptbeta]
        }
        for( tt in 1:pc$h[[hh]]$Th ){
          for( rr in 1:pc$h[[hh]]$Rh ){
            if( pc$h[[hh]]$pbetat[[tt]][rr] > 0 ){
              cat("\t")
              print(paste("Phenomenology: ",hh,"; Emplacement: ",tt,
                          "; Response: ",rr,sep=""))
              cat("\n")
              if( !is.null(levels) ){
                cat("\t")
                qbetatr = qbetat[,1:pc$h[[hh]]$pbetat[[tt]][rr]]
                if( "Phen" %in% pnamesh &&
                    pc$h[[hh]]$Phen == "Seismic" ){
                  qbetatr[,3] = -notExp(qbetatr[,3])
                }
                print(paste("POSTERIOR MEAN: ",round(apply(as.matrix(
                            qbetatr),2,mean),2),sep=""))
                cat("\n")
                qfin = apply(as.matrix(qbetatr),
                             2,quantile,probs=levels)
                for( qq in 1:nlevels ){
                  cat("\t")
                  print(paste("LEVEL ",100*levels[qq],"%: ",
                              round(qfin[qq,],2),sep=""))
                  cat("\n")
                }
                qbetat = as.matrix(qbetat[,
                                   -(1:pc$h[[hh]]$pbetat[[tt]][rr])])
              } else {
                cat("\t")
                betatr = betat[1:pc$h[[hh]]$pbetat[[tt]][rr]]
                if( "Phen" %in% pnamesh &&
                    pc$h[[hh]]$Phen == "Seismic" ){
                  betatr[3] = -notExp(betatr[3])
                }
                print(round(betatr,2))
                betat = betat[-(1:pc$h[[hh]]$pbetat[[tt]][rr])]
                cat("\n")
              }
            }
          }
        }
        if( !is.null(levels) ){
          sfin = as.matrix(sfin[,-(1:ptbeta)])
        } else {
          xfin = xfin[-(1:ptbeta)]
        }
      }
    }
  }

  # print variance components by response
  if( pc$pvc_1 > 0 ){
    print("SOURCE VARIANCE COMPONENTS")
    cat("\n")
    for( hh in 1:pc$H ){
      if( any(pc$h[[hh]]$pvc_1 > 0) ){
        for( rr in 1:pc$h[[hh]]$Rh ){
          if( pc$h[[hh]]$pvc_1[rr] > 0 ){
            cat("\t")
            print(paste("Phenomenology: ",hh,"; Response: ",rr,sep=""))
            cat("\n")
            if( !is.null(levels) ){
              cat("\t")
              print(paste("POSTERIOR MEAN: ",round(apply(as.matrix(
                    exp(sfin[,1:pc$h[[hh]]$pvc_1[rr]])),2,mean),4),
                    sep=""))
              cat("\n")
              qfin = apply(as.matrix(exp(sfin[,1:pc$h[[hh]]$pvc_1[rr]])),
                           2,quantile,probs=levels)
              for( qq in 1:nlevels ){
                cat("\t")
                print(paste("LEVEL ",100*levels[qq],"%: ",
                            round(qfin[qq,],4),sep=""))
                cat("\n")
              }
              sfin = as.matrix(sfin[,-(1:pc$h[[hh]]$pvc_1[rr])])
            } else {
              cat("\t")
              print(round(exp(xfin[1:pc$h[[hh]]$pvc_1[rr]]),4))
              xfin = xfin[-(1:pc$h[[hh]]$pvc_1[rr])]
              cat("\n")
            }
          }
        }
      }
    }
  }
  if( pc$pvc_2 > 0 ){
    print("PATH VARIANCE COMPONENTS")
    cat("\n")
    for( hh in 1:pc$H ){
      if( any(pc$h[[hh]]$pvc_2 > 0) ){
        for( rr in 1:pc$h[[hh]]$Rh ){
          if( pc$h[[hh]]$pvc_2[rr] > 0 ){
            cat("\t")
            print(paste("Phenomenology: ",hh,"; Response: ",rr,sep=""))
            cat("\n")
            if( !is.null(levels) ){
              cat("\t")
              print(paste("POSTERIOR MEAN: ",round(apply(as.matrix(
                    exp(sfin[,1:pc$h[[hh]]$pvc_2[rr]])),2,mean),4),
                    sep=""))
              cat("\n")
              qfin = apply(as.matrix(exp(sfin[,1:pc$h[[hh]]$pvc_2[rr]])),
                           2,quantile,probs=levels)
              for( qq in 1:nlevels ){
                cat("\t")
                print(paste("LEVEL ",100*levels[qq],"%: ",
                            round(qfin[qq,],4),sep=""))
                cat("\n")
              }
              sfin = as.matrix(sfin[,-(1:pc$h[[hh]]$pvc_2[rr])])
            } else {
              cat("\t")
              print(round(exp(xfin[1:pc$h[[hh]]$pvc_2[rr]]),4))
              xfin = xfin[-(1:pc$h[[hh]]$pvc_2[rr])]
              cat("\n")
            }
          }
        }
      }
    }
  }

  # print observational error covariance parameters
  print("OBSERVATIONAL ERROR COVARIANCE PARAMETERS")
  cat("\n")
  for( hh in 1:pc$H ){
    print(paste("Phenomenology ",hh,sep=""))
    cat("\n")
    ncpar = pc$h[[hh]]$Rh*(pc$h[[hh]]$Rh+1)/2
    if( !is.null(levels) ){
      qfin = NULL
      for( qq in 1:nrow(sfin) ){
        cfin = obs_err_cov(as.matrix(sfin[qq,1:ncpar]),pc$h[[hh]]$Rh)
        qfin = rbind(qfin,cfin)
      }
      qfin_mu = apply(as.matrix(qfin),2,mean)
      print("POSTERIOR MEAN:")
      cat("\n")
      print("Variances")
      cat("\n")
      print(round(qfin_mu[1:pc$h[[hh]]$Rh],4))
      cat("\n")
      if( pc$h[[hh]]$Rh > 1 ){
        print("Correlations")
        cat("\n")
        Corr_h = diag(1,pc$h[[hh]]$Rh)
        Corr_h[upper.tri(Corr_h)] = qfin_mu[-(1:pc$h[[hh]]$Rh)]
        print(round(Corr_h,2))
        cat("\n")
      }
      qfin = apply(as.matrix(qfin),2,quantile,probs=levels)
      sfin = as.matrix(sfin[,-(1:ncpar)])
      for( qq in 1:nlevels ){
        print("Variances")
        cat("\n")
        cat("\t")
        print(paste("LEVEL ",100*levels[qq],"%: ",
                    round(qfin[qq,1:pc$h[[hh]]$Rh],4),sep=""))
        cat("\n")
        if( pc$h[[hh]]$Rh > 1 ){
          print("Correlations")
          cat("\n")
          cat("\t")
          Corr_h = diag(1,pc$h[[hh]]$Rh)
          Corr_h[upper.tri(Corr_h)] = qfin[qq,-(1:pc$h[[hh]]$Rh)]
          print(paste("LEVEL ",100*levels[qq],"%:",sep=""))
          cat("\t")
          print(round(Corr_h,2))
          cat("\n")
        }
      }
    } else {
      cfin = obs_err_cov(xfin[1:ncpar],pc$h[[hh]]$Rh)
      xfin = xfin[-(1:ncpar)]
      print("Variances")
      cat("\n")
      print(round(cfin[1:pc$h[[hh]]$Rh],4))
      cat("\n")
      if( pc$h[[hh]]$Rh > 1 ){
        print("Correlations")
        cat("\n")
        Corr_h = diag(1,pc$h[[hh]]$Rh)
        Corr_h[upper.tri(Corr_h)] = cfin[-(1:pc$h[[hh]]$Rh)]
        print(round(Corr_h,2))
        cat("\n")
      }
    }
  }

  if( pc$nsource > 0 && pc$iPrior ){
    # print FGSN prior parameters
    print("FGSN PRIOR PARAMETERS")
    cat("\n")
    if( !is.null(levels) ){
      cfin = cbind(sfin[,1],exp(sfin[,2]),sfin[,3:pc$p_fgsn])
      qfin_mu = apply(as.matrix(cfin),2,mean)
      print("POSTERIOR MEAN:")
      cat("\n")
      print(paste("Alpha = ",round(qfin_mu[1],2),sep=""))
      print(paste("Lambda squared = ",round(qfin_mu[2],2),sep=""))
      print(paste("Omega = ",round(qfin_mu[-(1:2)],2),sep=""))
      cat("\n")
      qfin = apply(as.matrix(cfin),2,quantile,probs=levels)
      print("Alpha:")
      for( qq in 1:nlevels ){
        print(paste("LEVEL ",100*levels[qq],"%: ",
                    round(qfin[qq,1],2),sep=""))
        cat("\n")
      }
      print("Lambda squared:")
      for( qq in 1:nlevels ){
        print(paste("LEVEL ",100*levels[qq],"%: ",
                    round(qfin[qq,2],2),sep=""))
        cat("\n")
      }
      print("Omega:")
      for( qq in 1:nlevels ){
        print(paste("LEVEL ",100*levels[qq],"%: ",
                    round(qfin[qq,-(1:2)],2),sep=""))
        cat("\n")
      }
      sfin = as.matrix(sfin[,-(1:pc$p_fgsn)])
    } else {
      print(paste("Alpha = ",round(xfin[1],2),sep=""))
      print(paste("Lambda squared = ",round(exp(xfin[2]),2),sep=""))
      print(paste("Omega = ",round(xfin[3:pc$p_fgsn],2),sep=""))
      xfin = xfin[-(1:pc$p_fgsn)]
      cat("\n")
    }
  }

  if( pc$p_A > 0 && pc$iPrior ){
    # print scale parameters for variance component priors
    print("VARIANCE COMPONENT PRIOR SCALE PARAMETERS")
    cat("\n")
    ipvc_1 = FALSE; ipvc_2 = FALSE;
    if( pc$pvc_1 > 0 ){ ipvc_1 = TRUE }
    if( pc$pvc_2 > 0 ){ ipvc_2 = TRUE }
    for( hh in 1:pc$H ){
      if( (ipvc_1 && any(pc$h[[hh]]$pvc_1 > 0)) ||
          (ipvc_2 && any(pc$h[[hh]]$pvc_2 > 0)) ){
        for( rr in 1:pc$h[[hh]]$Rh ){
          if( (ipvc_1 && pc$h[[hh]]$pvc_1[rr] > 0) ||
              (ipvc_2 && pc$h[[hh]]$pvc_2[rr] > 0) ){
            cat("\t")
            print(paste("Phenomenology: ",hh,"; Response: ",rr,
                  sep=""))
            cat("\n")
            if( !is.null(levels) ){
              cat("\t")
              print(paste("POSTERIOR MEAN: ",
                          round(mean(exp(sfin[,1])),2),sep=""))
              cat("\n")
              qfin = quantile(exp(sfin[,1]),probs=levels)
              for( qq in 1:nlevels ){
                cat("\t")
                print(paste("LEVEL ",100*levels[qq],"%: ",
                      round(qfin[qq],2),sep=""))
                cat("\n")
              }
              sfin = as.matrix(sfin[,-1])
            } else {
              cat("\t")
              print(round(exp(xfin[1]),2))
              xfin = xfin[-1]
              cat("\n")
            }
          }
        }
      }
    }
  }
  return(pc)
}
