########################################################################
#                                                                      #
# This file contains code for organizing and printing the new event    #
# parameter estimates provided to the function. If a parameter matrix  #
# is provided along with quantile levels, the parameter quantiles are  #
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

print_ss_0 = function(xfin, pc, ci=NULL, levels=NULL)
{
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
  print("NEW EVENT INFERENCE PARAMETERS")
  cat("\n")
  if( !is.vector(xfin) ){
    stheta0 = as.matrix(sfin)
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
    pc$tmpi = stheta0
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
  } else {
    theta0_it = xfin
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
    if( !pc$iPrior ){ pc$tmle = theta0
    } else { pc$tmap = theta0 }
    print("ESTIMATE: ")
    cat("\n")
    print(round(theta0,2))
    cat("\n")
    if( pc$Sigma_mle$acov && !pc$iPrior ){
      print("STANDARD DEVIATION: ")
      cat("\n")
      theta0_sd = sqrt(diag(pc$Sigma_mle$II_nev))
      names(theta0_sd) = pc$theta_names
      print(round(theta0_sd,2))
      cat("\n")
      if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
        print("STANDARD DEVIATION FIXED MODEL PARAMETERS: ")
        cat("\n")
        theta0_sd_0 = sqrt(diag(pc$Sigma_mle$II_nev_0))
        names(theta0_sd_0) = pc$theta_names
        print(round(theta0_sd_0,2))
        cat("\n")
      }
      print("CORRELATION MATRIX: ")
      cat("\n")
      if( pc$ntheta0 > 1 ){ ISD_nev = diag(1/theta0_sd)
      } else { ISD_nev = 1/theta0_sd }
      C_nev = ISD_nev %*% pc$Sigma_mle$II_nev %*% ISD_nev
      rownames(C_nev) = pc$theta_names
      colnames(C_nev) = pc$theta_names
      print(round(as.matrix(C_nev),2))
      cat("\n")
      if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
        print("CORRELATION MATRIX FIXED MODEL PARAMETERS: ")
        cat("\n")
        if( pc$ntheta0 > 1 ){ ISD_nev_0 = diag(1/theta0_sd_0)
        } else { ISD_nev_0 = 1/theta0_sd_0 }
        C_nev_0 = ISD_nev_0 %*% pc$Sigma_mle$II_nev_0 %*% ISD_nev_0
        rownames(C_nev_0) = pc$theta_names
        colnames(C_nev_0) = pc$theta_names
        print(round(as.matrix(C_nev_0),2))
        cat("\n")
      }
      if( !is.null(ci) ){
        for( qq in 1:length(ci) ){
          z_alpha = qnorm((1-ci[qq])/2,lower.tail=FALSE)
          if( iit ){
            lb = theta0_it - z_alpha*sqrt(diag(pc$Sigma_mle$II_nev_it))
            ub = theta0_it + z_alpha*sqrt(diag(pc$Sigma_mle$II_nev_it))
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
            lb = theta0 - z_alpha*sqrt(diag(pc$Sigma_mle$II_nev))
            ub = theta0 + z_alpha*sqrt(diag(pc$Sigma_mle$II_nev))
          }
          print(paste(100*ci[qq],"%: ","CONFIDENCE INTERVAL:",sep=""))
          cat("\n")
          ci_mat = rbind(lb,ub)
          colnames(ci_mat) = pc$theta_names
          print(round(ci_mat,2))
          cat("\n")
          if( exists("rapid",where=pc,inherits=FALSE) && pc$rapid ){
            if( iit ){
              lb_0 = theta0_it -
                     z_alpha*sqrt(diag(pc$Sigma_mle$II_nev_0_it))
              ub_0 = theta0_it +
                     z_alpha*sqrt(diag(pc$Sigma_mle$II_nev_0_it))
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
              lb_0 = theta0 - z_alpha*sqrt(diag(pc$Sigma_mle$II_nev_0))
              ub_0 = theta0 + z_alpha*sqrt(diag(pc$Sigma_mle$II_nev_0))
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
        }
      }
    }
  }
  return(pc)
}
