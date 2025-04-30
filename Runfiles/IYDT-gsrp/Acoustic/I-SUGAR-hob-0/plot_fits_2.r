########################################################################
#                                                                      #
# This file plots the fitted vs. the observed data, where the fits are #
# obtained by evaluating the forward model using the coefficient MLEs  #
# and the observed data are adjusted for source effects to highlight   #
# the presence of path bias.                                           #
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

source("../../../Code/predict.r")

pred_mle = predict(c(p_cal$mle,p_cal$mle_cal), p_cal)

npg = 3 # number of path groups
ms_paths = vector("list",npg)
psy = vector("list",npg)
ms_paths[[1]] = vector("list",2)
psy[[1]] = vector("list",2)
ms_paths[[2]] = vector("list",2)
psy[[2]] = vector("list",2)
ms_paths[[3]] = vector("list",2)
psy[[3]] = vector("list",2)

ms_paths[[1]][[1]] = c("HTA-11","HTA-8","HTA-15")
psy[[1]][[1]] = c(2,5,6)
ms_paths[[2]][[1]] = c("HRR-10","HR-95","HRR-33")
psy[[2]][[1]] = c(7,9,10)
ms_paths[[3]][[1]] = c("HTP-6","HTP-16")
psy[[3]][[1]] = c(17,19)

ms_paths[[1]][[2]] = c("HR-103","HR-102","HR-101","HR-108")
psy[[1]][[2]] = c(2,4,5,6)
ms_paths[[2]][[2]] = c("HR-95","HRR-6","HRR-1")
psy[[2]][[2]] = c(7,9,10)
ms_paths[[3]][[2]] = "HTP-12"
psy[[3]][[2]] = 17

pcol = c("gray","pink","lightgreen")
pcol_ms = c("black","red","darkgreen")
for( hh in 1:p_cal$H ){
  obs = vector("list",p_cal$h[[hh]]$Rh)
  fit = vector("list",p_cal$h[[hh]]$Rh)
  for( rr in 1:p_cal$h[[hh]]$Rh ){
    for( ii in 1:p_cal$h[[hh]]$nsource ){
      # adjusted for source effects
      obs[[rr]] = c(obs[[rr]],
                    pred_mle$h[[hh]]$observed[[ii]][[rr]]-
                    pred_mle$h[[hh]]$source_effects[[ii]][[rr]])
      fit[[rr]] = c(fit[[rr]],pred_mle$h[[hh]]$fitted[[ii]][[rr]])
    }
  }
  for( rr in 1:p_cal$h[[hh]]$Rh ){
    pdf(paste("fit_vs_obs_h",hh,"_r",rr,"_2.pdf",sep=""))
    iType = FALSE
    for( ii in 1:p_cal$h[[hh]]$nsource ){
      if( p_cal$h[[hh]]$n[[ii]][rr] > 0 ){
        if( "Type" %in% names(p_cal$h[[hh]]$X[[ii]][[rr]]) ){
          tt = as.numeric(as.character(
               p_cal$h[[hh]]$X[[ii]][[rr]]$Type[1]))
          if( tt > 1 ){ iType = TRUE }
        } else { tt = 1 }
        ipath = vector("list",npg)
        ipath_all = NULL
        for( jj in 1:npg ){
          ipath[[jj]] = vector("list",length(ms_paths[[jj]][[rr]]))
          for( kk in 1:length(ms_paths[[jj]][[rr]]) ){
            ipath[[jj]][[kk]] =
              which(p_cal$h[[hh]]$X[[ii]][[rr]]$Path %in%
                    ms_paths[[jj]][[rr]][kk])
            ipath_all = c(ipath_all,ipath[[jj]][[kk]])
          }
        }
        ipath_all = sort(ipath_all)
        jpath = setdiff(1:p_cal$h[[hh]]$n[[ii]][[rr]],ipath_all)
        for( jj in 1:npg ){
          for( kk in 1:length(ms_paths[[jj]][[rr]]) ){
            if( ii == 1 && jj == 1 && kk == 1 ){
              if( length(jpath) > 0 ){
                plot(pred_mle$h[[hh]]$observed[[ii]][[rr]][jpath] -
                     pred_mle$h[[hh]]$source_effects[[ii]][[rr]][jpath],
                     pred_mle$h[[hh]]$fitted[[ii]][[rr]][jpath],
                     xlim=range(obs[[rr]]),ylim=range(fit[[rr]]),
                     xlab="observed",ylab="fitted",
                     pch=20,col=pcol[tt])
                if( length(ipath[[jj]][[kk]]) > 0 ){
                  points(pred_mle$h[[hh]]$observed[[ii]][[rr]][
                         ipath[[jj]][[kk]]] -
                         pred_mle$h[[hh]]$source_effects[[ii]][[rr]][
                         ipath[[jj]][[kk]]],
                         pred_mle$h[[hh]]$fitted[[ii]][[rr]][
                         ipath[[jj]][[kk]]],
                         pch=psy[[jj]][[rr]][kk],col=pcol_ms[jj])
                }
              } else if( length(ipath[[jj]][[kk]]) > 0 ){
                plot(pred_mle$h[[hh]]$observed[[ii]][[rr]][
                     ipath[[jj]][[kk]]] -
                     pred_mle$h[[hh]]$source_effects[[ii]][[rr]][
                     ipath[[jj]][[kk]]],
                     pred_mle$h[[hh]]$fitted[[ii]][[rr]][
                     ipath[[jj]][[kk]]],
                     xlim=range(obs[[rr]]),ylim=range(fit[[rr]]),
                     xlab="observed",ylab="fitted",
                     pch=psy[[jj]][[rr]][kk],col=pcol_ms[jj])
              }
            } else {
              if( length(jpath) > 0 ){
                points(pred_mle$h[[hh]]$observed[[ii]][[rr]][jpath] -
                       pred_mle$h[[hh]]$source_effects[[ii]][[rr]][
                       jpath],
                       pred_mle$h[[hh]]$fitted[[ii]][[rr]][jpath],
                       pch=20,col=pcol[tt])
                if( length(ipath[[jj]][[kk]]) > 0 ){
                  points(pred_mle$h[[hh]]$observed[[ii]][[rr]][
                         ipath[[jj]][[kk]]] -
                         pred_mle$h[[hh]]$source_effects[[ii]][[rr]][
                         ipath[[jj]][[kk]]],
                         pred_mle$h[[hh]]$fitted[[ii]][[rr]][
                         ipath[[jj]][[kk]]],
                         pch=psy[[jj]][[rr]][kk],col=pcol_ms[jj])
                }
              } else if( length(ipath[[jj]][[kk]]) > 0 ){
                points(pred_mle$h[[hh]]$observed[[ii]][[rr]][ipath[[jj]][[kk]]] -
                       pred_mle$h[[hh]]$source_effects[[ii]][[rr]][ipath[[jj]][[kk]]],
                       pred_mle$h[[hh]]$fitted[[ii]][[rr]][ipath[[jj]][[kk]]],
                       pch=psy[[jj]][[rr]][kk],col=pcol_ms[jj])
              }
            }
          }
        }
      }
    }
    abline(0,1,col="blue")
    if( iType ){
      legend("topleft",c("Type = Soft","Type = Hard","Type = Wet"),
             col=pcol,pch=20)
      legend("bottomright",c(ms_paths[[1]][[rr]],ms_paths[[2]][[rr]],
             ms_paths[[3]][[rr]]),col=rep(pcol_ms,
             times=c(length(psy[[1]][[rr]]),length(psy[[2]][[rr]]),
                     length(psy[[3]][[rr]]))),
             pch=c(psy[[1]][[rr]],psy[[2]][[rr]],psy[[3]][[rr]]))
    }
    graphics.off()
  }
}
