########################################################################
#                                                                      #
# This file plots the fitted vs. the observed data, where the fits are #
# obtained by evaluating the forward model using the coefficient MLEs  #
# and the observed data are adjusted for path effects to highlight the #
# presence of source bias.                                             #
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

ms_sources_1 = vector("list",2)
pch_1 = vector("list",2)
ms_sources_2 = vector("list",2)
pch_2 = vector("list",2)
ms_sources_3 = vector("list",2)
pch_3 = vector("list",2)

nsg = 3 # number of source groups
ms_sources_1[[1]] = c("HTA-5","HRI-7","HRIII-2")
pch_1[[1]] = c(7,9,10)
ms_sources_2[[1]] = c("HRIII-3")
pch_2[[1]] = 2
ms_sources_3[[1]] = c("HTP-9","HTP-8","HTP-7")
pch_3[[1]] = c(15,17,19)

ms_sources_1[[2]] = c("HRIII-2","HTA-2","HTA-4","HTA-1")
pch_1[[2]] = c(7,9,10,14)
ms_sources_2[[2]] = c("SAY-1")
pch_2[[2]] = 2
ms_sources_3[[2]] = c("HTP-5","HTP-6")
pch_3[[2]] = c(15,17)

pcol = c("gray","pink","lightgreen")
pcol_ms = c("black","red","darkgreen")
for( hh in 1:p_cal$H ){
  obs = vector("list",p_cal$h[[hh]]$Rh)
  fit = vector("list",p_cal$h[[hh]]$Rh)
  for( rr in 1:p_cal$h[[hh]]$Rh ){
    for( ii in 1:p_cal$h[[hh]]$nsource ){
      # adjusted for path effects
      obs[[rr]] = c(obs[[rr]],pred_mle$h[[hh]]$observed[[ii]][[rr]]-
                              pred_mle$h[[hh]]$path_effects[[ii]][[rr]])
      fit[[rr]] = c(fit[[rr]],pred_mle$h[[hh]]$fitted[[ii]][[rr]])
    }
  }
  for( rr in 1:p_cal$h[[hh]]$Rh ){
    pdf(paste("fit_vs_obs_h",hh,"_r",rr,"_1.pdf",sep=""))
    iType = FALSE
    for( ii in 1:p_cal$h[[hh]]$nsource ){
      if( p_cal$h[[hh]]$n[[ii]][rr] > 0 ){
        if( "Type" %in% names(p_cal$h[[hh]]$X[[ii]][[rr]]) ){
          tt = as.numeric(as.character(
               p_cal$h[[hh]]$X[[ii]][[rr]]$Type[1]))
          if( tt > 1 ){ iType = TRUE }
        } else { tt = 1 }
        ppch = 20
        ind_ms = vector("list",nsg)
        for( jj in 1:nsg ){
          eval(parse(text=paste("ind_ms[[",jj,"]] = ",
            "which(ms_sources_",jj,"[[",rr,"]] %in% p_cal$h[[",hh,
                   "]]$X[[",ii,"]][[",rr,"]]$Source[1])",
            sep="")))
          if( length(ind_ms[[jj]]) > 0 ){
            eval(parse(text=paste("ppch = pch_",jj,"[[",rr,"]][ind_ms[[",
                                  jj,"]]]",sep="")))
            pcolor = pcol_ms[jj]
          }
        }
        if( ppch == 20 ){ pcolor = pcol[tt] }
        if( ii == 1 ){
          plot(pred_mle$h[[hh]]$observed[[ii]][[rr]] -
               pred_mle$h[[hh]]$path_effects[[ii]][[rr]],
               pred_mle$h[[hh]]$fitted[[ii]][[rr]],
               xlim=range(obs[[rr]]),ylim=range(fit[[rr]]),
               xlab="observed",ylab="fitted",
               pch=ppch,col=pcolor)
        } else {
          points(pred_mle$h[[hh]]$observed[[ii]][[rr]] -
                 pred_mle$h[[hh]]$path_effects[[ii]][[rr]],
                 pred_mle$h[[hh]]$fitted[[ii]][[rr]],
                 pch=ppch,col=pcolor)
        }
      }
    }
    abline(0,1,col="blue")
    if( iType ){
      legend("topleft",c("Type = Soft","Type = Hard","Type = Wet"),
             col=pcol,pch=20)
      legend("bottomright",c(ms_sources_1[[rr]],ms_sources_2[[rr]],
             ms_sources_3[[rr]]),col=rep(pcol_ms,
             times=c(length(pch_1[[rr]]),length(pch_2[[rr]]),
                     length(pch_3[[rr]]))),
             pch=c(pch_1[[rr]],pch_2[[rr]],pch_3[[rr]]))
    }
    graphics.off()
  }
}
