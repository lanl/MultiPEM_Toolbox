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

source("../../../../Code/predict.r")

pred_mle = predict(c(p_cal$mle,p_cal$mle_cal), p_cal)

pcol = c("black","red","green")
for( hh in 1:p_cal$H ){
  obs = vector("list",p_cal$h[[hh]]$Rh)
  fit = vector("list",p_cal$h[[hh]]$Rh)
  for( rr in 1:p_cal$h[[hh]]$Rh ){
    for( ii in 1:p_cal$h[[hh]]$nsource ){
      obs[[rr]] = c(obs[[rr]],pred_mle$h[[hh]]$observed[[ii]][[rr]])
      fit[[rr]] = c(fit[[rr]],pred_mle$h[[hh]]$fitted[[ii]][[rr]])
    }
  }
  for( rr in 1:p_cal$h[[hh]]$Rh ){
    pdf(paste("fit_vs_obs_h",hh,"_r",rr,".pdf",sep=""))
    iType = FALSE
    for( ii in 1:p_cal$h[[hh]]$nsource ){
      if( p_cal$h[[hh]]$n[[ii]][rr] > 0 ){
        if( "Type" %in% names(p_cal$h[[hh]]$X[[ii]][[rr]]) ){
          tt = as.numeric(as.character(
               p_cal$h[[hh]]$X[[ii]][[rr]]$Type[1]))
          if( tt > 1 ){ iType = TRUE }
        } else { tt = 1 }
        if( ii == 1 ){
          plot(pred_mle$h[[hh]]$observed[[ii]][[rr]],
               pred_mle$h[[hh]]$fitted[[ii]][[rr]],
               xlim=range(obs[[rr]]),ylim=range(fit[[rr]]),
               xlab="observed",ylab="fitted",
               pch=20,col=pcol[tt])
        } else {
          if( ii == p_cal$h[[hh]]$nsource ){
            points(pred_mle$h[[hh]]$observed[[ii]][[rr]],
                   pred_mle$h[[hh]]$fitted[[ii]][[rr]],
                   pch=17,col=pcol[tt])
          } else {
            points(pred_mle$h[[hh]]$observed[[ii]][[rr]],
                   pred_mle$h[[hh]]$fitted[[ii]][[rr]],
                   pch=20,col=pcol[tt])
          }
        }
      }
    }
    abline(0,1,col="blue")
    if( iType ){
      legend("topleft",c("Type = Soft","Type = Hard","Type = Wet"),
             col=c("black","red","green"),pch=20)
    }
    graphics.off()
  }
}
