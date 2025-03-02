########################################################################
#                                                                      #
# This file plots summaries of the observed data.                      #
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

data_cal <- read.csv("../../Data/seismic_cal.csv")
data_new <- read.csv("../../Data/seismic_new.csv")
f1 <- 10^9; f2 <- 10^6;
y1_all <- exp(c(data_cal$Y1,data_new$Y1))*f1
sy1_all <- y1_all/exp(c(data_cal$W,data_new$W))^(1/3)
y2_all <- exp(c(data_cal$Y2,data_new$Y2))*f2
sy2_all <- y2_all

pdf("data_all_seismic.pdf")
iType = FALSE
pcol = c("black","red","green")
par(mgp=c(2,1,0),pty="s")
for( ii in 1:p_cal$h[[1]]$nsource ){
  if( "Type" %in% names(p_cal$h[[1]]$X[[ii]][[1]]) ){
    tt = as.numeric(as.character(
         p_cal$h[[1]]$X[[ii]][[1]]$Type[1]))
    if( tt > 1 ){ iType = TRUE }
  } else { tt = 1 }
  y1 = exp(p_cal$h[[1]]$Y[[ii]][[1]])*f1
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  if( ii == 1 ){
    plot(y1,y2,
         xlab=expression(paste("Displacement [nm]",sep="")),
         ylab=expression(paste("Velocity [",mu,"m/s]",sep="")),
         xlim=range(y1_all),ylim=range(y2_all),
         pch=20,log="xy",col=pcol[tt],main="Seismic")
  } else {
    if( ii == p_cal$h[[hh]]$nsource ){
      points(y1,y2,
             pch=17,col=pcol[tt])
    } else {
      points(y1,y2,
             pch=20,col=pcol[tt])
    }
  }
}
if( iType ){
  legend("bottomright",c("Type = Soft","Type = Hard","Type = Wet"),
         col=c("black","red","green"),pch=20)
}
graphics.off()

pdf("data_red_seismic.pdf")
pcol = c("black","red","green")
par(mgp=c(2,1,0),pty="s")
init_plot = FALSE
for( ii in 1:p_cal$h[[1]]$nsource ){
  if( "Type" %in% names(p_cal$h[[1]]$X[[ii]][[1]]) ){
    tt = as.numeric(as.character(
         p_cal$h[[1]]$X[[ii]][[1]]$Type[1]))
  } else { tt = 1 }
  y1 = exp(p_cal$h[[1]]$Y[[ii]][[1]])*f1
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  if( ii == p_cal$h[[hh]]$nsource ){
    HOB = data_new$HOB[1]
  } else {
    HOB = p_cal$h[[1]]$X[[ii]][[1]]$HOB[1]
  }
  if( tt == 1 && HOB >= 0 ){
    if( !init_plot ){
      plot(y1,y2,
           xlab=expression(paste("Displacement [nm]",sep="")),
           ylab=expression(paste("Velocity [",mu,"m/s]",sep="")),
           xlim=range(y1_all),ylim=range(y2_all),
           pch=20,log="xy",col=pcol[tt],main="Seismic")
      init_plot = TRUE
    } else {
      if( ii == p_cal$h[[hh]]$nsource ){
        points(y1,y2,
               pch=17,col=pcol[tt])
      } else {
        points(y1,y2,
               pch=20,col=pcol[tt])
      }
    }
  }
}
graphics.off()

pdf("data_all_scaled_seismic.pdf")
iType = FALSE
pcol = c("black","red","green")
par(mgp=c(2,1,0),pty="s")
for( ii in 1:p_cal$h[[1]]$nsource ){
  if( "Type" %in% names(p_cal$h[[1]]$X[[ii]][[1]]) ){
    tt = as.numeric(as.character(
         p_cal$h[[1]]$X[[ii]][[1]]$Type[1]))
    if( tt > 1 ){ iType = TRUE }
  } else { tt = 1 }
  y1 = exp(p_cal$h[[1]]$Y[[ii]][[1]])*f1
  if( ii == p_cal$h[[hh]]$nsource ){
    sy1 = y1/exp(data_new$W[1])^(1/3)
  } else {
    sy1 = y1/exp(p_cal$h[[1]]$X[[ii]][[1]]$W[1])^(1/3)
  }
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  sy2 = y2
  if( ii == 1 ){
    plot(sy1,sy2,
         xlab=expression(paste("Scaled Displacement [nm/kg"^{1/3},"]",
                               sep="")),
         ylab=expression(paste("Velocity [",mu,"m/s]",sep="")),
         xlim=range(sy1_all),ylim=range(sy2_all),
         pch=20,log="xy",col=pcol[tt],main="Seismic")
  } else { 
    if( ii == p_cal$h[[hh]]$nsource ){
      points(sy1,sy2,
             pch=17,col=pcol[tt])
    } else {
      points(sy1,sy2,
             pch=20,col=pcol[tt])
    }
  }
}
if( iType ){
  legend("bottomright",c("Type = Soft","Type = Hard","Type = Wet"),
         col=c("black","red","green"),pch=20)
}
graphics.off()

pdf("data_red_scaled_seismic.pdf")
pcol = c("black","red","green")
par(mgp=c(2,1,0),pty="s")
init_plot = FALSE
for( ii in 1:p_cal$h[[1]]$nsource ){
  if( "Type" %in% names(p_cal$h[[1]]$X[[ii]][[1]]) ){
    tt = as.numeric(as.character(
         p_cal$h[[1]]$X[[ii]][[1]]$Type[1]))
  } else { tt = 1 }
  y1 = exp(p_cal$h[[1]]$Y[[ii]][[1]])*f1
  if( ii == p_cal$h[[hh]]$nsource ){
    sy1 = y1/exp(data_new$W[1])^(1/3)
    HOB = data_new$HOB[1]
  } else {
    sy1 = y1/exp(p_cal$h[[1]]$X[[ii]][[1]]$W[1])^(1/3)
    HOB = p_cal$h[[1]]$X[[ii]][[1]]$HOB[1]
  }
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  sy2 = y2
  if( tt == 1 && HOB >= 0 ){
    if( !init_plot ){
      plot(sy1,sy2,
           xlab=expression(paste("Scaled Displacement [nm/kg"^{1/3},"]",
                                 sep="")),
           ylab=expression(paste("Velocity [",mu,"m/s]",sep="")),
           xlim=range(sy1_all),ylim=range(sy2_all),
           pch=20,log="xy",col=pcol[tt],main="Seismic")
      init_plot = TRUE
    } else {
      if( ii == p_cal$h[[hh]]$nsource ){
        points(sy1,sy2,
               pch=17,col=pcol[tt])
      } else {
        points(sy1,sy2,
               pch=20,col=pcol[tt])
      }
    }
  }
}
graphics.off()
