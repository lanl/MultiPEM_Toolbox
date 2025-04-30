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

y_eff <- function(h_s,tmin,beta)
{
  cf <- 1+beta[3]*exp(-(h_s/beta[4])^2)
  tmin_eff <- tmin*cf
  y_eff <- exp((log(tmin_eff) - beta[1])/beta[2])/10^6
  return(y_eff)
}

notExp <- function(x)
{
  f <- x
  ind <- x > 1
  f[ind] <- exp(1) * (x[ind]^2 + 1)/2
  ind <- (x <= 1) & (x > -1)
  f[ind] <- exp(x[ind])
  ind <- (x <= -1)
  x[ind] <- -x[ind]
  f[ind] <- exp(1) * (x[ind]^2 + 1)/2
  f[ind] <- 1/f[ind]
  f
}

data_cal <- read.csv("../../Data/optical_cal.csv")
data_new <- read.csv("../../Data/optical_new.csv")
f1 <- 10^3; f2 <- 10^3;
w_all <- c(data_cal$W,data_new$W)
y_all <- exp(w_all)/10^6
y1_all <- exp(c(data_cal$Y1,data_new$Y1))*f1
sy1_all <- y1_all/exp(w_all/3)/
           exp(c(data_cal$logPressureSc,data_new$logPressureSc))^(1/3)*
           exp(c(data_cal$logTempSc,data_new$logTempSc))^(1/2)
y2_all <- exp(c(data_cal$Y2,data_new$Y2))*f2
sy2_all <- y2_all/exp(w_all/3)/
           exp(c(data_cal$logPressureSc,data_new$logPressureSc))^(1/3)*
           exp(c(data_cal$logTempSc,data_new$logTempSc))^(1/2)

# Whitaker-Symbalisty coefficient estimates
beta <- c(-10.5121285,0.3674431,-0.4621065,-0.8703267)
beta[3] = notExp(beta[3])
beta[4] = notExp(beta[4])

tmin_all <- exp(beta[1]+beta[2]*w_all)
h_all <- c(data_cal$HOB,data_new$HOB)
h_s <- h_all/exp(w_all/3)
yeff_all <- y_eff(h_s,tmin_all,beta)
si_all <- (yeff_all >= 2*y_all)

pdf("data_all_optical.pdf")
iType = FALSE
pcol = "black"
pcol_eff = "red"
par(mgp=c(2,1,0),pty="s")
for( ii in 1:p_cal$h[[1]]$nsource ){
  if( "Type" %in% names(p_cal$h[[1]]$X[[ii]][[1]]) ){
    tt = as.numeric(as.character(
         p_cal$h[[1]]$X[[ii]][[1]]$Type[1]))
    if( tt > 1 ){ iType = TRUE }
  } else { tt = 1 }
  y1 = exp(p_cal$h[[1]]$Y[[ii]][[1]])*f1
  si = si_all[1:length(y1)]
  si_all = si_all[-(1:length(si))]
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  if( ii == 1 ){
    if( any(si) ){
      plot(y1[si],y2[si],
           xlab=expression('t'[min]*' [ms]'),
           ylab=expression('t'[max]^(2)*' [ms]'),
           xlim=range(y1_all),ylim=range(y2_all),
           pch=20,log="xy",col=pcol_eff[tt],main="Optical")
      if( !all(si) ){ points(y1[!si],y2[!si],pch=20,col=pcol[tt]) }
    } else {
      plot(y1[!si],y2[!si],
           xlab=expression('t'[min]*' [ms]'),
           ylab=expression('t'[max]^(2)*' [ms]'),
           xlim=range(y1_all),ylim=range(y2_all),
           pch=20,log="xy",col=pcol[tt],main="Optical")
    }
  } else {
    if( ii == p_cal$h[[hh]]$nsource ){
      if( any(si) ){
        points(y1[si],y2[si],pch=17,col=pcol_eff[tt])
        if( !all(si) ){ points(y1[!si],y2[!si],pch=17,col=pcol[tt]) }
      } else { points(y1[!si],y2[!si],pch=17,col=pcol[tt]) }
    } else {
      if( any(si) ){
        points(y1[si],y2[si],pch=20,col=pcol_eff[tt])
        if( !all(si) ){ points(y1[!si],y2[!si],pch=20,col=pcol[tt]) }
      } else { points(y1[!si],y2[!si],pch=20,col=pcol[tt]) }
    }
  }
}
if( iType ){
  legend("bottomright",c("Type = Soft","Type = Hard","Type = Wet"),
         col=c("black","red","green"),pch=20)
} else {
  legend("bottomright",c(expression(y[eff] >= '2y'),
                         expression(y[eff] < '2y')),
         pch=20,col=c(pcol_eff,pcol),bty="n")
}
graphics.off()

si_all <- (yeff_all >= 2*y_all)

pdf("data_red_optical.pdf")
init_plot = FALSE
pcol = "black"
pcol_eff = "red"
par(mgp=c(2,1,0),pty="s")
for( ii in 1:p_cal$h[[1]]$nsource ){
  if( "Type" %in% names(p_cal$h[[1]]$X[[ii]][[1]]) ){
    tt = as.numeric(as.character(
         p_cal$h[[1]]$X[[ii]][[1]]$Type[1]))
  } else { tt = 1 }
  y1 = exp(p_cal$h[[1]]$Y[[ii]][[1]])*f1
  si = si_all[1:length(y1)]
  si_all = si_all[-(1:length(si))]
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  if( ii == p_cal$h[[hh]]$nsource ){
    HOB = data_new$HOB[1]
  } else {
    HOB = p_cal$h[[1]]$X[[ii]][[1]]$HOB[1]
  }
  if( tt == 1 && HOB >= 0 ){
    if( !init_plot ){
      if( any(si) ){
        plot(y1[si],y2[si],
             xlab=expression('t'[min]*' [ms]'),
             ylab=expression('t'[max]^(2)*' [ms]'),
             xlim=range(y1_all),ylim=range(y2_all),
             pch=20,log="xy",col=pcol_eff[tt],main="Optical")
        if( !all(si) ){ points(y1[!si],y2[!si],pch=20,col=pcol[tt]) }
      } else {
        plot(y1[!si],y2[!si],
             xlab=expression('t'[min]*' [ms]'),
             ylab=expression('t'[max]^(2)*' [ms]'),
             xlim=range(y1_all),ylim=range(y2_all),
             pch=20,log="xy",col=pcol[tt],main="Optical")
      }
      init_plot = TRUE
    } else {
      if( ii == p_cal$h[[hh]]$nsource ){
        if( any(si) ){
          points(y1[si],y2[si],pch=17,col=pcol_eff[tt])
          if( !all(si) ){ points(y1[!si],y2[!si],pch=17,col=pcol[tt]) }
        } else { points(y1[!si],y2[!si],pch=17,col=pcol[tt]) }
      } else {
        if( any(si) ){
          points(y1[si],y2[si],pch=20,col=pcol_eff[tt])
          if( !all(si) ){ points(y1[!si],y2[!si],pch=20,col=pcol[tt]) }
        } else { points(y1[!si],y2[!si],pch=20,col=pcol[tt]) }
      }
    }
  }
}
legend("bottomright",c(expression(y[eff] >= '2y'),
                       expression(y[eff] < '2y')),
       pch=20,col=c(pcol_eff,pcol),bty="n")
graphics.off()

si_all <- (yeff_all >= 2*y_all)

pdf("data_all_scaled_optical.pdf")
iType = FALSE
pcol = "black"
pcol_eff = "red"
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
  sy1 = sy1/exp(p_cal$h[[1]]$X[[ii]][[1]]$logPressureSc[1])^(1/3)*
        exp(p_cal$h[[1]]$X[[ii]][[1]]$logTempSc[1])^(1/2)
  si = si_all[1:length(y1)]
  si_all = si_all[-(1:length(si))]
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  if( ii == p_cal$h[[hh]]$nsource ){
    sy2 = y2/exp(data_new$W[1])^(1/3)
  } else {
    sy2 = y2/exp(p_cal$h[[1]]$X[[ii]][[2]]$W[1])^(1/3)
  }
  sy2 = sy2/exp(p_cal$h[[1]]$X[[ii]][[2]]$logPressureSc[1])^(1/3)*
        exp(p_cal$h[[1]]$X[[ii]][[2]]$logTempSc[1])^(1/2)
  if( ii == 1 ){
    if( any(si) ){
      plot(sy1[si],sy2[si],
           xlab=expression('Scaled t'[min]*' [ms/kg'^{1/3}*']'),
           ylab=expression('Scaled t'[max]^(2)*' [ms/kg'^{1/3}*']'),
           xlim=range(sy1_all),ylim=range(sy2_all),
           pch=20,log="xy",col=pcol_eff[tt],main="Optical")
      if( !all(si) ){ points(sy1[!si],sy2[!si],pch=20,col=pcol[tt]) }
    } else {
      plot(sy1[!si],sy2[!si],
           xlab=expression('Scaled t'[min]*' [ms/kg'^{1/3}*']'),
           ylab=expression('Scaled t'[max]^(2)*' [ms/kg'^{1/3}*']'),
           xlim=range(sy1_all),ylim=range(sy2_all),
           pch=20,log="xy",col=pcol[tt],main="Optical")
    }
  } else { 
    if( ii == p_cal$h[[hh]]$nsource ){
      if( any(si) ){
        points(sy1[si],sy2[si],pch=17,col=pcol_eff[tt])
        if( !all(si) ){ points(sy1[!si],sy2[!si],pch=17,col=pcol[tt]) }
      } else { points(sy1[!si],sy2[!si],pch=17,col=pcol[tt]) }
    } else {
      if( any(si) ){
        points(sy1[si],sy2[si],pch=20,col=pcol_eff[tt])
        if( !all(si) ){ points(sy1[!si],sy2[!si],pch=20,col=pcol[tt]) }
      } else { points(sy1[!si],sy2[!si],pch=20,col=pcol[tt]) }
    }
  }
}
if( iType ){
  legend("bottomright",c("Type = Soft","Type = Hard","Type = Wet"),
         col=c("black","red","green"),pch=20)
} else {
#  legend("bottomright",c(expression(y[eff] >= '2y'),
#                         expression(y[eff] < '2y')),
#         pch=20,col=c(pcol_eff,pcol),bty="n")
}
graphics.off()

si_all <- (yeff_all >= 2*y_all)

pdf("data_red_scaled_optical.pdf")
init_plot = FALSE
pcol = "black"
pcol_eff = "red"
par(mgp=c(2,1,0),pty="s")
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
  sy1 = sy1/exp(p_cal$h[[1]]$X[[ii]][[1]]$logPressureSc[1])^(1/3)*
        exp(p_cal$h[[1]]$X[[ii]][[1]]$logTempSc[1])^(1/2)
  si = si_all[1:length(y1)]
  si_all = si_all[-(1:length(si))]
  y2 = exp(p_cal$h[[1]]$Y[[ii]][[2]])*f2
  if( ii == p_cal$h[[hh]]$nsource ){
    sy2 = y2/exp(data_new$W[1])^(1/3)
  } else {
    sy2 = y2/exp(p_cal$h[[1]]$X[[ii]][[2]]$W[1])^(1/3)
  }
  sy2 = sy2/exp(p_cal$h[[1]]$X[[ii]][[2]]$logPressureSc[1])^(1/3)*
        exp(p_cal$h[[1]]$X[[ii]][[2]]$logTempSc[1])^(1/2)
  if( tt == 1 && HOB >= 0 ){
    if( !init_plot ){
      if( any(si) ){
        plot(sy1[si],sy2[si],
             xlab=expression('Scaled t'[min]*' [ms/kg'^{1/3}*']'),
             ylab=expression('Scaled t'[max]^(2)*' [ms/kg'^{1/3}*']'),
             xlim=range(sy1_all),ylim=range(sy2_all),
             pch=20,log="xy",col=pcol_eff[tt],main="Optical")
        if( !all(si) ){ points(sy1[!si],sy2[!si],pch=20,col=pcol[tt]) }
      } else {
        plot(sy1[!si],sy2[!si],
             xlab=expression('Scaled t'[min]*' [ms/kg'^{1/3}*']'),
             ylab=expression('Scaled t'[max]^(2)*' [ms/kg'^{1/3}*']'),
             xlim=range(sy1_all),ylim=range(sy2_all),
             pch=20,log="xy",col=pcol[tt],main="Optical")
      }
      init_plot = TRUE
    } else {
      if( ii == p_cal$h[[hh]]$nsource ){
        if( any(si) ){
          points(sy1[si],sy2[si],pch=17,col=pcol_eff[tt])
          if( !all(si) ){
            points(sy1[!si],sy2[!si],pch=17,col=pcol[tt])
          }
        } else { points(sy1[!si],sy2[!si],pch=17,col=pcol[tt]) }
      } else {
        if( any(si) ){
          points(sy1[si],sy2[si],pch=20,col=pcol_eff[tt])
          if( !all(si) ){
            points(sy1[!si],sy2[!si],pch=20,col=pcol[tt])
          }
        } else { points(sy1[!si],sy2[!si],pch=20,col=pcol[tt]) }
      }
    }
  }
}
#legend("bottomright",c(expression(y[eff] >= '2y'),
#                       expression(y[eff] < '2y')),
#       pch=20,col=c(pcol_eff,pcol),bty="n")
graphics.off()
