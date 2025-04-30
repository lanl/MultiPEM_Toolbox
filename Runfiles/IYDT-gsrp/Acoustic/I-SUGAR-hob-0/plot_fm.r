########################################################################
#                                                                      #
# This file plots contours of the forward model for one sensor         #
# observing the new event.                                             #
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

# set up list "pm" with known quantities passed to
# forward model
# named objects
pnames = names(p_cal$h[[1]])
pm = new.env(hash=TRUE)
if( exists("notExp",where=p_cal,inherits=FALSE) ){
  pm$notExp = p_cal$notExp
}
if( !is.null(p_cal$h[[1]]$llpars[["yield_scaling"]]) ){
  pm$yield_scaling = p_cal$h[[1]]$llpars$yield_scaling
}
if( !is.null(p_cal$h[[1]]$llpars[["pressure_scaling"]]) ){
  pm$pressure_scaling = p_cal$h[[1]]$llpars$pressure_scaling
}
if( !is.null(p_cal$h[[1]]$llpars[["temp_scaling"]]) ){
 pm$temp_scaling = p_cal$h[[1]]$llpars$temp_scaling
}
# extract forward model parameters
if( p_cal$pbeta > 0 && any(p_cal$h[[1]]$pbeta > 0) ){
  pbeta = sum(p_cal$h[[1]]$pbeta)
  beta = p_cal$h[[1]]$beta
} else { pbeta = 0 }
if( p_cal$ptbeta > 0 && "ptbeta" %in% pnames ){
  ptbeta = sum(p_cal$h[[1]]$ptbeta)
  tbeta = p_cal$h[[1]]$tbeta
} else { ptbeta = 0 }
# number of responses
Rh = p_cal$h[[1]]$Rh
# number of responses for source "0"
n_h0 = p_cal$h[[1]]$n0
# named parameters in forward model call
if( "theta0_names" %in% pnames ){
  pm$theta_names = p_cal$h[[1]]$theta0_names
}
# extract forward model parameters for emplacement condition "tt"
# associated with source "0"
if( ptbeta > 0 ){
  for( rr in 1:Rh ){
    if( n_h0[rr] > 0 ){
      tt = as.numeric(as.character(p_cal$h[[1]]$X0[[rr]]$Type[1]))
      break
    } else { next }
  }
  st_betat = 0
  if( tt > 1 ){ st_betat = sum(p_cal$h[[1]]$ptbeta[1:(tt-1)]) }
  if( p_cal$h[[1]]$ptbeta[tt] > 0 ){
    betat = tbeta[st_betat+(1:p_cal$h[[1]]$ptbeta[tt])]
  }
} else { tt = 1 }
ng = 101
w = seq(11,17,length=ng)
h = seq(-10,160,length=ng)
y = vector("list",Rh)
f0 = p_cal$ffm$f0_a
rnames = c("Acoustic Impulse","Positive Phase Duration")
p_mle <- p_cal$tmle_0
p_true <- c(log(1200000),1.0668)
lY = seq(11.5,16.5,1)
Y = round(exp(lY)/10^6,1)
# iterate over responses
for(rr in 1:2){
  # covariate matrix for source "0"
  pm$X = p_cal$h[[1]]$X0[[rr]][1,]
  # extract forward model parameters for response "rr"
  if( pbeta > 0 && p_cal$h[[1]]$pbeta[rr] > 0 ){
    st_beta = 0
    if( rr > 1 ){ 
      st_beta = sum(p_cal$h[[1]]$pbeta[1:(rr-1)])
    }
    betar = beta[st_beta+(1:p_cal$h[[1]]$pbeta[rr])]
  } else { betar = NULL }
  if( ptbeta > 0 && p_cal$h[[1]]$pbetat[[tt]][rr] > 0 ){
    st_betatr = 0
    if( rr > 1 ){
      st_betatr = sum(p_cal$h[[1]]$pbetat[[tt]][1:(rr-1)])
    }
    betatr = betat[st_betatr+(1:p_cal$h[[1]]$pbetat[[tt]][rr])]
  } else { betatr = NULL }
  if( is.null(betatr) ){ beta_t = betar }
  if( is.null(betar) ){ beta_t = betatr }
  if( !is.null(betar) && !is.null(betatr) ){
    beta_t = numeric(p_cal$h[[1]]$pbeta[rr]+
                     p_cal$h[[1]]$pbetat[[tt]][rr])
    beta_t[p_cal$h[[1]]$ibetar[[(tt-1)*Rh+rr]]] = betar
    beta_t[p_cal$h[[1]]$ibetatr[[(tt-1)*Rh+rr]]] = betatr
  }
  pm$beta = beta_t
  if( "iResponse" %in% pnames ){
    pm$iresp = p_cal$h[[1]]$iResponse[rr]
  }
  y[[rr]] = matrix(0,ng,ng)
  for(ii in 1:ng){
    for(jj in 1:ng){
      y[[rr]][ii,jj] = f0(c(w[ii],h[jj]),pm)
    }
  }
  pdf(paste(rnames[rr],".pdf",sep=""))
  par(mgp=c(2,1,0),pty="s")
  contour(w,h,y[[rr]],xlab="Yield (kt)",ylab="HOB/DOB (m)",lwd=2,
          xlim=c(11,17),
          axes=FALSE)
  box()
  axis(1,at=lY,labels=as.character(Y))
  axis(2)
  points(p_mle[1],p_mle[2],pch=20,col="green",cex=2)
  points(p_true[1],p_true[2],pch=20,col="gold",cex=2)
  title(rnames[rr],line=1)
  graphics.off()
}
