########################################################################
#                                                                      #
# This file plots the profile log-likelihood function.                 #
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

require(doFuture)
require(iterators)

ngrid <- 51
w <- seq(11,17,length=ngrid)
th0 <- w
nth0 <- length(th0)

# parallel sampling plan
pl = "multicore"
ncor = 30
if( ncor == 1 ){ pl = "sequential"; plan(pl);
} else {
  if( ncor > nth0 ){ ncor = nth0 }
  if( pl != "sequential" ){ plan(pl,workers=ncor)
  } else { plan(pl) }
}

f_pll <- function(x, th0, pc=p_cal)
{
  return(pc$ll_full(c(th0,x), pc))
}

g_pll <- function(x, th0, pc=p_cal)
{
  q <- length(th0)
  g <- pc$gll_full(c(th0,x), pc)
  return(g[-(1:q)])
}

pll_opt = function(th0,x0,pc=p_cal)
{
  t0 <- pc$inv_transform(th0,pc=pc)
  tmp <- optim(x0, fn=f_pll, gr=g_pll, th0=t0, pc=pc,
               method="BFGS",
               control=list(fnscale = -1,maxit=1000))
  c(tmp$value,tmp$convergence)
}
save.image()

x0 <- p_cal$mle_cal
pll_mat <- NULL
pll_mat <- foreach(t0=iter(th0, by='row'),
             .combine='rbind') %dofuture% {
             pll_opt(t0,x0,p_cal)
           } %seed% TRUE
plan(sequential)
save.image()
pll <- pll_mat[,1]
pll_conv <- pll_mat[,2]
for( ii in 1:nth0 ){
  cat(paste("Iteration = ",ii,"; Convergence = ",
            pll_conv[ii],"\n",sep=""))
}

# parallel sampling plan
pl = "multicore"
ncor = 30
if( ncor == 1 ){ pl = "sequential"; plan(pl);
} else {
  if( ncor > nth0 ){ ncor = nth0 }
  if( pl != "sequential" ){ plan(pl,workers=ncor)
  } else { plan(pl) }
}

pll_0 <- NULL
pll_0 <- foreach(t0=iter(th0, by='row'),
           .combine='rbind') %dofuture% {
           t0 = p_cal$inv_transform(t0,pc=p_cal)
           return(p_cal$ll_0(t0,p_cal))
         } %seed% TRUE
plan(sequential)
save.image()

p_mle <- p_cal$tmle_0
p_true <- log(1200000)

lY = seq(11.5,16.5,1)
Y = round(exp(lY)/10^6,1)
pdf("profile_ll.pdf")
par(mgp=c(2,1,0),pty="s")
plot(w,pll,lwd=2,type="l",axes=FALSE,
     xlab="Yield (kt)",ylab="Profile Log-Likelihood",
     main="Profile Log-Likelihood")
box()
axis(1,at=lY,labels=as.character(Y))
axis(2)
abline(v=p_mle,col="green",lwd=2,lty=3)
abline(v=p_true,col="gold",lwd=2,lty=3)
graphics.off()

pdf("profile_ll_w.pdf")
par(mgp=c(2,1,0),pty="s")
plot(w,pll,lwd=2,type="l",
     xlab="log(Yield (kt))",ylab="Profile Log-Likelihood",
     main="Profile Log-Likelihood")
abline(v=p_mle,col="green",lwd=2,lty=3)
abline(v=p_true,col="gold",lwd=2,lty=3)
graphics.off()

pdf("profile_ll_0.pdf")
par(mgp=c(2,1,0),pty="s")
plot(w,pll_0,lwd=2,type="l",axes=FALSE,
     xlab="Yield (kt)",ylab="New Event Log-Likelihood",
     main="New Event Log-Likelihood")
box()
axis(1,at=lY,labels=as.character(Y))
axis(2)
abline(v=p_mle,col="green",lwd=2,lty=3)
abline(v=p_true,col="gold",lwd=2,lty=3)
graphics.off()

t0 <- p_cal$inv_transform(log(1200000),p_cal)
tmp <- optim(x0, fn=f_pll, gr=g_pll, th0=t0, pc=p_cal,
             method="BFGS",
             control=list(fnscale = -1,maxit=1000))
source(paste(gen_dir,"/make_opt.r",sep=""),local=TRUE)
pll_true = make_opt(c(t0,tmp$par), p_cal)
rm(make_opt)
source(paste(gen_dir,"/info_likelihood_0.r",sep=""),local=TRUE)
Sigma_pll_true = info_ll_0(pll_true, p_cal)
rm(info_ll_0)
