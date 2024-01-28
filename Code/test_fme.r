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

set.seed(100)

ll <- function(x)
{
  return(dnorm(y,x,sig,log=TRUE))
}

lp <- function(x)
{
  return(dnorm(x,mu_0,sig_0,log=TRUE))
}

mu_0 = -3 
sig_0 = 4

sig = 2
y = rnorm(1,mu_0,sig)
cat(paste("y = ",y,"\n",sep=""))

x0 = 0

require(FME)
ll_mod = function(x){ -2*ll(x) }
lp_mod = function(x){ -2*lp(x) }
burnin = 100000
niter = burnin+1000000
samp_po <- modMCMC(ll_mod,x0,prior=lp_mod,
                   jump=0.1,
                   niter=niter,burninlength=burnin,
                   updatecov=100,ntrydr=2,
                   verbose=TRUE)

pdf("posterior.pdf")
par(mgp=c(2,1,0),pty="s")
hist(samp_po$pars,prob=TRUE,xlab="x",
     main="Histogram of Sampled Posterior")
sig2_p = sig^2*sig_0^2/(sig^2+sig_0^2)
mu_p = sig2_p*(y/sig^2+mu_0/sig_0^2)
sig_p = sqrt(sig2_p)
xx = seq(mu_p-4*sig_p,mu_p+4*sig_p,length=1001)
lines(xx,dnorm(xx,mu_p,sig_p),lwd=2,col="green")
legend("topleft","true posterior",lwd=2,col="green")
graphics.off()

cat(paste("Sample posterior mean = ",
          mean(samp_po$pars),"\n",sep=""))
cat(paste("Posterior mean = ",mu_p,"\n",sep=""))
cat(paste("Sample posterior SD = ",
          sd(samp_po$pars),"\n",sep=""))
cat(paste("Posterior SD = ",sig_p,"\n",sep=""))
