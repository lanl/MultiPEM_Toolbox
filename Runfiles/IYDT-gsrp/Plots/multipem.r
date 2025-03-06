########################################################################
#                                                                      #
# This file plots 95% confidence intervals for yield based on its      #
# MLE for various analyses, as well as 95% credible intervals for      #
# yield based on posterior samples from various analyses.              #
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

# Plot of CIs based on MLE
pdf("ci_plot_mle.pdf",height=4,width=8)
plot(0,0,xlim=c(0.5,5.5),ylim=c(0,17.75),
     xlab="Phenomenology",ylab="Yield (kt)",
     type="n",axes=F,frame.plot=TRUE)
axis(1,at=c(1,1.5,2,3,3.5,4,5),
     labels=c("S","A","S + A","SE","O","SE + O","MultiPEM"),
     las=2,cex.axis=0.75)
axis(2)
segments(1,0.1,1,4.67,lwd=4,col="red")
segments(0.9,0.1,1.1,0.1,lwd=4,col="red")
segments(0.9,4.67,1.1,4.67,lwd=4,col="red")
segments(1.5,0.23,1.5,3.60,lwd=4,col="blue")
segments(1.4,0.23,1.6,0.23,lwd=4,col="blue")
segments(1.4,3.60,1.6,3.60,lwd=4,col="blue")
segments(2,0.32,2,1.67,lwd=4,col="gray50")
segments(1.9,0.32,2.1,0.32,lwd=4,col="gray50")
segments(1.9,1.67,2.1,1.67,lwd=4,col="gray50")
segments(3,0.5,3,2.5,lwd=4,col="violet")
segments(2.9,0.5,3.1,0.5,lwd=4,col="violet")
segments(2.9,2.5,3.1,2.5,lwd=4,col="violet")
segments(3.5,0.15,3.5,13.34,lwd=4,col="green")
segments(3.4,0.15,3.6,0.15,lwd=4,col="green")
segments(3.4,13.34,3.6,13.34,lwd=4,col="green")
segments(4,0.53,4,2.42,lwd=4,col="gray50")
segments(3.9,0.53,4.1,0.53,lwd=4,col="gray50")
segments(3.9,2.42,4.1,2.42,lwd=4,col="gray50")
segments(5,0.56,5,1.58,lwd=4,col="black")
segments(4.9,0.56,5.1,0.56,lwd=4,col="black")
segments(4.9,1.58,5.1,1.58,lwd=4,col="black")
abline(h = 1.2,lty=2,lwd=2,col="orange")
legend("topright",c("Seismic (S)","Acoustic (A)",
                    "Surface Effects (SE)","Optical (O)",
                    "Sugar (1.2 kt)"),
       col=c("red","blue","violet","green","orange"),
       lwd=2,lty=c(1,1,1,1,2),cex=0.5)
graphics.off()

# Plot of CIs based on MLE
pdf("ci_plot_mle_v.pdf",height=8,width=4)
plot(0,0,xlim=c(0,17.75),ylim=c(0.5,5.5),
     xlab="Yield (kt)",ylab="Phenomenology",
     type="n",axes=F,frame.plot=TRUE)
axis(1)
axis(2,at=c(1,1.5,2,3,3.5,4,5),
     labels=c("S","A","S + A","SE","O","SE + O","MultiPEM"),
     las=2,cex.axis=0.75)
segments(0.1,1,4.67,1,lwd=4,col="red")
segments(0.1,0.9,0.1,1.1,lwd=4,col="red")
segments(4.67,0.9,4.67,1.1,lwd=4,col="red")
segments(0.23,1.5,3.60,1.5,lwd=4,col="blue")
segments(0.23,1.4,0.23,1.6,lwd=4,col="blue")
segments(3.60,1.4,3.60,1.6,lwd=4,col="blue")
segments(0.32,2,1.67,2,lwd=4,col="gray50")
segments(0.32,1.9,0.32,2.1,lwd=4,col="gray50")
segments(1.67,1.9,1.67,2.1,lwd=4,col="gray50")
segments(0.5,3,2.5,3,lwd=4,col="violet")
segments(0.5,2.9,0.5,3.1,lwd=4,col="violet")
segments(2.5,2.9,2.5,3.1,lwd=4,col="violet")
segments(0.15,3.5,13.34,3.5,lwd=4,col="green")
segments(0.15,3.4,0.15,3.6,lwd=4,col="green")
segments(13.34,3.4,13.34,3.6,lwd=4,col="green")
segments(0.53,4,2.42,4,lwd=4,col="gray50")
segments(0.53,3.9,0.53,4.1,lwd=4,col="gray50")
segments(2.42,3.9,2.42,4.1,lwd=4,col="gray50")
segments(0.56,5,1.58,5,lwd=4,col="black")
segments(0.56,4.9,0.56,5.1,lwd=4,col="black")
segments(1.58,4.9,1.58,5.1,lwd=4,col="black")
abline(v = 1.2,lty=2,lwd=2,col="orange")
legend("topright",c("Seismic (S)","Acoustic (A)",
                    "Surface Effects (SE)","Optical (O)",
                    "Sugar (1.2 kt)"),
       col=c("red","blue","violet","green","orange"),
       lwd=2,lty=c(1,1,1,1,2),cex=0.5)
graphics.off()

# Plot of CIs based on posterior samples
pdf("ci_plot_bayes.pdf",height=4,width=8)
plot(0,0,xlim=c(0.5,5.5),ylim=c(0,17.75),
     xlab="Phenomenology",ylab="Yield (kt)",
     type="n",axes=F,frame.plot=TRUE)
axis(1,at=c(1,1.5,2,3,3.5,4,5),
     labels=c("S","A","S + A","SE","O","SE + O","MultiPEM"),
     las=2,cex.axis=0.75)
axis(2)
segments(1,0.26,1,3.00,lwd=4,col="red")
segments(0.9,0.26,1.1,0.26,lwd=4,col="red")
segments(0.9,3.00,1.1,3.00,lwd=4,col="red")
segments(1.5,0.24,1.5,2.33,lwd=4,col="blue")
segments(1.4,0.24,1.6,0.24,lwd=4,col="blue")
segments(1.4,2.33,1.6,2.33,lwd=4,col="blue")
segments(2,0.32,2,1.70,lwd=4,col="gray50")
segments(1.9,0.32,2.1,0.32,lwd=4,col="gray50")
segments(1.9,1.70,2.1,1.70,lwd=4,col="gray50")
segments(3,0.27,3,4.63,lwd=4,col="violet")
segments(2.9,0.27,3.1,0.27,lwd=4,col="violet")
segments(2.9,4.63,3.1,4.63,lwd=4,col="violet")
segments(3.5,0.15,3.5,17.46,lwd=4,col="green")
segments(3.4,0.15,3.6,0.15,lwd=4,col="green")
segments(3.4,17.46,3.6,17.46,lwd=4,col="green")
segments(4,0.41,4,3.69,lwd=4,col="gray50")
segments(3.9,0.41,4.1,0.41,lwd=4,col="gray50")
segments(3.9,3.69,4.1,3.69,lwd=4,col="gray50")
segments(5,0.45,5,1.79,lwd=4,col="black")
segments(4.9,0.45,5.1,0.45,lwd=4,col="black")
segments(4.9,1.79,5.1,1.79,lwd=4,col="black")
abline(h = 1.2,lty=2,lwd=2,col="orange")
legend("topright",c("Seismic (S)","Acoustic (A)",
                    "Surface Effects (SE)","Optical (O)",
                    "Sugar (1.2 kt)"),
       col=c("red","blue","violet","green","orange"),
       lwd=2,lty=c(1,1,1,1,2),cex=0.5)
graphics.off()
