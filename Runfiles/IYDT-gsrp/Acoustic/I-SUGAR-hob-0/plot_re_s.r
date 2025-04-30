########################################################################
#                                                                      #
# This file plots 95% credible intervals of source random effects.     #
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

# critical value
cv <- qnorm(.975)

# sources
s_names = names(pred_mle$h[[1]]$mu_b1[[1]])
len_s = length(s_names)

# lower bound
lb_s = pred_mle$h[[1]]$mu_b1[[1]] - cv*pred_mle$h[[1]]$sd_b1[[1]]

# upper bound
ub_s = pred_mle$h[[1]]$mu_b1[[1]] + cv*pred_mle$h[[1]]$sd_b1[[1]]

# labels
lab_s = NULL
for( ii in 1:len_s ){
  if( ub_s[ii] <= 0 || lb_s[ii] >= 0 ){ lab_s = c(lab_s,s_names[ii])
  } else { lab_s = c(lab_s,"") }
}

# Plot of CIs based on MLE
pdf("ci_plot_bs.pdf",height=4,width=8)
plot(0,0,xlim=c(0,len_s+1),ylim=c(min(lb_s),max(ub_s)),
     xlab="Source",ylab="Log Acoustic Impulse",
     type="n",axes=F,frame.plot=TRUE)
axis(1,at=1:length(s_names),
     labels=lab_s,
     las=2,cex.axis=0.75)
axis(2)
for( ii in 1:len_s ){
  if( lab_s[ii] %in% s_names ){ col_s = "red"
  } else { col_s = "pink" }
  segments(ii,lb_s[ii],ii,ub_s[ii],lwd=2,col=col_s)
}
abline(h = 0,lty=2,lwd=2)
graphics.off()
