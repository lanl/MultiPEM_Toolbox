########################################################################
#                                                                      #
# This file plots 95% credible intervals of path random effects.       #
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

# paths
p_names = names(pred_mle$h[[1]]$mu_b2[[1]])
len_p = length(p_names)

# lower bound
lb_p = pred_mle$h[[1]]$mu_b2[[1]] - cv*pred_mle$h[[1]]$sd_b2[[1]]

# upper bound
ub_p = pred_mle$h[[1]]$mu_b2[[1]] + cv*pred_mle$h[[1]]$sd_b2[[1]]

# labels
lab_p = NULL
for( ii in 1:len_p ){
  if( ub_p[ii] <= 0 || lb_p[ii] >= 0 ){ lab_p = c(lab_p,p_names[ii])
  } else { lab_p = c(lab_p,"") }
}
names(lab_p) = p_names

# thinning
ind_p = sort(sample(p_names,round(0.7*len_p)))
while( any( lab_p[ind_p] %in% p_names ) ){
  ind_p = sort(sample(p_names,round(0.7*len_p)))
}
nind_p = !(p_names %in% ind_p)
p_names = p_names[nind_p]
len_p = length(p_names)
lb_p = lb_p[nind_p]
ub_p = ub_p[nind_p]
lab_p = lab_p[nind_p]

# Plot of CIs based on MLE
pdf("ci_plot_bp.pdf",height=4,width=8)
plot(0,0,xlim=c(0,len_p+1),ylim=c(min(lb_p),max(ub_p)),
     xlab="Path",ylab="Log Acoustic Impulse",
     type="n",axes=F,frame.plot=TRUE)
axis(1,at=1:length(p_names),
     labels=lab_p,
     las=2,cex.axis=0.5)
axis(2)
for( ii in 1:len_p ){
  if( lab_p[ii] %in% p_names ){ col_p = "blue"
  } else { col_p = "cyan" }
  segments(ii,lb_p[ii],ii,ub_p[ii],lwd=2,col=col_p)
}
abline(h = 0,lty=2,lwd=2)
graphics.off()
