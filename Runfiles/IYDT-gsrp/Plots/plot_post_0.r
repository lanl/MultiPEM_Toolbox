########################################################################
#                                                                      #
# This file plots posterior summaries of yield and height-of-burst     #
# based on various analyses.                                           #
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

samp_s = as.matrix(read.table(
         "../Seismic/I-SUGAR-hob-0/postsamp_seismic_0.txt"))
samp_a = as.matrix(read.table(
       "../Acoustic/I-SUGAR-hob-0/postsamp_acoustic_0.txt"))
samp_o = as.matrix(read.table(
         "../Optical/I-SUGAR-hob-0/postsamp_optical_0.txt"))
samp_o_eiv = as.matrix(read.table(
             paste("../Optical/I-EIV-SUGAR-hob-0/",
                   "postsamp_optical_eiv_0.txt",sep="")))
samp_c = as.matrix(read.table(
               "../Crater/I-SUGAR-0/postsamp_crater_0.txt"))
samp_c_eiv = as.matrix(read.table(
   "../Crater/I-EIV-SUGAR-0/postsamp_crater_eiv_0.txt"))
samp_4 = as.matrix(read.table(
            "../4-Phen/I-SUGAR-hob-0/postsamp_4phen_0.txt"))
samp_4_eiv = as.matrix(read.table(
"../4-Phen/I-EIV-SUGAR-hob-0/postsamp_4phen_eiv_0.txt"))

w_lb = min(samp_s[,1],samp_a[,1],samp_o[,1],samp_c,samp_4[,1])
w_ub = max(samp_s[,1],samp_a[,1],samp_o[,1],samp_c,samp_4[,1])

hob_lb = -10
hob_ub = 160

ijoint <- seq(1,nrow(samp_4),length=1000)
samp_s_y = samp_s[ijoint,]; samp_s_y[,1] = exp(samp_s_y[,1])/10^6
samp_a_y = samp_a[ijoint,]; samp_a_y[,1] = exp(samp_a_y[,1])/10^6
samp_o_y = samp_o[ijoint,]; samp_o_y[,1] = exp(samp_o_y[,1])/10^6
samp_o_eiv_y = samp_o_eiv[ijoint,]
samp_o_eiv_y[,1] = exp(samp_o_eiv_y[,1])/10^6
samp_c_y = samp_c[ijoint,]; samp_c_y = exp(samp_c_y)/10^6
samp_c_eiv_y = samp_c_eiv[ijoint,]
samp_c_eiv_y = exp(samp_c_eiv_y)/10^6
samp_4_y = samp_4[ijoint,]; samp_4_y[,1] = exp(samp_4_y[,1])/10^6
samp_4_eiv_y = samp_4_eiv[ijoint,]
samp_4_eiv_y[,1] = exp(samp_4_eiv_y[,1])/10^6

ss_bayes <- NULL
ss_bayes <- c(apply(samp_s_y,2,mean), 2*apply(samp_s[ijoint,],2,sd),
              cor(samp_s[ijoint,1],samp_s[ijoint,2]))
ss_bayes <- rbind(ss_bayes, c(apply(samp_a_y,2,mean),
                  2*apply(samp_a[ijoint,],2,sd),
                  cor(samp_a[ijoint,1],samp_a[ijoint,2])))
ss_bayes <- rbind(ss_bayes, c(apply(samp_o_y,2,mean),
                  2*apply(samp_o[ijoint,],2,sd),
                  cor(samp_o[ijoint,1],samp_o[ijoint,2])))
ss_bayes <- rbind(ss_bayes, c(mean(samp_c_y),NA,
                  2*sd(samp_c[ijoint,]),NA,NA))
ss_bayes <- rbind(ss_bayes, c(apply(samp_4_y,2,mean),
                  2*apply(samp_4[ijoint,],2,sd),
                  cor(samp_4[ijoint,1],samp_4[ijoint,2])))
ss_bayes[,3] <- exp(ss_bayes[,3])
rownames(ss_bayes) <- c("seismic","acoustic","optical","crater",
                        "multipem")
colnames(ss_bayes) <- c("mean_yield","mean_hob","re_yield","e_hob",
                        "corr")

ss_bayes_eiv <- NULL
ss_bayes_eiv <- c(apply(samp_o_eiv_y,2,mean),
                  2*apply(samp_o_eiv[ijoint,],2,sd),
                  cor(samp_o_eiv[ijoint,1],samp_o_eiv[ijoint,2]))
ss_bayes_eiv <- rbind(ss_bayes_eiv, c(mean(samp_c_eiv_y),NA,
                      2*sd(samp_c_eiv[ijoint,]),NA,NA))
ss_bayes_eiv <- rbind(ss_bayes_eiv, c(apply(samp_4_eiv_y,2,mean),
                      2*apply(samp_4_eiv[ijoint,],2,sd),
                      cor(samp_4_eiv[ijoint,1],samp_4_eiv[ijoint,2])))
ss_bayes_eiv[,3] <- exp(ss_bayes_eiv[,3])
rownames(ss_bayes_eiv) <- c("optical","crater","multipem")
colnames(ss_bayes_eiv) <- c("mean_yield","mean_hob","re_yield","e_hob",
                            "corr")

lY = seq(10,18,2)
Y = round(exp(lY)/10^6,1)
pdf("post_th0_0.pdf")
par(mgp=c(2,1,0),pty="s")
plot(samp_s[ijoint,1],samp_s[ijoint,2],col="red",pch=24,
     xlab="Yield (kt)",ylab="HOB/DOB (m)",
     xlim=c(9,19),ylim=c(hob_lb,hob_ub),
     cex=0.5,bg="red",axes=FALSE)
box()
axis(1,at=lY,labels=as.character(Y))
axis(2)
points(samp_a[ijoint,1],samp_a[ijoint,2],col="blue",pch=25,
       cex=0.5,bg="blue")
points(samp_o[ijoint,1],samp_o[ijoint,2],col="green",pch=22,
       cex=0.5,bg="green")
points(samp_4[ijoint,1],samp_4[ijoint,2],pch=20)
points(13.998,1.0668,pch=20,cex=2,col="orange")
abline(h = 0, lwd=2, lty=2, col="violet")
legend("topleft",legend=c("seismic","acoustic",
       "optical","multiPEM","sugar"),
       col=c("red","blue","green","black","orange"),
       pch=c(24,25,22,20,20),
       pt.bg=c("red","blue","green","black","orange"),
       pt.cex=c(0.5,0.5,0.5,1,1))
graphics.off()

y_s = exp(samp_s[,1])/10^6
y_a = exp(samp_a[,1])/10^6
y_o = exp(samp_o[,1])/10^6
y_c = exp(samp_c)/10^6
y_4 = exp(samp_4[,1])/10^6

y_lb = min(y_s,y_a,y_o,y_4)
y_ub = max(y_s,y_a,y_o,y_4)

pdf("post_th0_ykt_0.pdf")
par(mgp=c(2,1,0),pty="s")
plot(y_s[ijoint],samp_s[ijoint,2],col="red",pch=24,
     xlab="Yield (kt)",ylab="HOB/DOB (m)",
     xlim=c(0.745,3.75),ylim=c(hob_lb,hob_ub),
     cex=0.5,bg="red")
points(y_a[ijoint],samp_a[ijoint,2],col="blue",pch=25,
       cex=0.5,bg="blue")
points(y_o[ijoint],samp_o[ijoint,2],col="green",pch=22,
       cex=0.5,bg="green")
points(y_4[ijoint],samp_4[ijoint,2],pch=20)
points(1.2,1.0668,pch=20,col="orange")
abline(h = 0, lwd=2, lty=2, col="violet")
legend("topleft",legend=c("seismic","acoustic",
       "optical","multiPEM","sugar"),
       col=c("red","blue","green","black","orange"),
       pch=c(24,25,22,20,20),
       pt.bg=c("red","blue","green","black","orange"),
       pt.cex=c(0.5,0.5,0.5,1,1))
graphics.off()

w_s_den = density(samp_s[,1],adjust=2,from=min(samp_s[,1]),
                  to=max(samp_s[,1]))
w_a_den = density(samp_a[,1],adjust=2,from=min(samp_a[,1]),
                  to=max(samp_a[,1]))
w_o_den = density(samp_o[,1],adjust=2,from=min(samp_o[,1]),
                  to=max(samp_o[,1]))
w_c_den = density(samp_c,adjust=2,from=min(samp_c),
                  to=max(samp_c))
w_4_den = density(samp_4[,1],adjust=2,from=min(samp_4[,1]),
                  to=max(samp_4[,1]))

w_lb = min(w_s_den$x,w_a_den$x,w_o_den$x,w_c_den$x,w_4_den$x)
w_ub = max(w_s_den$x,w_a_den$x,w_o_den$x,w_c_den$x,w_4_den$x)
w_den_ub = max(w_s_den$y,w_a_den$y,w_o_den$y,w_c_den$y,w_4_den$y)

pdf("post_w_0.pdf")
par(mgp=c(2,1,0),pty="s")
plot(w_s_den$x,w_s_den$y,type="l",col="red",
     xlab="log(Yield (kt))",ylab="density",
     xlim=c(11,17),ylim=c(0,w_den_ub),lwd=2)
lines(w_a_den$x,w_a_den$y,col="blue",lwd=2)
lines(w_o_den$x,w_o_den$y,col="green",lwd=2)
lines(w_c_den$x,w_c_den$y,col="violet",lwd=2)
lines(w_4_den$x,w_4_den$y,lwd=2)
abline(h=0,col="gray")
abline(v=13.998,col="orange",lwd=2,lty=2)
legend("topright",legend=c("seismic","acoustic",
       "optical","surface","multiPEM","sugar"),
       col=c("red","blue","green","violet",
             "black","orange"),lty=c(1,1,1,1,1,2))
graphics.off()

y_max = 3

if(is.null(y_max)){ ul = max(y_s) } else { ul = y_max }
y_s_den = density(y_s,adjust=2,from=min(y_s),to=ul)
if(is.null(y_max)){ ul = max(y_a) } else { ul = y_max }
y_a_den = density(y_a,adjust=2,from=min(y_a),to=ul)
if(is.null(y_max)){ ul = max(y_o) } else { ul = y_max }
y_o_den = density(y_o,adjust=2,from=min(y_o),to=ul)
if(is.null(y_max)){ ul = max(y_c) } else { ul = y_max }
y_c_den = density(y_c,adjust=2,from=min(y_c),to=ul)
if(is.null(y_max)){ ul = max(y_4) } else { ul = y_max }
y_4_den = density(y_4,adjust=2,from=min(y_4),to=ul)

y_lb = min(y_s_den$x,y_a_den$x,y_o_den$x,y_c_den$x,y_4_den$x)
y_ub = max(y_s_den$x,y_a_den$x,y_o_den$x,y_c_den$x,y_4_den$x)
y_den_ub = max(y_s_den$y,y_a_den$y,y_o_den$y,y_c_den$y,y_4_den$y)

pdf("post_y_0.pdf")
par(mgp=c(2,1,0),pty="s")
plot(y_s_den$x,y_s_den$y,type="l",col="red",
     xlab="Yield (kt)",ylab="density",
     xlim=c(y_lb,3),ylim=c(0,y_den_ub),lwd=2)
lines(y_a_den$x,y_a_den$y,col="blue",lwd=2)
lines(y_o_den$x,y_o_den$y,col="green",lwd=2)
lines(y_c_den$x,y_c_den$y,col="violet",lwd=2)
lines(y_4_den$x,y_4_den$y,lwd=2)
abline(h=0,col="gray")
abline(v=1.2,col="orange",lwd=2,lty=2)
legend("topright",legend=c("seismic","acoustic",
       "optical","surface","multiPEM","sugar"),
       col=c("red","blue","green","violet",
             "black","orange"),lty=c(1,1,1,1,1,2))
graphics.off()

hob_s_den = density(samp_s[,2],adjust=5,from=min(samp_s[,2]),
                    to=max(samp_s[,2]))
hob_a_den = density(samp_a[,2],adjust=5,from=min(samp_a[,2]),
                    to=max(samp_a[,2]))
hob_o_den = density(samp_o[,2],adjust=5,from=min(samp_o[,2]),
                    to=max(samp_o[,2]))
hob_4_den = density(samp_4[,2],adjust=5,from=min(samp_4[,2]),
                    to=max(samp_4[,2]))

hob_lb = min(hob_s_den$x,hob_a_den$x,hob_o_den$x,hob_4_den$x)
hob_ub = max(hob_s_den$x,hob_a_den$x,hob_o_den$x,hob_4_den$x)
hob_den_ub = max(hob_s_den$y,hob_a_den$y,hob_o_den$y,hob_4_den$y)

pdf("post_hob_0.pdf")
par(mgp=c(2,1,0),pty="s")
plot(hob_s_den$x,hob_s_den$y,type="l",col="red",
     xlab="HOB/DOB (m)",ylab="density",
     xlim=c(hob_lb,hob_ub),ylim=c(0,hob_den_ub),lwd=2)
lines(hob_a_den$x,hob_a_den$y,col="blue",lwd=2)
lines(hob_o_den$x,hob_o_den$y,col="green",lwd=2)
lines(hob_4_den$x,hob_4_den$y,lwd=2)
abline(h=0,col="gray")
abline(v=1.0668,col="orange",lwd=2,lty=2)
legend("topright",legend=c("seismic","acoustic",
       "optical","multiPEM","sugar"),
       col=c("red","blue","green",
             "black","orange"),lty=c(1,1,1,1,2))
graphics.off()
