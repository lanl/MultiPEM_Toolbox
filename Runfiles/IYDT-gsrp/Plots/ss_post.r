########################################################################
#                                                                      #
# This file prints the posterior means of yield and height-of-burst,   #
# their joint posterior covariance matrix, and 95% credible intervals  #
# for new event yield and height-of-burst.                             #
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

# Posterior Summaries
post = p_cal$tmpi_0

pmu_0 = apply(post,2,mean)
pmu_0[1] = exp(pmu_0[1])/10^6
print(pmu_0)

ps_0 = 2*apply(post,2,sd)
ps_0[1] = exp(ps_0[1])
print(ps_0)
print(signif((ps_0[1]-1)*100,2))

if( length(pmu_0) > 1 ){
  pCorr = cor(post)
  print(pCorr)
}

# Credible Intervals
post[,1] = exp(post[,1])/10^6
lb_p = apply(post,2,quantile,probs=0.025)
ub_p = apply(post,2,quantile,probs=0.975)
print(lb_p)
print(ub_p)
