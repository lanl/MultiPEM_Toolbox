########################################################################
#                                                                      #
# This file contains code for Sequential Monte Carlo (SMC) sampling of #
# posterior distributions. This implementation of SMC is designed to   #
# only be used for second stage new event device parameter inference.  #
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

require(abind)
require(ramcmc)
require(doFuture)
require(iterators)

# effective sample size reduction factor
fEss = 0.95

# adaptive MCMC parameters
gamma0 = 2/3
alphat = 0.234

# effective SS for adaptation
beta = 0.5
ness = (1-beta)^(-1/gamma0)-1

# initial sampling from the prior distribution
# (default: uniform distribution on the hypercube)
unrestricted = function(N, range, psamp) { # initial sample
  if( !is.null(range) ){ D = nrow(range)
  } else { D = nrow(psamp) }
  samp = NULL
  if( !is.null(range) ){
    for (d in 1:D) samp = rbind(samp, runif(N, range[d,1], range[d,2]))
  } else {
    samp = psamp
    L = ncol(samp)
    if( N < L ){
      ind = sample(L,N)
      samp = samp[,ind]
    }
  }
  return(samp)
}

# log-prior distribution
# (default: uniform distribution)
log_P = function(x, p_par) {
  return(p_par$lprior_0(x, p_par))
}

# log-likelihood function
log_L = function(x, l_par) {
  return(l_par$ll_0(x, l_par))
}

# sampling at each time step
Gibbs = function(x, D, p_par, l_par, nu_t, gamma, alphat) {
  Lsigma = matrix(x[(D+3):(D*D+D+2)], nrow=D)
  U = rnorm(D,0,1)
  Z = Lsigma %*% U
  newx = x[1:D] + Z
  Z = Z / norm(as.matrix(U), type="2")
  lp = log_P(newx, p_par)
  ll = log_L(newx, l_par)
  ratio = lp - x[D+1] + nu_t * (ll - x[D+2])
  alpha = min(1, exp(ratio))
  a = 0
  u = runif(1)
  if (log(u) <= ratio) {
    x[1:(D+2)] = c(newx, lp, ll)
    a = 1
  }
  diffalpha = alpha - alphat
  fact = sqrt(gamma * abs(diffalpha))
  if (diffalpha >= 0){
    Lsigma = chol_update(Lsigma, fact*Z)
  } else {
    Lsigma = chol_downdate(Lsigma, fact*Z)
  }
  x[(D+3):(D*D+D+2)] = as.vector(Lsigma)
  return(list(x = x, a = a))
}

# weight calculator
iw = function(x, D, nu, nu0) {
  if (nu == nu0) return(0)
  else return((nu - nu0) * x[D+1])
}

# adaptive specification of the constraint parameter
adapt_seq = function(nu, nu0, D, N, ref, r_ess) {
  wt = foreach(x=iter(ref$sample, by='column'),
         .combine='c') %dofuture% {
         iw(x, D, nu, nu0)
       } %seed% TRUE
  Wt = ref$Wt * exp(wt)
  Wt = Wt / sum(Wt)
  ESS = ifelse(sum(is.na(Wt)) == N, 0, 1 / sum(Wt ^ 2))
  return(ESS - r_ess * N)
}

# list to matrix function
ltomat = function(y) {
  return(c(y$x,y$a))
}

# reference classes
rc1 = setRefClass("rc1", fields = list(sample = "matrix",
                                       Wt = "vector"))

# main function: sampling from a log-posterior distribution
SMC = function(N=1000, M=10, nuseq_T=1, range=NULL, psamp=NULL,
               pc=p_cal) {
  if( !is.null(range) ){ D = nrow(range)
  } else { D = nrow(psamp) }
  t = 1
  nuseq = c(0)
  a = c(0)
  Wt = array(dim=c(N, 1))
  samplet = unrestricted(N, range=range, psamp=psamp)
  lpdent = foreach(x=iter(samplet, by='column'),
             .combine='cbind') %dofuture% {
             c(log_P(x, p_par=pc), log_L(x, l_par=pc))
           } %seed% TRUE
  samplet = array(samplet, dim=c(D, 1, N))
  lpdent = array(lpdent, dim=c(2, 1, N))
  Wt[,1] = rep(1 / N, N)
  repeat {
    t = t+1
    newsample = rbind(samplet[,t-1,], lpdent[2,t-1,])
    newWt = Wt[,t-1]
    newref = rc1(sample = newsample, Wt = newWt)
    nuseq[t] = ifelse(adapt_seq(nu = nuseq_T, nu0 = nuseq[t-1], D = D,
                                N = N, ref = newref, r_ess = fEss) > 0,
                      nuseq_T,
                      uniroot(adapt_seq,
                              interval = c(nuseq[t-1], nuseq_T),
                              nu0 = nuseq[t-1], D = D, N = N,
                              ref = newref, r_ess = fEss)$root)
    if (nuseq[t] == nuseq[t-1]) wt = rep(0, N)
    else wt = (nuseq[t] - nuseq[t-1]) * newsample[D+1,]
    newsample = rbind(samplet[,t-1,], lpdent[,t-1,])
    newWt = newWt * exp(wt)
    newWt = newWt / sum(newWt)
    index = sample(1:N, N, prob = newWt, replace = T)
    newsample = newsample[,index]
    newWt = rep(1/N, N)
    Lsigma = t(chol(cov(t(newsample[1:D,,drop=FALSE]))))
    vLsigma = as.vector(Lsigma)
    newsample = rbind(newsample, matrix(rep(vLsigma,N),nrow=D*D))
    if (t > 2) a = c(a,0)
    for(i in 1:M) {
      gamma = min(1, D*((i+ness)^(-gamma0)))
      out = foreach(x=iter(newsample, by = 'column')) %dofuture% {
              Gibbs(x, D, pc, pc, nuseq[t], gamma, alphat)
            } %seed% TRUE
      mout = sapply(out,ltomat)
      newsample = mout[1:(D*D+D+2),]
      a[t-1] = a[t-1] + sum(mout[D*D+D+3,])
    }
    samplet = abind(samplet, newsample[1:D,,drop=FALSE], along = 2)
    lpdent = abind(lpdent,newsample[(D+1):(D+2),],along = 2)
    Wt = abind(Wt, newWt, along = 2)
    a[t-1] = a[t-1] / (M*N)
    if (nuseq[t] >= nuseq_T) break
  }
  t_final = t
  sample = matrix(samplet[,t_final,],nrow=D)
  return(list(sample=t(sample),acr=a))
}
