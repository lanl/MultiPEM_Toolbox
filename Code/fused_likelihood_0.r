########################################################################
#                                                                      #
# Fused log-likelihood value+gradient for the NEW EVENT parameters      #
# (the "_0" stack used by calc_bayes_0.r).                              #
#                                                                      #
# As in fused_likelihood_cal.r, every NUTS leapfrog leaf needs BOTH the #
# log-density and its gradient at the same point. gll_0 already forms    #
# the residuals and IOmega %*% resid needed by the log-likelihood        #
#   ll = -logdet_cOmega - resid . (IOmega resid)/2 - n_h0_tot*log(2pi)/2 #
# (see log_likelihood_0.r); here the model covariance Omega is FIXED     #
# (IOmega and logdet_cOmega are precomputed at calibration), so there is #
# no per-call Cholesky -- fusing simply avoids the duplicate forward-    #
# model + residual pass that a separate ll_0 call would repeat.          #
#                                                                      #
# llg_0 is built by splicing the LIVE glog_likelihood_0.r source at      #
# load time (make_llg_0), so the fused gradient cannot drift from the    #
# production gradient. If the expected anchors are ever not found,       #
# make_llg_0 returns NULL and callers fall back to the separate          #
# lpost_0/glpost_0 path.                                                #
#                                                                      #
# The fused ll reproduces ll_0's exact NaN forward-model guard           #
# (yhat NaN -> logp = -Inf), so a fused NUTS run is bitwise identical to  #
# a separate-callback run with the same seed.                            #
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

# Build the fused llg_0 by transforming the live glog_likelihood_0.r text.
# Returns the function on success, or NULL if any anchor fails to match.
make_llg_0 <- function(gdir)
{
  src <- tryCatch(
    readLines(paste(gdir, "/glog_likelihood_0.r", sep = ""), warn = FALSE),
    error = function(e) NULL)
  if (is.null(src)) return(NULL)
  txt <- paste(src, collapse = "\n")

  ok <- TRUE
  splice <- function(txt, anchor, replacement) {
    if (!ok) return(txt)
    hits <- gregexpr(anchor, txt, fixed = TRUE)[[1]]
    if (length(hits) != 1L || hits[1] == -1L) { ok <<- FALSE; return(txt) }
    sub(anchor, replacement, txt, fixed = TRUE)
  }

  # (A) rename the function
  txt <- splice(txt, "gll_0 = function(x, pc)", "llg_0 = function(x, pc)")

  # (B) initialize the log-likelihood accumulator alongside the gradient init.
  txt <- splice(txt,
    "  # new event parameters\n  gr_th0 = numeric(pc$ntheta0)\n",
    paste0("  # new event parameters\n  gr_th0 = numeric(pc$ntheta0)\n",
           "  # log-likelihood accumulator (built from the same residuals /\n",
           "  # IOmega the gradient uses; see log_likelihood_0.r).\n",
           "  .ll_acc = 0\n"))

  # (C) NaN forward-model guard -> list with logp -Inf (matches ll_0's -Inf).
  txt <- splice(txt,
    "        if( any(is.nan(yhat)) ){ return(rep(NaN,pc$ntheta0)) }",
    "        if( any(is.nan(yhat)) ){ return(list(ll=-Inf, grad=NaN)) }")

  # (D) Accumulate the log-likelihood using ll_0's exact arithmetic expression
  # and leave the gradient's g_th0 line unchanged. Using the same
  # t(resid) %*% IOmega %*% resid expression as ll_0 (rather than factoring
  # IOmega %*% resid out of the g_th0 product) makes the fused log-density agree
  # bit-for-bit with the separate ll_0. ll_0's contribution
  # (log_likelihood_0.r):
  #   ll = ll - logdet_cOmega - t(resid)%*%IOmega%*%resid/2 - n_h0_tot*log(2pi)/2
  txt <- splice(txt,
    "    g_th0 = t(Jac_th0) %*% pc$h[[hh]]$IOmega %*% resid\n",
    paste0(
    "    .ll_acc = .ll_acc - pc$h[[hh]]$logdet_cOmega -\n",
    "         as.numeric(t(resid) %*% pc$h[[hh]]$IOmega %*% resid)/2 -\n",
    "         n_h0_tot*log(2*pi)/2\n",
    "    g_th0 = t(Jac_th0) %*% pc$h[[hh]]$IOmega %*% resid\n"))

  # (E) final return -> list(ll=, grad=)
  txt <- splice(txt,
    "  return(as.vector(gr_th0))",
    "  return(list(ll=as.numeric(.ll_acc), grad=as.vector(gr_th0)))")

  if (!ok) return(NULL)
  expr <- tryCatch(parse(text = txt), error = function(e) NULL)
  if (is.null(expr)) return(NULL)
  env <- new.env(parent = environment(make_llg_0))
  eval(expr, envir = env)
  get("llg_0", envir = env)
}

# Fused full log-POSTERIOR value+gradient for the _0 stack. Mirrors
# lpost_0 (= ll_0 + lprior_0) and glpost_0 (= gll_0 + glprior_0).
# Returns list(logp=, grad=) -- the contract nuts_step_fused_callback expects.
lpg_0 <- function(x, pc = p_cal)
{
  r <- pc$llg_0(x, pc)             # list(ll=, grad=)
  list(logp = r$ll  + pc$lprior_0(x, pc),
       grad = r$grad + pc$glprior_0(x, pc))
}
