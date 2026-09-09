########################################################################
#                                                                      #
# Fused log-likelihood value+gradient for the calibration parameters.  #
#                                                                      #
# In NUTS, every leapfrog leaf evaluates BOTH the log-density and its  #
# gradient at the SAME parameter point. Computed separately (ll_full   #
# then gll_full) they each rebuild the identical residuals, model      #
# covariance matrices Omega, and Cholesky factor chol(Omega) -- the    #
# dominant per-leaf cost. The gradient routine (gll_cal) already forms  #
# chol(Omega), its inverse, and resid_io = Omega^{-1} resid; the        #
# log-likelihood is then a near-free byproduct                          #
#   ll = -sum(log(diag(cOmega))) - resid . resid_io / 2                 #
#        - n_Omega*log(2*pi)/2                                          #
# (see log_likelihood_cal.r). Fusing the two returns list(ll=, grad=)   #
# from ONE pass, ~halving the likelihood cost per leaf.                 #
#                                                                      #
# llg_cal is built by splicing the LIVE glog_likelihood_cal.r source at #
# load time (make_llg_cal), so the fused gradient can never silently    #
# drift from the production gradient. If the expected anchors are ever  #
# not found (e.g. gll_cal is restructured), make_llg_cal returns NULL   #
# and callers fall back to the separate ll_full/gll_full path.          #
#                                                                      #
# The fused ll reproduces ll_cal's EXACT accept/reject guards           #
# (infinite Omega, non-SPD Omega, kappa(cOmega) > 1000 -> logp = -Inf), #
# so a fused NUTS run is bitwise identical to a separate-callback run    #
# with the same seed -- fusing changes only HOW logp/grad are obtained. #
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

# Build the fused llg_cal by transforming the live glog_likelihood_cal.r text.
# Returns the function on success, or NULL if any anchor fails to match
# (signalling the caller to fall back to the non-fused path).
make_llg_cal <- function(gdir)
{
  src <- tryCatch(
    readLines(paste(gdir, "/glog_likelihood_cal.r", sep = ""), warn = FALSE),
    error = function(e) NULL)
  if (is.null(src)) return(NULL)
  txt <- paste(src, collapse = "\n")

  # Splice helper: the anchor MUST occur exactly once, else we bail (NULL).
  ok <- TRUE
  splice <- function(txt, anchor, replacement) {
    if (!ok) return(txt)
    hits <- gregexpr(anchor, txt, fixed = TRUE)[[1]]
    if (length(hits) != 1L || hits[1] == -1L) { ok <<- FALSE; return(txt) }
    sub(anchor, replacement, txt, fixed = TRUE)
  }

  # (A) rename the function
  txt <- splice(txt, "gll_cal = function(x, pc)", "llg_cal = function(x, pc)")

  # (B) initialize the log-likelihood accumulator. NOTE: `ll` is already used as
  # a loop index in the gradient body, so the accumulator must be named .ll_acc.
  txt <- splice(txt,
    "  # observational error covariance parameters\n  gr_eps = NULL\n",
    paste0("  # observational error covariance parameters\n  gr_eps = NULL\n",
           "  # log-likelihood accumulator (built from the same cOmega/resid_io\n",
           "  # the gradient computes). Named .ll_acc to avoid collision with the\n",
           "  # `ll` loop-index used in the variance-component loops.\n",
           "  .ll_acc = 0\n"))

  # (C) Replace the gradient's unguarded Cholesky with ll_cal's EXACT guarded
  # version and accumulate the log-likelihood from the same factor. Returning
  # logp = -Inf at precisely ll_cal's reject points keeps a fused NUTS run
  # bitwise identical to a separate-callback run. Where logp = -Inf the leaf is
  # divergent and its gradient is never reused, so grad = NaN there is safe.
  # The gradient source now factors Omega densely (as.matrix + base chol);
  # match that anchor and keep the dense chol here, adding ll_cal's guards.
  txt <- splice(txt,
    paste0(
    "      Omega = as.matrix(Omega)\n",
    "      cOmega = chol(Omega)\n",
    "      IOmega = chol2inv(cOmega)\n",
    "      resid_io = as.numeric(IOmega %*% resid)\n"),
    paste0(
    "      if( any(is.infinite(Omega)) ){ return(list(ll=-Inf, grad=NaN)) }\n",
    "      Omega = as.matrix(Omega)\n",
    "      cOmegaCatch = pc$tryCatch.W.E(chol(Omega))\n",
    "      if( is.matrix(cOmegaCatch$value) ){ cOmega = cOmegaCatch$value\n",
    "      } else { return(list(ll=-Inf, grad=NaN)) }\n",
    "      if( kappa(cOmega) > 1000 ){ return(list(ll=-Inf, grad=NaN)) }\n",
    "      IOmega = chol2inv(cOmega)\n",
    "      resid_io = as.numeric(IOmega %*% resid)\n",
    "      .ll_acc = .ll_acc - sum(log(diag(cOmega))) -\n",
    "           as.numeric(crossprod(resid, resid_io))/2 - n_Omega*log(2*pi)/2\n"))

  # (D) NaN forward-model guard -> list with logp -Inf (matches ll_cal's -Inf).
  txt <- splice(txt,
    "            if( any(is.nan(yhat)) ){ return(rep(NaN,pc$nmpars)) }",
    "            if( any(is.nan(yhat)) ){ return(list(ll=-Inf, grad=NaN)) }")

  # (E) final return -> list(ll=, grad=)
  txt <- splice(txt,
    "  return(c(gr_th0,gr_cp,gr_eiv,gr_beta0,gr_betat,gr_vc1,gr_vc2,gr_eps))",
    paste0("  return(list(ll=as.numeric(.ll_acc),\n",
           "              grad=c(gr_th0,gr_cp,gr_eiv,gr_beta0,gr_betat,",
           "gr_vc1,gr_vc2,gr_eps)))"))

  if (!ok) return(NULL)
  expr <- tryCatch(parse(text = txt), error = function(e) NULL)
  if (is.null(expr)) return(NULL)
  env <- new.env(parent = environment(make_llg_cal))
  eval(expr, envir = env)
  get("llg_cal", envir = env)
}

# Fused full log-likelihood value+gradient. Mirrors ll_full + gll_full exactly
# (same EIV index arithmetic) but calls the fused llg_cal once. list(ll=, grad=).
llg_full <- function(x, pc = p_cal)
{
  r   <- pc$llg_cal(x, pc)          # list(ll=, grad=)
  ll  <- r$ll
  gll <- r$grad

  if (exists("eiv", where = pc, inherits = FALSE) && pc$eiv) {
    # gradient EIV block -- identical index arithmetic to gll_full
    xg <- x
    st_geiv <- 0
    if (exists("nev", where = pc, inherits = FALSE) && pc$nev) {
      xg <- xg[-(1:pc$ntheta0)]
      st_geiv <- st_geiv + pc$ntheta0
    }
    if (pc$ncalp > 0) {
      st_geiv <- st_geiv + pc$ncalp
      xg <- xg[-(1:pc$ncalp)]
    }
    xg <- xg[1:pc$nsource]
    igeiv <- st_geiv + (1:pc$nsource)
    gll[igeiv] <- gll[igeiv] + pc$gll_eiv(xg, pc)

    # value EIV block -- identical extraction to ll_full
    xl <- x
    if (exists("nev", where = pc, inherits = FALSE) && pc$nev) {
      xl <- xl[-(1:pc$ntheta0)]
    }
    if (pc$ncalp > 0) {
      xl <- xl[-(1:pc$ncalp)]
    }
    xl <- xl[1:pc$nsource]
    ll <- ll + pc$ll_eiv(xl, pc)
  }

  list(ll = ll, grad = gll)
}

# Fused full log-POSTERIOR value+gradient. Mirrors lpost_full + glpost_full.
# Returns list(logp=, grad=) -- the exact contract nuts_step_fused_callback
# expects. The gradient prior padding matches glog_posterior_full.r.
lpg_full <- function(x, pc = p_cal)
{
  r   <- pc$llg_full(x, pc)
  ll  <- r$ll
  gll <- r$grad
  g   <- c(gll, numeric(pc$p_fgsn + pc$p_A)) + pc$glprior(x, pc)
  lp  <- ll + pc$lprior(x, pc)
  list(logp = lp, grad = g)
}
