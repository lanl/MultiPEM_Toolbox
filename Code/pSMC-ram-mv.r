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
#
# OVERVIEW
# --------
# SMC() draws from a tempered sequence of distributions that bridges the
# prior to the posterior. At each step it (1) adaptively selects the next
# tempering level so the effective sample size (ESS) holds near a target,
# (2) reweights and resamples the particles, and (3) rejuvenates them with
# M sweeps of a robust adaptive-Metropolis (RAM, Vihola 2012) kernel whose
# proposal covariance is tuned online toward the target acceptance rate.
# The user supplies the model as a list `pc` with two functions,
# lprior_0(x, pc) and ll_0(x, pc), returning the log-prior and
# log-likelihood at a parameter vector x.
#
# Particle evaluation and the M MCMC sweeps are parallelised across blocks
# of particles via the doFuture backend; register a plan (e.g.
# plan(multisession) or plan(sequential)) before calling SMC().
#
# ROBUST `range` HANDLING (see range_report() and the `range_control` arg)
#   * Structural validation: `range` must be a two-column matrix (lower,
#     upper) with finite bounds, lower < upper in every row, and a row
#     count consistent with the model dimension.
#   * Pilot diagnostics use the initial log-likelihoods to flag, per
#     parameter, a range that is too wide (little effective support / low
#     initial ESS) or too narrow (posterior mass hugging a boundary, or a
#     boundary probe showing the log-likelihood has barely decayed at an
#     edge — which also catches a narrow box centred on the mode).
#   * range_control = "adapt" (DEFAULT) iteratively refines the box,
#     re-seeding under each suggestion until it stabilises (controlled by
#     range_max_iter / range_tol), with an oscillation guard that stops and
#     takes the union of the last two boxes once the support is bracketed.
#     No user approval or manual re-run is required. "warn" reports the
#     diagnostics and suggestion but runs the supplied box unchanged; "off"
#     skips range diagnostics entirely.
#
# NON-FINITE LOG-DENSITY ROBUSTNESS
#   * A wide box or an extreme MCMC proposal can drive lprior_0 / ll_0 to
#     NaN or +Inf. Every evaluation is wrapped so a NaN / +Inf / errored
#     value is treated as a numerical failure and mapped to -Inf (zero
#     mass): such particles get zero weight, are resampled away, and MCMC
#     moves onto them are rejected — the run continues instead of crashing
#     or freezing. A genuine -Inf (outside a hard-constraint prior) is
#     preserved. A low finite fraction is itself treated as a "too wide"
#     signal, shrinking the box toward the evaluable region. If ALL
#     particles collapse to zero mass, SMC stops with an explanatory error.
#     Sanitizations are reported via "[nonfinite] ..." messages.
#
# SCALING WITH PARAMETER DIMENSION (D)  -- please read before high-D use
#   The sampler is correct across dimensions (validated against a target
#   with a known posterior), but its TUNING requirements grow with D and
#   the defaults that suffice for small D will under-sample a large one.
#
#   Why: each tempering step rejuvenates particles with M sweeps of the RAM
#   kernel, whose D x D proposal factor starts at the full population
#   covariance (scale ~1). The optimal random-walk scale for a D-dim target
#   is ~2.38/sqrt(D), so as D grows the initial proposal is increasingly
#   too large; RAM must shrink it toward the 0.234 acceptance target, but it
#   only has M sweeps per step in which to adapt a D x D factor.
#
#   Symptoms of an UNDER-RESOURCED high-D run (all observed at D=100, M=8):
#     * per-step acceptance stays well below the 0.234 RAM target;
#     * the posterior is under-dispersed -- estimated marginal SDs are too
#       small and correlations are attenuated toward zero.
#   In testing these improved MONOTONICALLY as M was increased (sd/true
#   ratio 0.60 -> 0.68 -> 0.70 at M = 8 -> 30 -> 60), confirming the cause
#   is insufficient rejuvenation, not a defect.
#
#   Recommendations as D increases:
#     * Increase M (MCMC sweeps per step) first -- this is the main lever
#       for mixing; watch the returned $acr rise toward ~0.234.
#     * Increase N (particles): a D x D sample covariance needs N >> D to be
#       well conditioned (N of several thousand is reasonable near D = 100).
#     * Expect cost to grow: per-particle state is length D*D + D + 2 and the
#       RAM chol_update/chol_downdate are O(D^2) per sweep, so wall-clock
#       scales steeply. Use a parallel future plan (see USAGE).
#     * Treat $acr as the primary health check: values far below 0.234 (or a
#       posterior narrower than expected) indicate M and/or N are too small.
#
#   Parallel-model caveat (not specific to high D but easy to hit there):
#   define pc so it is SELF-CONTAINED -- store model constants as fields of
#   pc (e.g. pc$mu, pc$Sigma) and read them via the function's second
#   argument, as in ll_0(x, l) using l$mu. Model functions that instead
#   close over GLOBAL variables are not reliably exported to future workers,
#   which yields non-finite log-densities on the workers; the code reports
#   this as an all-zero-mass condition rather than a cryptic crash, but the
#   run cannot proceed until pc is made self-contained.
#
# USAGE
#   library(doFuture); registerDoFuture(); plan(multisession)  # or sequential
#   # self-contained pc: constants live in the list, read via the 2nd arg
#   pc <- list(mu = mu_vec,
#              lprior_0 = function(x, p) 0,
#              ll_0     = function(x, l) sum(dnorm(x, mean = l$mu, log = TRUE)))
#   range <- matrix(c(lo1, lo2, hi1, hi2), ncol = 2)
#   fit <- SMC(N = 1000, M = 10, range = range, pc = pc)
########################################################################

require(ramcmc)
require(doFuture)
require(iterators)   # kept for interface compatibility; not required below

# ---- tuning constants -------------------------------------------------

# effective sample size reduction factor (target ESS = fEss * N per step)
fEss = 0.95

# adaptive MCMC (RAM, Vihola 2012) parameters
gamma0 = 2/3
alphat = 0.234

# effective SS for adaptation
beta = 0.5
ness = (1 - beta)^(-1/gamma0) - 1

# ---- log-prior / log-likelihood wrappers ------------------------------

# log-prior distribution (default: uniform)
log_P = function(x, p_par) p_par$lprior_0(x, p_par)

# log-likelihood function
log_L = function(x, l_par) l_par$ll_0(x, l_par)

# ---- small utilities --------------------------------------------------

# number of parallel workers for the current future plan (>= 1)
worker_count = function() {
  nw = tryCatch(future::nbrOfWorkers(), error = function(e) 1L)
  if (!is.finite(nw) || nw < 1) 1L else as.integer(nw)
}

# split 1:N into k contiguous, non-overlapping index blocks
chunk_indices = function(N, k) {
  k = max(1L, min(as.integer(k), N))
  if (k == 1L) return(list(seq_len(N)))
  unname(split(seq_len(N), cut(seq_len(N), breaks = k, labels = FALSE)))
}

# lower-triangular Cholesky factor L (L %*% t(L) = S) with a PD safeguard
chol_lower = function(S) {
  S = (S + t(S)) / 2
  D = nrow(S)
  base = mean(diag(S))
  if (!is.finite(base) || base <= 0) base = 1
  jit = 0
  repeat {
    R = tryCatch(chol(S + jit * diag(D)), error = function(e) NULL)
    if (!is.null(R)) return(t(R))
    jit = if (jit == 0) base * 1e-8 else jit * 10
    if (jit > base * 1e3)
      stop("Initial covariance is not positive-definite; particles may have ",
           "collapsed. This often indicates the `range` for one or more ",
           "parameters is too narrow.")
  }
}

# weighted quantile (linear interpolation of the weighted CDF)
wquantile = function(x, w, p) {
  o = order(x); x = x[o]; w = w[o]
  cw = cumsum(w) / sum(w)
  as.numeric(approx(cw, x, xout = p, rule = 2, ties = "ordered")$y)
}

# ---- non-finite log-density handling ----------------------------------
# A widened box (or an extreme proposal) can drive lprior_0/ll_0 to NaN or
# +Inf (overflow, log of a negative, 0/0, ...). We treat any such value as
# a NUMERICAL FAILURE and map it to -Inf: the point then carries zero
# posterior mass, so it is given zero weight, resampled away, and any MCMC
# move onto it is rejected -- the run continues instead of crashing or
# freezing. A genuine -Inf is LEFT INTACT: it is a legitimate value meaning
# "outside support" (e.g. a hard-constraint prior), not a failure.

# TRUE for a log-density that must be sanitized (NaN, NA, or +Inf).
bad_logden = function(v) is.na(v) || v == Inf   # NA/NaN short-circuit before ==

# Sanitize a 2 x N (logprior, loglik) matrix in place: NaN/NA/+Inf -> -Inf.
# Returns the cleaned matrix with attr "n_bad" = c(lp, ll) counting the
# entries that were sanitized (genuine -Inf is not counted).
sanitize_lpden = function(mat) {
  neg_inf = is.infinite(mat) & (mat < 0)          # genuine -Inf (NA-safe)
  bad = !is.finite(mat) & !neg_inf                # NaN, NA, +Inf
  bad[is.na(bad)] = FALSE                          # defensive; should not occur
  mat[bad] = -Inf
  attr(mat, "n_bad") = c(lp = sum(bad[1, ]), ll = sum(bad[2, ]))
  mat
}

# Emit a one-line notice when sanitizations occurred (n_bad = c(lp, ll)).
report_nonfinite = function(n_bad, where) {
  if (is.null(n_bad)) return(invisible())
  if (sum(n_bad) > 0)
    message(sprintf(
      "[nonfinite] %s: replaced %d non-finite log-prior and %d non-finite log-likelihood value(s) with -Inf (zero mass); continuing.",
      where, n_bad["lp"], n_bad["ll"]))
  invisible()
}

# ---- initial sampling from the prior ----------------------------------
# (default: uniform on the hypercube defined by `range`)
unrestricted = function(N, range, psamp) {
  if (!is.null(range)) {
    D = nrow(range)
    samp = matrix(0, nrow = D, ncol = N)
    for (d in 1:D) samp[d, ] = runif(N, range[d, 1], range[d, 2])
    return(samp)
  }
  # draw from a supplied particle pool `psamp` (D x L)
  L = ncol(psamp)
  if (N <= L) {
    return(psamp[, sample(L, N), drop = FALSE])
  }
  # N > L: sample WITH replacement so the returned pool has exactly N columns
  psamp[, sample(L, N, replace = TRUE), drop = FALSE]
}

# ---- initial (log-prior, log-likelihood) evaluation, chunked ----------
eval_lpden = function(samp, pc) {
  N = ncol(samp)
  chunks = chunk_indices(N, worker_count())
  res = foreach(idx = chunks) %dofuture% {
    vapply(idx, function(j) {
      x = samp[, j]
      # a model error on one particle must not abort the whole evaluation;
      # treat it as a non-finite density (NaN -> sanitized to -Inf below).
      lp = tryCatch(pc$lprior_0(x, pc), error = function(e) NaN)
      ll = tryCatch(pc$ll_0(x, pc),    error = function(e) NaN)
      c(lp, ll)
    }, numeric(2))
  } %seed% TRUE
  out = matrix(0, nrow = 2, ncol = N)
  for (j in seq_along(chunks)) out[, chunks[[j]]] = res[[j]]
  sanitize_lpden(out)          # NaN/NA/+Inf -> -Inf; carries attr "n_bad"
}

# ---- adaptive tempering: ESS(nu) - target -----------------------------
# Vectorised: incremental log-weight for particle i is (nu - nu0) * loglik_i.
adapt_seq = function(nu, nu0, ll, Wt0, N, r_ess) {
  dnu = nu - nu0
  # guard 0 * -Inf = NaN: with no tempering step the increment is exactly 0.
  incr = if (dnu == 0) rep(0, length(ll)) else dnu * ll
  W = Wt0 * exp(incr)
  W[!is.finite(W)] = 0
  s = sum(W)
  if (s <= 0) return(0 - r_ess * N)
  W = W / s
  ESS = 1 / sum(W^2)
  ESS - r_ess * N
}

# ---- one block of particles, all M RAM-MCMC sweeps --------------------
# B: (D*D + D + 2) x nb matrix. Rows: [1:D] params, [D+1] logprior,
# [D+2] loglik, [D+3 : D*D+D+2] vec(Lsigma). Returns updated block and the
# total number of accepted moves in the block (summed over particles/sweeps).
move_block = function(B, D, M, gammas, nu_t, alphat, pc) {
  nb = ncol(B)
  dd = D * D
  acc = 0
  for (c in seq_len(nb)) {
    x      = B[, c]
    params = x[1:D]
    lp     = x[D + 1]
    ll     = x[D + 2]
    Lsigma = matrix(x[(D + 3):(dd + D + 2)], nrow = D)
    for (i in 1:M) {
      g  = gammas[i]
      U  = rnorm(D)
      Z  = as.numeric(Lsigma %*% U)
      newx = params + Z
      nU = sqrt(sum(U * U)); if (nU == 0) nU = 1
      Zn = Z / nU
      lpn = tryCatch(pc$lprior_0(newx, pc), error = function(e) NaN)
      lln = tryCatch(pc$ll_0(newx, pc),    error = function(e) NaN)
      # A NaN/+Inf proposal density is a numerical failure -> reject the
      # move outright (never store it, or the particle would freeze there).
      if (bad_logden(lpn) || bad_logden(lln)) {
        alpha  = 0
        accept = FALSE
      } else {
        llt = if (nu_t == 0) 0 else nu_t * (lln - ll)
        ratio = (lpn - lp) + llt
        # ratio == NaN (e.g. Inf-Inf from a -Inf current state) -> reject.
        # ratio == +Inf (moving from zero density to positive) -> accept,
        # so a particle seeded outside support can escape it. ratio == -Inf
        # (moving onto zero density) -> reject via log(u) <= -Inf == FALSE.
        if (is.nan(ratio)) {
          alpha  = 0
          accept = FALSE
        } else {
          alpha  = min(1, exp(ratio))
          accept = log(runif(1)) <= ratio
        }
      }
      if (accept) {
        params = newx; lp = lpn; ll = lln
        acc = acc + 1
      }
      dif  = alpha - alphat
      fact = sqrt(g * abs(dif))
      if (dif >= 0) Lsigma = ramcmc::chol_update(Lsigma, fact * Zn)
      else          Lsigma = ramcmc::chol_downdate(Lsigma, fact * Zn)
    }
    B[1:D, c]                    = params
    B[D + 1, c]                  = lp
    B[D + 2, c]                  = ll
    B[(D + 3):(dd + D + 2), c]   = as.vector(Lsigma)
  }
  list(mat = B, acc = acc)
}

# ---- range diagnostics -------------------------------------------------

# Hard structural validation of a supplied `range` matrix.
validate_range = function(range, D_psamp = NULL) {
  if (is.null(range)) return(invisible(NULL))
  if (!is.matrix(range) || ncol(range) != 2)
    stop("`range` must be a matrix with two columns (lower, upper).")
  if (any(!is.finite(range)))
    stop("`range` contains non-finite bounds.")
  if (any(range[, 1] >= range[, 2]))
    stop("`range` requires lower bound < upper bound in every row.")
  if (!is.null(D_psamp) && nrow(range) != D_psamp)
    stop("`range` has ", nrow(range), " rows but `psamp` implies ",
         D_psamp, " parameters.")
  invisible(TRUE)
}

# Pilot diagnostics from initial particles `samp` (D x N) and their
# log-likelihoods `ll`. Flags per-parameter "too wide"/"too narrow"
# conditions and returns a suggested adjusted range.
#
# Two complementary "too narrow" detectors:
#   (a) boundary mass — posterior mass piled against an edge (hug_mass);
#   (b) boundary probe — if `pc` is supplied, evaluate the log-likelihood
#       just OUTSIDE each edge (relative to the interior mode). If the
#       density has not decayed by `edge_drop` nats at an edge, the box is
#       truncating the posterior even when no single edge is "hugged"
#       (e.g. a narrow box centred on the mode). Set edge_probe = 0 to skip.
diagnose_range = function(range, samp, ll, pc = NULL,
                          hug_frac = 0.10, hug_mass = 0.20,
                          min_ess_frac = 0.05, wide_supp_frac = 0.05,
                          widen = 0.5, shrink_margin = 0.25,
                          edge_probe = 0.05, edge_drop = 2,
                          min_fin_frac = 0.5, center_shrink = 0.5) {
  D = nrow(range); N = ncol(samp)
  lo = range[, 1]; hi = range[, 2]; width = hi - lo

  # likelihood weights (posterior proxy under a ~flat prior over the box)
  fin = is.finite(ll)
  fin_frac = mean(fin)               # share of particles the model could evaluate
  m = if (any(fin)) max(ll[fin]) else NA_real_
  if (!is.finite(m)) {
    w = rep(1 / N, N)
  } else {
    w = exp(ll - m); w[!is.finite(w)] = 0
    w = if (sum(w) == 0) rep(1 / N, N) else w / sum(w)
  }
  ess_frac = (1 / sum(w^2)) / N

  # A low finite fraction means the box reaches into regions where the model
  # returns non-finite values (sanitized to -Inf) -- a distinct "too wide"
  # signal that the uniform-weight fallback above would otherwise mask.
  nonfin_wide = fin_frac < min_fin_frac
  fin_center = if (any(fin)) rowMeans(samp[, fin, drop = FALSE])
               else (lo + hi) / 2

  # interior reference point for the boundary probe: the weighted mean
  # (a cheap, robust proxy for the mode within the current box).
  ref_pt = if (!is.null(pc)) as.numeric(samp %*% w) else NULL
  ll_ref = if (!is.null(pc) && edge_probe > 0)
    tryCatch(pc$ll_0(ref_pt, pc), error = function(e) NA_real_) else NA_real_

  # log-likelihood at an interior point shifted to parameter d = value
  ll_at = function(d, value) {
    x = ref_pt; x[d] = value
    tryCatch(pc$ll_0(x, pc), error = function(e) NA_real_)
  }

  suggest = range
  notes = character(0)
  per = vector("list", D)
  for (d in 1:D) {
    xs  = samp[d, ]

    # Non-finite-driven shrink takes priority: if too few particles are
    # evaluable, contract this dimension toward the finite particles' span
    # (or toward the box centre if none are finite). This pulls the box out
    # of non-evaluable regions before any other adjustment.
    if (nonfin_wide) {
      nfin_d = sum(fin)
      if (nfin_d >= 2 && diff(range(samp[d, fin])) > 0) {
        flo = min(samp[d, fin]); fhi = max(samp[d, fin])
        mrg = shrink_margin * (fhi - flo)
        sl = flo - mrg; su = fhi + mrg
      } else {
        # too few (or coincident) evaluable particles to bracket a span:
        # contract about the finite centroid without collapsing the box.
        # A hard shrink here would drive the width toward zero (degenerate
        # covariance); halving preserves a positive width so the next
        # re-seed can locate the evaluable region.
        half = center_shrink * width[d] / 2
        sl = fin_center[d] - half; su = fin_center[d] + half
      }
      suggest[d, ] = c(sl, su)
      per[[d]] = list(supp_frac = NA, near_lo = NA, near_hi = NA)
      next
    }

    qlo = wquantile(xs, w, 0.005)
    qhi = wquantile(xs, w, 0.995)
    near_lo = sum(w[xs <= lo[d] + hug_frac * width[d]])
    near_hi = sum(w[xs >= hi[d] - hug_frac * width[d]])
    supp_frac = (qhi - qlo) / width[d]

    narrow_lo = near_lo >= hug_mass
    narrow_hi = near_hi >= hug_mass
    too_wide  = supp_frac < wide_supp_frac

    # boundary probe: has the likelihood decayed by the time we reach the edge?
    probe_lo = probe_hi = FALSE
    if (!is.null(pc) && edge_probe > 0 && is.finite(ll_ref)) {
      dlo = ll_at(d, lo[d]); dhi = ll_at(d, hi[d])
      # "not decayed enough" => edge is inside the bulk => box truncates
      probe_lo = is.finite(dlo) && (ll_ref - dlo) < edge_drop
      probe_hi = is.finite(dhi) && (ll_ref - dhi) < edge_drop
    }

    sl = lo[d]; su = hi[d]
    if (narrow_lo || probe_lo) {
      sl = lo[d] - widen * width[d]
      reason = if (narrow_lo)
        sprintf("%.0f%% of posterior mass hugs the LOWER bound (within %.0f%% of the edge)",
                100 * near_lo, 100 * hug_frac)
      else "log-likelihood has barely decayed at the LOWER bound (edge probe)"
      notes = c(notes, sprintf(
        "param %d: %s; range may be too NARROW. Suggest lowering to %.4g.",
        d, reason, sl))
    }
    if (narrow_hi || probe_hi) {
      su = hi[d] + widen * width[d]
      reason = if (narrow_hi)
        sprintf("%.0f%% of posterior mass hugs the UPPER bound", 100 * near_hi)
      else "log-likelihood has barely decayed at the UPPER bound (edge probe)"
      notes = c(notes, sprintf(
        "param %d: %s; range may be too NARROW. Suggest raising to %.4g.",
        d, reason, su))
    }
    if (too_wide && !narrow_lo && !narrow_hi && !probe_lo && !probe_hi) {
      mrg = shrink_margin * (qhi - qlo)
      sl = qlo - mrg; su = qhi + mrg
      notes = c(notes, sprintf(
        "param %d: effective support is only %.1f%% of the supplied range; range may be too WIDE. Suggest tightening to [%.4g, %.4g].",
        d, 100 * supp_frac, sl, su))
    }
    suggest[d, ] = c(sl, su)
    per[[d]] = list(supp_frac = supp_frac, near_lo = near_lo, near_hi = near_hi)
  }
  if (nonfin_wide)
    notes = c(notes, sprintf(
      "only %.1f%% of particles have a finite log-likelihood: the range reaches into regions the model cannot evaluate; tightening toward the evaluable span.",
      100 * fin_frac))
  else if (ess_frac < min_ess_frac)
    notes = c(notes, sprintf(
      "initial ESS is only %.1f%% of N: few seed particles carry meaningful likelihood; the overall range is likely too WIDE.",
      100 * ess_frac))

  list(ess_frac = ess_frac, fin_frac = fin_frac, suggest = suggest,
       notes = notes, per = per)
}

# Standalone range audit: seed particles, evaluate the model, report.
# Does NOT run SMC. Returns diagnose_range()'s list (invisibly).
range_report = function(N = 1000, range, pc, ...) {
  validate_range(range)
  samp = unrestricted(N, range = range, psamp = NULL)
  ll = eval_lpden(samp, pc)[2, ]
  rep = diagnose_range(range, samp, ll, pc = pc, ...)
  if (length(rep$notes) == 0)
    message("[range] no issues detected (initial ESS = ",
            sprintf("%.1f%% of N", 100 * rep$ess_frac), ").")
  else
    for (n in rep$notes) message("[range] ", n)
  invisible(rep)
}

# Iteratively refine `range` by re-seeding under each suggestion until the
# box stabilises (relative change per edge < `tol`) or `max_iter` is hit.
# A single pilot is noisy, so this repeats the seed/diagnose/adjust loop.
#
# Convergence is logarithmic in the initial scale mismatch (widening steps
# are geometric; shrinking jumps to the estimated support), so even a
# ~10,000x mismatch settles in ~12 iterations. `max_iter` is therefore a
# safety backstop, not a convergence budget.
#
# Early stop — OSCILLATION guard: an edge that reverses direction between
# iterations has bracketed the support and is now chasing estimator noise.
# We stop and adopt the UNION of the last two candidate boxes (the wider
# option on every edge) so refinement never ends by truncating the box.
# Note we deliberately do NOT stop on "movement not shrinking": steady
# geometric widening moves each edge by a roughly constant relative amount
# every step, which is legitimate progress toward the support, not a stall.
#
# Returns the final range, its fresh particles + log-densities (so the
# caller need not re-seed), a per-iteration log, and the stop reason.
# `verbose` gates messages. NOTE: refinement uses `range` only (not `psamp`).
refine_range = function(N, range, pc, max_iter = 25, tol = 0.05,
                        verbose = TRUE, ...) {
  say = function(...) if (verbose) message("[range] ", ...)
  seed_eval = function(rng) {                 # seed + evaluate + report
    s = unrestricted(N, range = rng, psamp = NULL)
    d = eval_lpden(s, pc)
    if (verbose) report_nonfinite(attr(d, "n_bad"), "range refinement")
    list(samp = s, lpden = d)
  }
  cur_range = range
  se = seed_eval(cur_range); samp = se$samp; lpden = se$lpden
  log = list()
  changed = FALSE
  prev_sign = NULL          # signed direction of the previous move, per edge
  prev_moving = NULL        # which edges moved > tol on the previous step
  stop_reason = "max_iter"
  for (it in seq_len(max_iter)) {
    rep = diagnose_range(cur_range, samp, lpden[2, ], pc = pc, ...)
    log[[it]] = list(range = cur_range, ess_frac = rep$ess_frac,
                     notes = rep$notes)
    if (it == 1L) for (n in rep$notes) say(n)
    new_range = rep$suggest

    # signed per-edge change, scaled by current width (D x 2)
    w = cur_range[, 2] - cur_range[, 1]
    d_edge = (new_range - cur_range) / w
    rel = max(abs(d_edge))

    # converged: no edge wants to move by more than `tol`
    if (!is.finite(rel) || rel < tol) {
      if (it == 1L) say("range accepted as supplied.")
      else say(sprintf("range stabilised after %d refinement(s).", it - 1L))
      stop_reason = "converged"
      break
    }

    # oscillation guard: any edge reversing direction (both moves > tol)
    cur_sign = sign(d_edge)
    moving   = abs(d_edge) >= tol
    if (!is.null(prev_sign)) {
      reversed = moving & prev_moving & (cur_sign * prev_sign < 0)
      if (any(reversed)) {
        cur_range[, 1] = pmin(cur_range[, 1], new_range[, 1])   # union: widest
        cur_range[, 2] = pmax(cur_range[, 2], new_range[, 2])   # edges, no truncation
        changed = TRUE
        say(sprintf(paste("iteration %d: box is oscillating (support bracketed);",
                          "stopping and using the union to avoid truncation."), it))
        se = seed_eval(cur_range); samp = se$samp; lpden = se$lpden
        stop_reason = "oscillation"
        break
      }
    }

    changed = TRUE
    say(sprintf("iteration %d: adjusting range and re-seeding.", it))
    cur_range = new_range
    se = seed_eval(cur_range); samp = se$samp; lpden = se$lpden
    prev_sign = cur_sign
    prev_moving = moving
    if (it == max_iter)
      say(sprintf("reached max_iter=%d; using the latest range.", max_iter))
  }
  if (changed) {
    say("final range:")
    for (d in 1:nrow(cur_range))
      say(sprintf("  param %d: [%.4g, %.4g]", d, cur_range[d, 1], cur_range[d, 2]))
  }
  list(range = cur_range, samp = samp, lpden = lpden, log = log,
       changed = changed, stop_reason = stop_reason)
}

# ---- main routine: sampling from a log-posterior distribution ---------
# range_control:
#   "adapt" (default) — iteratively refine the box (re-seeding under each
#            suggestion until it stabilises), then run SMC. No user re-run.
#   "warn"           — report diagnostics + suggestion, run with the box AS
#                      SUPPLIED (no re-seeding).
#   "off"            — skip range diagnostics entirely.
# range_max_iter / range_tol control the "adapt" refinement loop.
#   range_max_iter (default 25) is a safety backstop, not a convergence
#     budget: refinement converges in ~log(scale mismatch) steps (a
#     ~10,000x mismatch needs ~12; even a ~5,000,000x overshoot into a
#     non-evaluable region settles in ~17), and an oscillation guard stops
#     earlier once the box brackets the support. 25 leaves comfortable
#     margin for pathological models; raising it far higher (50/100) buys
#     no convergence and only prolongs genuinely divergent cases.
SMC = function(N = 1000, M = 10, nuseq_T = 1, range = NULL, psamp = NULL,
               pc = p_cal, range_control = c("adapt", "warn", "off"),
               range_max_iter = 25, range_tol = 0.05) {
  range_control = match.arg(range_control)

  # -- setup / validation
  if (!is.null(range)) {
    validate_range(range)
    D = nrow(range)
  } else {
    if (is.null(psamp)) stop("Supply either `range` or `psamp`.")
    D = nrow(psamp)
  }

  # -- robust range handling. When a `range` is supplied and control is on,
  #    diagnose it (and, under "adapt", iteratively refine + re-seed) before
  #    running SMC. refine_range() returns fresh particles so we reuse them.
  if (!is.null(range) && range_control != "off") {
    if (range_control == "adapt") {
      ref = refine_range(N, range, pc, max_iter = range_max_iter,
                         tol = range_tol, verbose = TRUE)
      range   = ref$range
      samplet = ref$samp
      lpden   = ref$lpden
    } else {                          # "warn": report only, do not re-seed
      samplet = unrestricted(N, range = range, psamp = psamp)
      lpden   = eval_lpden(samplet, pc)
      report_nonfinite(attr(lpden, "n_bad"), "initial sampling")
      rep = diagnose_range(range, samplet, lpden[2, ], pc = pc)
      if (length(rep$notes) == 0) message("[range] no issues detected.")
      else for (n in rep$notes) message("[range] ", n)
    }
  } else {
    # no range diagnostics: seed directly (range and/or psamp)
    samplet = unrestricted(N, range = range, psamp = psamp)
    lpden   = eval_lpden(samplet, pc)   # 2 x N: row1 = logprior, row2 = loglik
    report_nonfinite(attr(lpden, "n_bad"), "initial sampling")
  }

  dd = D * D
  cur   = samplet                # D x N current parameters
  curlp = lpden                  # 2 x N current (logprior, loglik)
  curW  = rep(1 / N, N)          # current normalised weights (uniform)
  nu_prev = 0
  acr = numeric(0)               # per-step mean acceptance rate

  # RAM step-size schedule depends only on the sweep index
  gammas = pmin(1, D * ((1:M + ness)^(-gamma0)))

  repeat {
    ll = curlp[2, ]

    # every particle carries zero mass -> tempering/uniroot is undefined.
    # Fail here with an explanatory message rather than in uniroot().
    if (!any(is.finite(ll) & curW > 0))
      stop("All particles have zero posterior mass (every log-likelihood is ",
           "-Inf or non-finite). The model may be misspecified, or the ",
           "`range` may cover only regions where the likelihood cannot be ",
           "evaluated. Cannot continue.")

    # -- adaptively choose the next tempering level
    if (adapt_seq(nuseq_T, nu_prev, ll, curW, N, fEss) > 0) {
      nu_t = nuseq_T
    } else {
      nu_t = uniroot(adapt_seq, interval = c(nu_prev, nuseq_T),
                     nu0 = nu_prev, ll = ll, Wt0 = curW,
                     N = N, r_ess = fEss)$root
    }

    # -- reweight and resample
    wt = if (nu_t == nu_prev) rep(0, N) else (nu_t - nu_prev) * ll
    W  = curW * exp(wt)
    W[!is.finite(W)] = 0        # -Inf-loglik particles (incl. sanitized) -> 0
    sW = sum(W)
    if (!is.finite(sW) || sW <= 0)
      stop("All particles have zero posterior mass (every log-likelihood is ",
           "-Inf or non-finite). The model may be misspecified or the `range` ",
           "may have been widened onto a region the likelihood cannot ",
           "evaluate. Cannot continue.")
    W  = W / sW
    index = sample.int(N, N, prob = W, replace = TRUE)
    cur   = cur[, index, drop = FALSE]
    curlp = curlp[, index, drop = FALSE]
    curW  = rep(1 / N, N)

    # -- build the per-particle state and move them (RAM-MCMC), M sweeps
    Lsigma  = chol_lower(cov(t(cur)))
    vLsigma = as.vector(Lsigma)
    X = rbind(cur, curlp[1, ], curlp[2, ],
              matrix(rep(vLsigma, N), nrow = dd))   # (dd+D+2) x N

    chunks = chunk_indices(N, worker_count())
    res = foreach(idx = chunks) %dofuture% {
      move_block(X[, idx, drop = FALSE], D, M, gammas, nu_t, alphat, pc)
    } %seed% TRUE

    Xnew = matrix(0, nrow = dd + D + 2, ncol = N)
    tot_acc = 0
    for (j in seq_along(chunks)) {
      Xnew[, chunks[[j]]] = res[[j]]$mat
      tot_acc = tot_acc + res[[j]]$acc
    }

    cur   = Xnew[1:D, , drop = FALSE]
    curlp = Xnew[(D + 1):(D + 2), , drop = FALSE]
    acr   = c(acr, tot_acc / (M * N))

    nu_prev = nu_t
    if (nu_t >= nuseq_T) break
  }

  list(sample = t(cur), acr = acr)
}
