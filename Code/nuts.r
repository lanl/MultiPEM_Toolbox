########################################################################
# nuts.r
#
# A No-U-Turn Sampler (NUTS) for R.
#
# The trajectory machinery (leapfrog, recursive tree doubling, U-turn
# checks, proposal selection) runs in compiled C++ (nuts_core.cpp, via
# Rcpp/RcppEigen). It implements the modern algorithm used by Stan:
#
#   * multinomial sampling with biased progressive selection (Betancourt 2017);
#   * gradient-threaded leapfrog (each gradient reused across steps);
#   * the initial joint density is computed once per iteration, not per leaf;
#   * dual-averaging step-size adaptation (Hoffman & Gelman 2014);
#   * windowed diagonal mass-matrix (metric) adaptation via Welford variance.
#
# If the user does not supply grad_f, numerical gradients are used (numDeriv
# Richardson by default, or a central difference).
#
# References:
#   Hoffman & Gelman (2014) JMLR 15:1593-1623.
#   Betancourt (2017) arXiv:1701.02434.
#
# Usage:
#   source("helpers.r"); source("nuts.r")
#   fit <- NUTS(theta0, f, grad_f, n_iter = 2000)          # trace matrix
#   fit <- NUTS(theta0, f, grad_f, n_iter = 2000, diagnostics = TRUE)
#
# © 2023. Triad National Security, LLC. All rights reserved.           
# This program was produced under U.S. Government contract             
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which   
# is operated by Triad National Security, LLC for the U.S. Department  
# of Energy/National Nuclear Security Administration. All rights in    
# the program are reserved by Triad National Security, LLC, and the    
# U.S. Department of Energy/National Nuclear Security Administration.  
# The Government is granted for itself and others acting on its behalf 
# a nonexclusive, paid-up, irrevocable worldwide license in this       
# material to reproduce, prepare derivative works, distribute copies   
# to the public, perform publicly and display publicly, and to permit  
# others to do so.                                                     
########################################################################

if (!exists("compile_nuts_core")) {
  # Auto-source the helper if the user only sourced this file. Look next to
  # this file first (so it works from any working directory), then fall back
  # to the working directory.
  .this_dir <- tryCatch(
    normalizePath(dirname(sys.frame(1)$ofile), mustWork = FALSE),
    error = function(e) getwd())
  cand <- c(file.path(.this_dir, "helpers.r"), file.path(getwd(), "helpers.r"))
  hp <- cand[file.exists(cand)][1]
  if (!is.na(hp)) source(hp) else
    stop("nuts.r: please source('helpers.r') first.")
}

#' No-U-Turn Sampler.
#'
#' @param theta   Initial parameter vector.
#' @param f       Log-density (up to a constant), function of a numeric vector.
#' @param grad_f  Gradient of f. If NULL, computed numerically (see numeric_grad).
#' @param n_iter  Total number of iterations (including warmup).
#' @param warmup  Number of warmup/adaptation iterations. Default floor(n_iter/2).
#'                Warmup draws are discarded from the returned trace.
#' @param M_diag  Optional fixed diagonal mass matrix. If NULL and adapt_mass is
#'                TRUE it is adapted; if adapt_mass is FALSE it defaults to ones.
#' @param adapt_mass Adapt a diagonal mass matrix during warmup (default TRUE).
#' @param delta   Target mean Metropolis accept probability (default 0.8).
#' @param max_treedepth Maximum NUTS tree depth (default 10).
#' @param eps     Initial step-size guess (a reasonable value is found from it).
#' @param numeric_grad One of "none" (require grad_f), "richardson", "simple".
#'                Ignored when grad_f is supplied.
#' @param backend "callback" (R target functions) or "xptr" (compiled target;
#'                f_ptr/g_ptr must then be supplied instead of f/grad_f).
#' @param f_ptr,g_ptr XPtr externalptr's for the compiled target (backend="xptr").
#' @param seed    Optional RNG seed.
#' @param diagnostics If TRUE, return a list with the trace plus per-iteration
#'                diagnostics; otherwise return just the post-warmup trace matrix.
#' @param verbose Print adaptation progress.
#' @return A matrix (post-warmup draws x parameters), or a list if diagnostics.
NUTS <- function(theta, f = NULL, grad_f = NULL, n_iter,
                 warmup = floor(n_iter / 2),
                 M_diag = NULL, adapt_mass = TRUE,
                 delta = 0.8, max_treedepth = 10, eps = 1,
                 numeric_grad = c("richardson", "simple", "none"),
                 backend = c("callback", "xptr", "fused"),
                 f_ptr = NULL, g_ptr = NULL, fg = NULL,
                 seed = NULL, diagnostics = FALSE, verbose = TRUE) {

  numeric_grad <- match.arg(numeric_grad)
  # Auto-select the fused backend when a fused evaluator is supplied and the
  # caller did not explicitly pick a backend.
  if (missing(backend) && !is.null(fg)) backend <- "fused"
  backend <- match.arg(backend)
  if (!is.null(seed)) set.seed(seed)
  compile_nuts_core()

  d <- length(theta)
  warmup <- min(warmup, n_iter)

  # ---- Resolve the gradient (analytic or numerical) for the callback backend.
  if (backend == "callback") {
    if (is.null(f)) stop("NUTS: f is required for backend='callback'.")
    if (is.null(grad_f)) {
      if (numeric_grad == "none")
        stop("NUTS: grad_f is NULL and numeric_grad='none'.")
      grad_f <- make_numeric_grad(f, method = numeric_grad)
    }
  } else if (backend == "fused") {
    if (is.null(fg)) stop("NUTS: backend='fused' requires fg (a function ",
                          "returning list(logp=, grad=)).")
    # Derive single-purpose f / grad_f from fg for the rare paths that need
    # only one (find_reasonable_epsilon). These call fg and select one field,
    # so they return exactly the same values the hot path uses -- keeping the
    # RNG stream and hence the draws identical to a matched callback run.
    if (is.null(f))      f      <- function(x) fg(x)$logp
    if (is.null(grad_f)) grad_f <- function(x) fg(x)$grad
  } else {
    if (is.null(f_ptr) || is.null(g_ptr))
      stop("NUTS: backend='xptr' requires f_ptr and g_ptr.")
  }

  # ---- Initial metric (inv_M == diagonal of M^{-1} == momentum variances... )
  # Convention: momentum r ~ N(0, M); inv_M = 1/M_diag. We store inv_M.
  if (!is.null(M_diag)) {
    inv_M <- 1 / M_diag
    adapt_mass <- FALSE
  } else {
    inv_M <- rep(1, d)
  }

  # A small closure to run one C++ transition regardless of backend.
  # Returns list(theta, accept_stat, tree_depth, n_leapfrog, divergent, energy).
  one_step <- function(theta, eps, inv_M) {
    if (backend == "callback")
      .nuts_env$nuts_step_callback(theta, f, grad_f, eps, inv_M, max_treedepth)
    else if (backend == "fused")
      .nuts_env$nuts_step_fused_callback(theta, fg, eps, inv_M, max_treedepth)
    else
      # nuts_sample_xptr with n_iter=1 acts as a single step.
      local({
        out <- .nuts_env$nuts_sample_xptr(theta, f_ptr, g_ptr, 1L, eps,
                                          matrix(inv_M, nrow = 1), max_treedepth)
        list(theta = out$theta[1, ], accept_stat = out$accept_stat[1],
             tree_depth = out$tree_depth[1], n_leapfrog = out$n_leapfrog[1],
             divergent = out$divergent[1] > 0, energy = out$energy[1])
      })
  }

  # ---- Find a reasonable starting step size.
  if (backend == "callback" || backend == "fused") {
    eps <- find_reasonable_epsilon(theta, f, grad_f, inv_M, eps)
  } else {
    # For xptr we just run the dual-averaging from the supplied eps; a coarse
    # search is unnecessary because adaptation converges quickly.
  }
  if (verbose) message(sprintf("Initial step size eps = %.4g", eps))

  da <- dual_avg_init(eps, delta = delta)
  sched <- make_warmup_schedule(warmup)
  wf <- welford_init(d)

  # ---- Storage.
  trace <- matrix(0, n_iter, d)
  accept_stat <- numeric(n_iter)
  tree_depth <- integer(n_iter)
  n_leapfrog <- numeric(n_iter)
  divergent <- logical(n_iter)
  energy <- numeric(n_iter)
  eps_used <- numeric(n_iter)

  cur <- theta
  for (it in seq_len(n_iter)) {
    in_warmup <- it <= warmup
    eps_it <- if (in_warmup) da$eps else da$eps_bar

    step <- one_step(cur, eps_it, inv_M)
    cur <- as.numeric(step$theta)

    trace[it, ] <- cur
    accept_stat[it] <- step$accept_stat
    tree_depth[it] <- step$tree_depth
    n_leapfrog[it] <- step$n_leapfrog
    divergent[it] <- isTRUE(step$divergent)
    energy[it] <- step$energy
    eps_used[it] <- eps_it

    if (in_warmup) {
      # Step-size dual averaging every warmup iteration.
      da <- dual_avg_update(da, step$accept_stat)

      # Metric adaptation inside slow windows.
      if (adapt_mass && sched$is_slow[it]) {
        wf <- welford_update(wf, cur)
        if (sched$window_end[it]) {
          inv_M <- welford_variance(wf)
          wf <- welford_init(d)               # reset accumulator for next window
          # Re-find a reasonable step size under the new metric and restart DA.
          if (backend == "callback" || backend == "fused")
            eps <- find_reasonable_epsilon(cur, f, grad_f, inv_M, da$eps)
          else
            eps <- da$eps
          da <- dual_avg_init(eps, delta = delta)
          if (verbose)
            message(sprintf("  [warmup %d] metric updated; eps reset to %.4g",
                            it, eps))
        }
      }
      if (it == warmup) {
        da$eps <- da$eps_bar   # freeze at the dual-averaging optimum
        if (verbose)
          message(sprintf("Warmup complete: eps = %.4g, mean accept = %.3f",
                          da$eps_bar, mean(accept_stat[seq_len(warmup)])))
      }
    }
  }

  keep <- if (warmup < n_iter) (warmup + 1):n_iter else seq_len(n_iter)
  post <- trace[keep, , drop = FALSE]

  if (!diagnostics) return(post)

  list(
    theta = post,
    warmup_draws = trace[seq_len(warmup), , drop = FALSE],
    all_draws = trace,
    diagnostics = data.frame(
      iter = seq_len(n_iter),
      warmup = seq_len(n_iter) <= warmup,
      accept_stat = accept_stat,
      tree_depth = tree_depth,
      n_leapfrog = n_leapfrog,
      divergent = divergent,
      energy = energy,
      eps = eps_used
    ),
    eps = da$eps_bar,
    inv_M = inv_M,
    M_diag = 1 / inv_M,
    n_divergent = sum(divergent[keep]),
    mean_accept = mean(accept_stat[keep])
  )
}

# ---------------------------------------------------------------------------
# Convenience: run a whole *sampling* chain (fixed eps, fixed metric) entirely
# in C++. Use this when you already have an adapted eps/inv_M and want to draw
# the sampling phase without returning to R between iterations.
# ---------------------------------------------------------------------------
NUTS_sample_fixed <- function(theta, f = NULL, grad_f = NULL, n_iter, eps, inv_M,
                              max_treedepth = 10,
                              backend = c("callback", "xptr"),
                              f_ptr = NULL, g_ptr = NULL, seed = NULL) {
  backend <- match.arg(backend)
  if (!is.null(seed)) set.seed(seed)
  compile_nuts_core()
  eps_seq <- as.numeric(eps)
  imm <- matrix(inv_M, nrow = 1)
  if (backend == "callback")
    .nuts_env$nuts_sample_callback(theta, f, grad_f, n_iter, eps_seq, imm, max_treedepth)
  else
    .nuts_env$nuts_sample_xptr(theta, f_ptr, g_ptr, n_iter, eps_seq, imm, max_treedepth)
}
