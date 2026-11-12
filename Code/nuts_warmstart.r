########################################################################
# nuts_warmstart.r
#
# Quality-safe NUTS warmup speedup: resolve a DENSE warm-start metric (inv_M)
# that seeds the full run's warmup so it does not begin from the identity metric
# (which, on a correlated posterior, pins at the maximum tree depth for the
# first ~100 iterations -- the dominant cost). The seed only sets the STARTING
# metric; the full run re-adapts normally (NUTS M_init, adapt_mass stays TRUE),
# so it is quality-neutral -- it changes warmup speed, not the final adapted
# metric or step size. This mirrors Stan's metric_file / inv_metric warm-start.
#
# nuts_warmstart (string) selects the source:
#   "auto"        - DEFAULT: use the saved metric file if present in the run dir
#                   (near-zero cost), else run the pilot.
#   "pilot"       - always run a dense NUTS pilot to estimate the seed covariance.
#   "none"        - cold start from the identity metric (legacy behavior).
#   "<file path>" - load a specific saved adapted metric file.
# A missing / wrong-dimension / unreadable file WARNS and falls back to "pilot"
# (never aborts a multi-day run).
#
# Each stack uses its own saved-metric file name (calibration: metric.RData;
# event analysis / runMPEM_0.r: metric_0.RData) to avoid cross-stack confusion.
# The pilot size is fixed in code (not a user/deck knob); it favors seed QUALITY
# over speed (full tree depth, ample warmup), since a good seed is what makes the
# full production run's warmup cheap, and the full run re-adapts regardless.
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

# ---- fixed pilot configuration (intentionally NOT exposed to the deck) -----
# The pilot produces the dense covariance SEED for the full run. SEED QUALITY IS
# THE PRIORITY: the pilot runs only when no saved metric exists, and a good seed
# is what makes the full run's warmup cheap, so we do NOT starve the pilot for
# speed. It uses the full production tree depth (not capped) and an ample warmup
# so the adapted metric is representative. It is still substantially cheaper than
# the full run -- fewer warmup iterations and no production/sampling phase -- and
# the full run re-adapts from the seed regardless (quality-neutral).
.nuts_pilot_nburn <- 300L    # warmup iterations (ample adaptation windows)
.nuts_pilot_nsamp <- 0L      # no sampling phase needed (metric read from warmup)
.nuts_pilot_maxtd <- 10L     # full production tree depth (seed not depth-limited)

#' Run the throwaway dense pilot and return its adapted dense inv_M (d x d), or
#' NULL on failure. Reuses the caller's fused evaluator (lpg) when available,
#' else the separate f / grad_f callbacks -- i.e. the exact target the full run
#' will sample, started from the MAP (x0).
nuts_pilot_metric <- function(x0, lpg = NULL, f = NULL, grad_f = NULL,
                              use_fused = !is.null(lpg), seed = 1L) {
  tryCatch({
    n_iter <- .nuts_pilot_nburn + .nuts_pilot_nsamp
    if (use_fused) {
      fit <- NUTS(x0, fg = lpg, n_iter = n_iter, warmup = .nuts_pilot_nburn,
                  metric = "dense", max_treedepth = .nuts_pilot_maxtd,
                  seed = seed, diagnostics = TRUE, verbose = FALSE)
    } else {
      fit <- NUTS(x0, f = f, grad_f = grad_f, n_iter = n_iter,
                  warmup = .nuts_pilot_nburn, metric = "dense",
                  max_treedepth = .nuts_pilot_maxtd, seed = seed,
                  diagnostics = TRUE, verbose = FALSE)
    }
    M <- as.matrix(fit$inv_M)
    if (!all(is.finite(M)) || nrow(M) != length(x0)) return(NULL)
    (M + t(M)) / 2
  }, error = function(e) {
    warning("nuts_pilot_metric: pilot failed (", conditionMessage(e),
            "); cold start will be used.", call. = FALSE)
    NULL
  })
}

#' Load a saved dense metric from an .RData written by save_nuts_metric().
#' Returns the d x d inv_M matrix if present and correctly dimensioned, else
#' NULL. Prints provenance when available.
load_nuts_metric <- function(path, d) {
  tryCatch({
    if (!file.exists(path)) {
      warning("nuts_warmstart: metric file not found: ", path,
              " -- falling back to pilot.", call. = FALSE)
      return(NULL)
    }
    e <- new.env()
    load(path, envir = e)
    obj <- if (exists("nuts_metric_save", envir = e)) get("nuts_metric_save", envir = e)
           else NULL
    if (is.null(obj) || is.null(obj$inv_M)) {
      warning("nuts_warmstart: no saved metric object in ", path,
              " -- falling back to pilot.", call. = FALSE)
      return(NULL)
    }
    M <- as.matrix(obj$inv_M)
    if (nrow(M) != d || ncol(M) != d) {
      warning(sprintf(paste0("nuts_warmstart: saved metric is %dx%d but model ",
                             "dimension is %d -- falling back to pilot."),
                      nrow(M), ncol(M), d), call. = FALSE)
      return(NULL)
    }
    prov <- if (!is.null(obj$provenance)) obj$provenance else "unknown"
    print(paste0("NUTS warm-start: seeding from saved metric (", prov,
                 "), dim ", d))
    (M + t(M)) / 2
  }, error = function(e) {
    warning("nuts_warmstart: failed to load ", path, " (",
            conditionMessage(e), ") -- falling back to pilot.", call. = FALSE)
    NULL
  })
}

#' Save an adapted dense metric for reuse as a future warm-start seed.
#' @param inv_M d x d adapted covariance-scale metric.
#' @param path destination .RData (default "./metric.RData").
#' @param imcmc sampler label; provenance a free-text tag (e.g. case dir + date).
save_nuts_metric <- function(inv_M, path = "./metric.RData",
                             imcmc = "NUTS", provenance = "") {
  tryCatch({
    M <- as.matrix(inv_M)
    nuts_metric_save <- list(inv_M = (M + t(M)) / 2, dim = nrow(M),
                             imcmc = imcmc, provenance = provenance)
    save(nuts_metric_save, file = path)
    print(paste0("NUTS adapted metric saved to ", path,
                 " (dim ", nrow(M), ") for future warm-start."))
    invisible(TRUE)
  }, error = function(e) {
    warning("save_nuts_metric: could not write ", path, " (",
            conditionMessage(e), ")", call. = FALSE)
    invisible(FALSE)
  })
}

#' Resolve nuts_warmstart (string) to a dense warm-start metric (d x d) or NULL
#' (NULL => cold start).
#'   "auto"  -> if `metric_file` exists in the run dir and is dim-correct, load
#'              it (near-zero cost); otherwise run the pilot. THE DEFAULT.
#'   "none"  -> NULL (cold start).
#'   "pilot" -> always run the pilot (ignore any saved file).
#'   <path>  -> load that specific file; fall back to pilot on any failure.
#' @param metric_file name of the saved-metric file this stack uses (e.g.
#'        "metric.RData" for calibration, "metric_0.RData" for event analysis).
#'        Used only by the "auto" keyword.
resolve_warmstart_metric <- function(nuts_warmstart, x0, d,
                                      lpg = NULL, f = NULL, grad_f = NULL,
                                      use_fused = !is.null(lpg), seed = 1L,
                                      metric_file = "metric.RData") {
  if (is.null(nuts_warmstart) || identical(nuts_warmstart, "none"))
    return(NULL)
  if (identical(nuts_warmstart, "auto")) {
    if (file.exists(metric_file)) {
      print(paste0("NUTS warm-start (auto): found ", metric_file,
                   " -- loading saved metric."))
      M <- load_nuts_metric(metric_file, d)
      if (!is.null(M)) return(M)
      print("NUTS warm-start (auto): saved metric unusable; running pilot.")
    } else {
      print(paste0("NUTS warm-start (auto): no ", metric_file,
                   " in run dir -- running pilot."))
    }
    return(nuts_pilot_metric(x0, lpg, f, grad_f, use_fused, seed))
  }
  if (identical(nuts_warmstart, "pilot")) {
    print("NUTS warm-start: estimating dense seed metric from a pilot run...")
    return(nuts_pilot_metric(x0, lpg, f, grad_f, use_fused, seed))
  }
  # treat as an explicit file path; fall back to pilot on any failure
  M <- load_nuts_metric(nuts_warmstart, d)
  if (!is.null(M)) return(M)
  print("NUTS warm-start: falling back to pilot seed.")
  nuts_pilot_metric(x0, lpg, f, grad_f, use_fused, seed)
}
