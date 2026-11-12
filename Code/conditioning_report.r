########################################################################
# conditioning_report.r
#
# Report the numerical conditioning of an information (or covariance) matrix
# produced by MLE / MAP estimation, and -- crucially -- name the PARAMETER
# DIRECTIONS the data poorly constrain. An ill-conditioned information matrix
# means some linear combinations of parameters are nearly undetermined by the
# data; those directions produce near-infinite posterior variances, unreliable
# standard errors ("SE undefined"), and poor MCMC mixing. This report makes that
# explicit instead of leaving it to be discovered via NaN warnings downstream.
#
# For an INFORMATION matrix (which="information"): small eigenvalues are the
# ill-determined directions. For a COVARIANCE matrix (which="covariance"): large
# eigenvalues are the ill-determined directions. Either way we report the
# condition number and, for each flagged direction, the parameters with the
# largest loadings in that eigenvector -- mapped through human-readable labels.
#
# Pure side-effect print; returns its input invisibly; wrapped so it can never
# abort the run.
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

#' @param M matrix (base matrix or Matrix object) -- information or covariance.
#' @param which "information" (default) or "covariance" -- sets which end of the
#'        eigenvalue spectrum is the ill-determined one.
#' @param label short string naming the matrix (e.g. "MAP information (calp)").
#' @param param_names optional length-nrow(M) labels (from make_param_labels).
#' @param tol relative eigenvalue threshold flagging a near-null direction
#'        (default 1e-8 of the largest eigenvalue).
#' @param n_dir max number of ill-determined directions to detail (default 5).
#' @param n_load max parameter loadings to show per direction (default 4).
report_conditioning <- function(M, which = c("information", "covariance"),
                                label = "", param_names = NULL,
                                tol = 1e-8, n_dir = 5, n_load = 4) {
  which <- match.arg(which)
  tryCatch({
    A <- as.matrix(M)
    A <- (A + t(A)) / 2                              # enforce symmetry
    d <- nrow(A)
    if (is.null(param_names) || length(param_names) != d)
      param_names <- paste0("param", seq_len(d))

    ev <- eigen(A, symmetric = TRUE)                 # values descending
    val <- ev$values
    vec <- ev$vectors
    amax <- max(abs(val))
    if (!is.finite(amax) || amax <= 0) {
      cat(sprintf("\nCONDITIONING [%s]: degenerate matrix (max|eigenvalue|=%.3g)\n",
                  label, amax)); return(invisible(M))
    }
    # condition number (ratio of largest to smallest positive eigenvalue)
    pos <- val[val > 0]
    cond <- if (length(pos)) max(pos) / min(pos) else Inf
    nneg <- sum(val < -tol * amax)                   # negative = not PD

    # ill-determined directions: smallest eigenvalues for information,
    # largest for covariance.
    flagged <- which(abs(val) < tol * amax)
    if (which == "covariance") {
      # for covariance the degenerate directions are the LARGEST variances;
      # flag eigenvalues > (1/tol)*min as the analogue
      posv <- val[val > 0]
      if (length(posv)) flagged <- which(val > min(posv) / tol)
    }

    cat(sprintf("\nCONDITIONING [%s] (%s matrix, d=%d)\n", label, which, d))
    cat(sprintf("  condition number = %.3g; eigenvalue range [%.3g, %.3g]\n",
                cond, min(val), max(val)))
    if (nneg > 0)
      cat(sprintf("  WARNING: %d negative eigenvalue(s) -- matrix is NOT ",
                  nneg), "positive-definite (saddle/!converged MLE-MAP).\n", sep = "")
    if (!length(flagged)) {
      cat(sprintf("  no directions below tol (%.0e x max); well-conditioned.\n", tol))
      return(invisible(M))
    }

    # order flagged directions by how ill-determined (most degenerate first)
    ord <- if (which == "information") flagged[order(val[flagged])]
           else flagged[order(-val[flagged])]
    nshow <- min(n_dir, length(ord))
    cat(sprintf(paste0("  %d ill-determined direction(s) (showing %d); for ",
                       "each, the parameters with the largest loadings:\n"),
                length(ord), nshow))
    for (k in seq_len(nshow)) {
      j <- ord[k]
      loading <- vec[, j]
      top <- order(-abs(loading))[seq_len(min(n_load, d))]
      parts <- sprintf("%s (%.2f)", param_names[top], loading[top])
      cat(sprintf("    direction %d (eigenvalue %.3g): %s\n",
                  k, val[j], paste(parts, collapse = ", ")))
    }
    invisible(M)
  }, error = function(e)
    cat("\nCONDITIONING [", label, "]: report skipped (",
        conditionMessage(e), ")\n", sep = ""))
}

#' Convenience wrapper used at the obs_info / info_ll call sites. Resolves
#' human-readable parameter labels from `pc` via make_param_labels() when that
#' helper is available AND its length matches the matrix dimension; otherwise
#' passes NULL (the report then falls back to "paramK"). Fully guarded: if
#' report_conditioning is unavailable or anything fails, it is a silent no-op, so
#' the obs_info/info_ll functions are safe whether or not this file was sourced.
report_conditioning_pc <- function(M, pc, which = "information", label = "") {
  tryCatch({
    if (!exists("report_conditioning")) return(invisible(M))
    d <- nrow(as.matrix(M))
    pnm <- NULL
    if (exists("make_param_labels")) {
      cand <- tryCatch(make_param_labels(pc, target = d), error = function(e) NULL)
      if (!is.null(cand) && length(cand) == d) pnm <- cand
    }
    report_conditioning(M, which = which, label = label, param_names = pnm)
  }, error = function(e) invisible(M))
}
