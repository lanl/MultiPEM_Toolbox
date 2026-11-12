########################################################################
# mcmc_convergence.r
#
# Hand-rolled MCMC convergence diagnostics for the MultiPEM Bayesian samplers
# (RAM, FME, NUTS). No external package dependency. Implements the
# rank-normalized split-Rhat and bulk/tail effective sample size of
# Vehtari, Gelman, Simpson, Carpenter & Burkner (2021), "Rank-normalization,
# folding, and localization: An improved Rhat for assessing convergence of
# MCMC", Bayesian Analysis 16(2):667-718.
#
# Targets:
#   Rhat   <= 1.01  : converged (per that paper's recommendation)
#   ESS    >= 400   : adequate (>= ~100 per chain rule of thumb, bulk & tail)
#
# All reporting functions are wrapped so a diagnostic failure can never abort a
# production run -- they print a note and return invisibly on error.
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

# ---- low-level helpers ----------------------------------------------------

# Autocovariance of a single series at all lags via FFT (fast, exact).
.mc_acov <- function(x) {
  n <- length(x)
  x <- x - mean(x)
  nfft <- 2^ceiling(log2(2 * n))
  f <- stats::fft(c(x, rep(0, nfft - n)))
  ac <- Re(stats::fft(f * Conj(f), inverse = TRUE)) / nfft
  ac[seq_len(n)] / n
}

# Rank-normalize a matrix of draws (iter x chain): average rank across all
# values, then map to normal scores (the Blom transform used by Vehtari 2021).
.mc_rank_normalize <- function(mat) {
  r <- rank(mat, ties.method = "average")
  z <- stats::qnorm((r - 3/8) / (length(r) + 1/4))
  matrix(z, nrow = nrow(mat), ncol = ncol(mat))
}

# Split-Rhat on a (iter x chain) matrix of a single parameter (already
# rank-normalized if desired). Splits each chain in half -> 2*chains halves.
.mc_split_rhat <- function(mat) {
  n <- nrow(mat); m <- ncol(mat)
  if (n < 4L) return(NA_real_)
  h <- floor(n / 2)
  halves <- cbind(mat[1:h, , drop = FALSE], mat[(n - h + 1):n, , drop = FALSE])
  # halves is h x (2m)
  nn <- nrow(halves); mm <- ncol(halves)
  cm <- colMeans(halves)
  gm <- mean(cm)
  B <- nn / (mm - 1) * sum((cm - gm)^2)            # between
  W <- mean(apply(halves, 2, stats::var))          # within
  if (!is.finite(W) || W <= 0) return(NA_real_)
  varhat <- (nn - 1) / nn * W + B / nn
  sqrt(varhat / W)
}

# Bulk/tail ESS on a (iter x chain) matrix of a single parameter.
# Uses Geyer's initial monotone sequence estimator on the mean autocorrelation
# across chains (the standard multi-chain ESS; Vehtari 2021 eq. 10-13).
.mc_ess <- function(mat) {
  n <- nrow(mat); m <- ncol(mat)
  if (n < 4L) return(NA_real_)
  # per-chain autocovariance, averaged; combine with between-chain variance
  acov <- apply(mat, 2, .mc_acov)                  # n x m
  chain_mean_var <- apply(mat, 2, stats::var) * (n - 1) / n
  mean_var <- mean(chain_mean_var)
  if (!is.finite(mean_var) || mean_var <= 0) return(NA_real_)
  cm <- colMeans(mat)
  B <- if (m > 1) stats::var(cm) else 0
  var_plus <- mean_var * (n - 1) / n + B
  rho <- rep(0, n)
  mean_acov <- rowMeans(acov)                      # average over chains
  rho <- 1 - (mean_var - mean_acov) / var_plus     # combined autocorrelation
  # Geyer initial positive + monotone sequence on paired sums
  t_max <- n - 1
  # pair sums P_k = rho_{2k} + rho_{2k+1}
  np <- floor((t_max - 1) / 2)
  if (np < 1) return(m * n)
  P <- rho[2 * (1:np)] + rho[2 * (1:np) + 1]
  # initial positive sequence
  pos <- which(P < 0)
  k_stop <- if (length(pos)) pos[1] - 1 else np
  if (k_stop < 1) return(m * n)
  P <- P[1:k_stop]
  # initial monotone sequence (make non-increasing). Iterate with seq_along(P)[-1]
  # (empty when length(P)==1) rather than 2:length(P), which would count DOWN to
  # c(2,1) for a length-1 P and read the out-of-bounds P[2]==NA (R's 2:length
  # footgun). A length-1 P is common for a well-mixed parameter whose paired
  # autocorrelation sequence stops at the first lag.
  for (k in seq_along(P)[-1]) if (P[k] > P[k - 1]) P[k] <- P[k - 1]
  tau <- -1 + 2 * sum(P)                            # integrated autocorr time
  if (!is.finite(tau) || tau < 1) tau <- 1
  (m * n) / tau
}

# ---- public API -----------------------------------------------------------

#' Per-parameter rank-normalized split-Rhat and bulk/tail ESS from a list of
#' per-chain draw matrices (each iter x npars, all the same dims).
#' @return data.frame(param, rhat, ess_bulk, ess_tail)
rhat_ess <- function(chain_list) {
  m <- length(chain_list)
  if (m < 1L) stop("rhat_ess: need >= 1 chain")
  d <- ncol(chain_list[[1]])
  n <- nrow(chain_list[[1]])
  rhat <- ess_bulk <- ess_tail <- rep(NA_real_, d)
  for (j in seq_len(d)) {
    # iter x chain matrix for parameter j
    mat <- sapply(chain_list, function(ch) ch[, j])
    if (is.null(dim(mat))) mat <- matrix(mat, ncol = m)
    if (any(!is.finite(mat))) next
    # bulk: rank-normalized draws
    z <- .mc_rank_normalize(mat)
    rhat[j]     <- .mc_split_rhat(z)
    ess_bulk[j] <- .mc_ess(z)
    # tail: rank-normalized folded draws (|x - median|)
    fold <- abs(mat - stats::median(mat))
    zf <- .mc_rank_normalize(fold)
    ess_tail[j] <- .mc_ess(zf)
  }
  data.frame(param = seq_len(d), rhat = rhat,
             ess_bulk = ess_bulk, ess_tail = ess_tail)
}

#' Within-chain split-Rhat/ESS for a SINGLE chain (used for the nimp>1
#' multiple-imputation case, where each imputation is its own target). Splits
#' the one chain into halves and treats them as 2 pseudo-chains.
split_rhat_ess <- function(single_chain) {
  n <- nrow(single_chain); d <- ncol(single_chain)
  h <- floor(n / 2)
  if (h < 2L) return(data.frame(param = seq_len(d), rhat = NA_real_,
                                ess_bulk = NA_real_, ess_tail = NA_real_))
  rhat_ess(list(single_chain[1:h, , drop = FALSE],
                single_chain[(n - h + 1):n, , drop = FALSE]))
}

#' Print a compact convergence report. chain_list is a list of per-chain draw
#' matrices; labels is a length-npars character vector (from make_param_labels).
#' kind tags the header: "cross-chain", "within-chain", or "smc".
#' @param ndraw_note optional string appended to the header (e.g. chain count)
report_convergence <- function(chain_list, labels = NULL,
                               kind = c("cross-chain", "within-chain", "smc"),
                               top = 10, rhat_warn = 1.01, ess_warn = 400) {
  kind <- match.arg(kind)
  tryCatch({
    hdr <- switch(kind,
      "cross-chain"  = "MCMC CONVERGENCE (cross-chain rank-normalized split-R-hat / ESS)",
      "within-chain" = "MCMC CONVERGENCE (within-chain split-R-hat / ESS, per imputation)",
      "smc"          = "MCMC CONVERGENCE")
    cat("\n", hdr, "\n", sep = "")
    if (kind == "smc") {
      cat("  R-hat/ESS not applicable -- SMC is a particle sampler ",
          "(no Markov chains to compare).\n", sep = "")
      return(invisible())
    }
    if (length(chain_list) < 1L || is.null(chain_list[[1]])) {
      cat("  (no chains available to diagnose)\n"); return(invisible())
    }
    d <- ncol(chain_list[[1]])
    if (is.null(labels) || length(labels) != d)
      labels <- paste0("param", seq_len(d))

    diag <- rhat_ess(chain_list)
    nchain <- length(chain_list); niter <- nrow(chain_list[[1]])
    cat(sprintf("  %d chains x %d draws, %d parameters\n", nchain, niter, d))

    rh <- diag$rhat; eb <- diag$ess_bulk; et <- diag$ess_tail
    nbad_rhat <- c(sum(rh > 1.01, na.rm = TRUE),
                   sum(rh > 1.05, na.rm = TRUE),
                   sum(rh > 1.10, na.rm = TRUE))
    cat(sprintf("  max R-hat = %.3f  (# >1.01: %d, >1.05: %d, >1.10: %d of %d)\n",
                max(rh, na.rm = TRUE), nbad_rhat[1], nbad_rhat[2], nbad_rhat[3], d))
    cat(sprintf("  min bulk-ESS = %.0f, min tail-ESS = %.0f  (# bulk <%d: %d)\n",
                min(eb, na.rm = TRUE), min(et, na.rm = TRUE), ess_warn,
                sum(eb < ess_warn, na.rm = TRUE)))

    verdict <- if (max(rh, na.rm = TRUE) <= rhat_warn &&
                   min(eb, na.rm = TRUE) >= ess_warn)
                 "CONVERGED (all R-hat <= 1.01 and bulk-ESS >= 400)"
               else "NOT CONVERGED -- inspect the parameters below"
    cat("  verdict: ", verdict, "\n", sep = "")

    worst <- order(-rh)[seq_len(min(top, d))]
    worst <- worst[is.finite(rh[worst]) & rh[worst] > rhat_warn]
    if (length(worst)) {
      cat("  worst parameters (by R-hat):\n")
      for (i in worst)
        cat(sprintf("    R-hat=%.2f bulk-ESS=%.0f tail-ESS=%.0f  %s\n",
                    rh[i], eb[i], et[i], labels[i]))
    }
    invisible(diag)
  }, error = function(e)
    cat("  (convergence report skipped: ", conditionMessage(e), ")\n", sep = ""))
}
