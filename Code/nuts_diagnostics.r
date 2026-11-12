########################################################################
# nuts_diagnostics.r
#
# Read-only post-hoc diagnostics for a NUTS fit (the list returned by
# NUTS(..., diagnostics = TRUE)). Nothing here changes the sampler, the draws,
# or the RNG stream -- these functions only summarize data the sampler already
# produced.
#
# divergence_attribution(): locate WHERE in parameter space divergent
# transitions concentrate. Divergences are the fingerprint of geometry the
# leapfrog integrator cannot follow (classically, the neck of a funnel). If they
# cluster in a particular coordinate or block (e.g. the FGSN errors-in-variables
# yields or the variance-component scales), that block is the one to
# reparameterize (non-center) -- a constant mass matrix, diagonal or dense,
# cannot fully fix a funnel. If divergences are ~0 or scattered with no
# concentration, no reparameterization is indicated.
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

#' Attribute divergent transitions to parameters (and optional blocks).
#'
#' Model-agnostic: for each coordinate we robustly standardize the draws using
#' the NON-divergent post-warmup sample (median / MAD), then measure how extreme
#' that coordinate is, on average, at the DIVERGENT iterations. A large score
#' means divergences systematically occur where that coordinate takes extreme
#' values -- the signature of a funnel in that coordinate.
#'
#' @param fit    A list from NUTS(..., diagnostics = TRUE). Must contain
#'               all_draws (n_iter x d) and diagnostics$divergent / $warmup.
#' @param blocks Optional named list mapping block name -> integer coordinate
#'               indices (e.g. list(FGSN = 101:104, model = 1:100)). When given,
#'               per-coordinate scores are aggregated to a per-block summary.
#' @param top    Number of top coordinates to report (default 10).
#' @return A list with: n_divergent, n_post (post-warmup count), per_coord (a
#'         data.frame ranked by score), and per_block (data.frame or NULL).
#'         Returns gracefully (empty, with a message) when there are no
#'         post-warmup divergences.
divergence_attribution <- function(fit, blocks = NULL, top = 10) {
  if (is.null(fit$all_draws) || is.null(fit$diagnostics))
    stop("divergence_attribution: fit must come from NUTS(diagnostics = TRUE) ",
         "(needs all_draws and diagnostics).")
  dg <- fit$diagnostics
  draws <- fit$all_draws
  d <- ncol(draws)

  post <- !dg$warmup
  div  <- post & dg$divergent
  ok   <- post & !dg$divergent
  n_div  <- sum(div)
  n_post <- sum(post)

  empty <- list(n_divergent = n_div, n_post = n_post,
                per_coord = NULL, per_block = NULL)
  if (n_div == 0L) {
    message(sprintf("divergence_attribution: 0/%d post-warmup divergences; ",
                    n_post), "no attribution needed.")
    return(empty)
  }
  if (sum(ok) < 2L) {
    message("divergence_attribution: too few non-divergent draws to form a ",
            "reference scale; reporting raw divergent spread only.")
  }

  # Robust per-coordinate center/scale from the non-divergent reference sample
  # (falls back to all post-warmup draws if the non-divergent set is tiny).
  ref <- if (sum(ok) >= 2L) draws[ok, , drop = FALSE] else draws[post, , drop = FALSE]
  divd <- draws[div, , drop = FALSE]
  ok_med  <- apply(ref,  2, median)
  ok_scale <- apply(ref, 2, mad)
  # Guard constant / near-constant coordinates (mad == 0): a tiny positive scale
  # so they score ~0 rather than Inf.
  ok_scale[!(ok_scale > 0)] <- .Machine$double.eps^0.5

  # Causal funnel signal = robust LOCATION SHIFT: how far the divergent draws'
  # center sits from the non-divergent center, in reference-scale units. A
  # funnel's neck coordinate shifts strongly toward one tail at divergences; a
  # coordinate that merely BLOWS UP reactively (symmetric overshoot) keeps its
  # center and so scores low -- which is what we want (we seek the cause, not
  # the symptom). We also report the spread ratio as secondary context.
  div_med   <- apply(divd, 2, median)
  shift     <- (div_med - ok_med) / ok_scale           # signed, in MADs
  score     <- abs(shift)
  direction <- sign(shift)                             # +1 high tail / -1 low tail
  div_scale <- if (nrow(divd) >= 2L) apply(divd, 2, mad) else rep(NA_real_, d)
  spread_ratio <- div_scale / ok_scale                 # >1 => also blows up here

  per_coord <- data.frame(
    coord        = seq_len(d),
    score        = score,
    direction    = direction,
    spread_ratio = spread_ratio,
    stringsAsFactors = FALSE
  )
  per_coord <- per_coord[order(-per_coord$score), , drop = FALSE]
  rownames(per_coord) <- NULL
  per_coord_top <- head(per_coord, top)

  per_block <- NULL
  if (!is.null(blocks)) {
    per_block <- do.call(rbind, lapply(names(blocks), function(nm) {
      idx <- blocks[[nm]]
      idx <- idx[idx >= 1 & idx <= d]
      if (!length(idx)) return(NULL)
      data.frame(block = nm, n_coord = length(idx),
                 mean_score = mean(score[idx]),
                 max_score  = max(score[idx]),
                 stringsAsFactors = FALSE)
    }))
    if (!is.null(per_block)) {
      per_block <- per_block[order(-per_block$max_score), , drop = FALSE]
      rownames(per_block) <- NULL
    }
  }

  list(n_divergent = n_div, n_post = n_post,
       per_coord = per_coord_top, per_block = per_block)
}

#' Pretty-print a divergence_attribution() result to the transcript.
print_divergence_attribution <- function(att, label = "") {
  hdr <- if (nzchar(label)) paste0("DIVERGENCE ATTRIBUTION (", label, ")")
         else "DIVERGENCE ATTRIBUTION"
  cat("\n", hdr, "\n", sep = "")
  cat(sprintf("  %d / %d post-warmup divergences\n", att$n_divergent, att$n_post))
  if (att$n_divergent == 0L) { cat("  (no attribution needed)\n"); return(invisible()) }
  if (!is.null(att$per_block)) {
    cat("  by block (ranked by max location-shift; higher = divergences ",
        "concentrate in that block's tail):\n", sep = "")
    for (i in seq_len(nrow(att$per_block))) {
      b <- att$per_block[i, ]
      cat(sprintf("    %-14s n=%-3d  mean shift=%.2f  max shift=%.2f (MADs)\n",
                  b$block, b$n_coord, b$mean_score, b$max_score))
    }
  }
  cat("  top coordinates (robust location shift at divergences, in MADs; ",
      "spread = divergent/reference MAD):\n", sep = "")
  for (i in seq_len(nrow(att$per_coord))) {
    r <- att$per_coord[i, ]
    tail <- if (r$direction > 0) "high" else if (r$direction < 0) "low" else "none"
    sr <- if (is.finite(r$spread_ratio)) sprintf("%.1fx", r$spread_ratio) else "n/a"
    cat(sprintf("    coord %-4d  shift=%.2f (%s tail)  spread=%s\n",
                r$coord, r$score, tail, sr))
  }
  invisible()
}

#' Build a coordinate->block map from a MultiPEM p_cal. Robustly identifies the
#' FGSN errors-in-variables block (the last p_fgsn coordinates) and the scale
#' block (the p_A coordinates preceding it); everything else is "model". These
#' are the blocks most prone to funnel geometry. Returns a named list suitable
#' for the `blocks` argument of divergence_attribution().
nuts_blocks_from_pcal <- function(p_cal) {
  np <- p_cal$npars; nm <- p_cal$nmpars
  pf <- if (is.null(p_cal$p_fgsn)) 0L else p_cal$p_fgsn
  pA <- if (is.null(p_cal$p_A))    0L else p_cal$p_A
  # Parameter-vector layout is [ model (nmpars) | FGSN (p_fgsn) | A (p_A) ]:
  # log_prior.r / glog_prior.r strip A from the tail FIRST, then FGSN from the
  # new tail, so A is the final block and FGSN precedes it. The block index
  # ranges below must follow that order exactly.
  blocks <- list()
  if (nm > 0)  blocks[["model"]]   <- seq_len(nm)
  if (pf > 0)  blocks[["FGSN"]]    <- (nm + 1):(nm + pf)
  if (pA > 0)  blocks[["scale_A"]] <- (nm + pf + 1):(nm + pf + pA)
  blocks
}
