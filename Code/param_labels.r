########################################################################
# param_labels.r
#
# Build a human-readable label for every column of the MultiPEM parameter
# vector (the columns of the posterior-sample matrix `mpi`, equivalently the
# rows/cols of the MLE/MAP information matrix). The labels name the
# phenomenology, model block, emplacement/response, and within-block index --
# e.g. "Phen 1 | Emplacement 2 | Response 1 | coef 3" -- so the convergence
# (R-hat/ESS) and conditioning reports can point at the exact parameter(s)
# affected instead of a bare coordinate number.
#
# The column ORDER replicated here is exactly the order in which
# print_sumstats.r / print_sumstats_0.r consume `mpi` (they slice leading
# columns off left-to-right). Any change to that order must be mirrored here.
#   calibration layout (make_param_labels):
#     [new-event theta0 if nev] [calp if ncalp>0] [EIV sources] [common coeffs]
#     [emplacement coeffs] [source VC] [path VC] [obs-error cov] [FGSN]
#
# These functions are PURE (no printing, no side effects). They never error:
# if the reconstructed block counts do not match the expected total, the
# remaining columns are labelled "paramK" and the vector is padded/truncated to
# the requested length so callers always get a usable, correctly-sized vector.
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

# Field accessor that works whether pc is an environment or a plain list
# (p_cal is an environment in the live toolbox but may load as a list from an
# .RData). Returns NULL when absent.
.pl_get <- function(pc, name) {
  if (is.environment(pc)) {
    if (exists(name, where = pc, inherits = FALSE)) get(name, envir = pc) else NULL
  } else {
    pc[[name]]
  }
}

# Safe accessor: integer field or 0 when unset/NULL.
.pl_int <- function(pc, name) {
  v <- .pl_get(pc, name)
  if (is.null(v) || !length(v)) 0L else as.integer(v)
}

# Pad/truncate a label vector to exactly `target` entries, filling any shortfall
# with generic "paramK" labels. Guards block-count drift and other applications'
# layouts so the caller always gets a correctly-sized vector.
.pl_fit <- function(labels, target) {
  n <- length(labels)
  if (n == target) return(labels)
  if (n > target)  return(labels[seq_len(target)])
  c(labels, paste0("param", (n + 1):target))
}

#' Build column labels for the calibration parameter vector.
#' @param pc the p_cal environment (post-preprocessing; Bayesian fields optional)
#' @param target desired length; defaults to npars if the Bayesian padding is
#'        present, else nmpars. Pass explicitly to force a length.
#' @return character vector of labels, length == target.
make_param_labels <- function(pc, target = NULL) {
  H      <- .pl_int(pc, "H")
  ncalp  <- .pl_int(pc, "ncalp")
  nsrc   <- .pl_int(pc, "nsource")
  nev    <- isTRUE(.pl_get(pc, "nev"))
  ntheta0<- .pl_int(pc, "ntheta0")
  p_fgsn <- .pl_int(pc, "p_fgsn")   # 0 unless the Bayesian block has run
  p_A    <- .pl_int(pc, "p_A")      # 0 unless the Bayesian block has run
  nmpars <- .pl_int(pc, "nmpars")
  npars  <- .pl_int(pc, "npars")

  if (is.null(target)) target <- if (npars > 0) npars else nmpars
  if (target <= 0) return(character(0))

  lab <- character(0)
  add <- function(x) lab <<- c(lab, x)

  # [new-event theta0] -- only in the joint/new-event layout
  if (nev && ntheta0 > 0) {
    nm <- if (!is.null(pc$theta_names) && length(pc$theta_names) == ntheta0)
            pc$theta_names else paste0("theta0_", seq_len(ntheta0))
    add(paste0("NewEvent | ", nm))
  }

  # [calibration inference parameters]
  if (ncalp > 0) {
    nm <- if (!is.null(pc$cal_par_names) && length(pc$cal_par_names) == ncalp)
            pc$cal_par_names else paste0("calp_", seq_len(ncalp))
    add(paste0("CalParam | ", nm))
  }

  # [errors-in-variables source yields]
  if (nsrc > 0) {
    nm <- if (!is.null(pc$eiv_sources) && length(pc$eiv_sources) == nsrc)
            pc$eiv_sources else seq_len(nsrc)
    add(paste0("EIV yield | source ", nm))
  }

  # [common coefficients]  per phenomenology, per response
  for (hh in seq_len(H)) {
    pb <- pc$h[[hh]]$pbeta
    if (!is.null(pb)) for (rr in seq_along(pb)) {
      if (pb[rr] > 0) add(sprintf("Phen %d | Common coef | Response %d | coef %d",
                                  hh, rr, seq_len(pb[rr])))
    }
  }

  # [emplacement-condition coefficients] per phen, per emplacement, per response
  for (hh in seq_len(H)) {
    pbt <- pc$h[[hh]]$pbetat
    Th  <- .pl_int(pc$h[[hh]], "Th")
    Rh  <- .pl_int(pc$h[[hh]], "Rh")
    if (!is.null(pbt) && Th > 0) for (tt in seq_len(Th)) {
      ptt <- pbt[[tt]]
      if (!is.null(ptt)) for (rr in seq_len(Rh)) {
        if (!is.na(ptt[rr]) && ptt[rr] > 0)
          add(sprintf("Phen %d | Emplacement %d | Response %d | coef %d",
                      hh, tt, rr, seq_len(ptt[rr])))
      }
    }
  }

  # [source variance components]
  for (hh in seq_len(H)) {
    v <- pc$h[[hh]]$pvc_1
    if (!is.null(v)) for (rr in seq_along(v)) {
      if (v[rr] > 0) add(sprintf("Phen %d | Source variance | Response %d | comp %d",
                                 hh, rr, seq_len(v[rr])))
    }
  }

  # [path variance components]
  for (hh in seq_len(H)) {
    v <- pc$h[[hh]]$pvc_2
    if (!is.null(v)) for (rr in seq_along(v)) {
      if (v[rr] > 0) add(sprintf("Phen %d | Path variance | Response %d | comp %d",
                                 hh, rr, seq_len(v[rr])))
    }
  }

  # [observational error covariance]  Rh*(Rh+1)/2 params per phenomenology
  for (hh in seq_len(H)) {
    Rh <- .pl_int(pc$h[[hh]], "Rh")
    if (Rh > 0) {
      ncpar <- Rh * (Rh + 1) / 2
      add(sprintf("Phen %d | Obs-error cov | param %d", hh, seq_len(ncpar)))
    }
  }

  # [FGSN errors-in-variables prior parameters]  (Bayesian only)
  if (p_fgsn > 0) {
    fnm <- c("Alpha", "LambdaSq")
    if (p_fgsn >= 1) add("FGSN prior | Alpha")
    if (p_fgsn >= 2) add("FGSN prior | LambdaSq")
    if (p_fgsn > 2)  add(paste0("FGSN prior | Omega ", seq_len(p_fgsn - 2)))
  }

  # [variance-component prior scale params A]  (Bayesian only). print_sumstats.r
  # emits one scale per phenomenology x response where the source (pvc_1) or path
  # (pvc_2) variance component is present, labeled by phenomenology/response --
  # mirror that enumeration exactly so labels align with the sampled order.
  if (p_A > 0) {
    pvc1_on <- .pl_int(pc, "pvc_1") > 0
    pvc2_on <- .pl_int(pc, "pvc_2") > 0
    for (hh in seq_len(H)) {
      v1 <- pc$h[[hh]]$pvc_1; v2 <- pc$h[[hh]]$pvc_2
      Rh <- .pl_int(pc$h[[hh]], "Rh")
      for (rr in seq_len(Rh)) {
        on1 <- pvc1_on && !is.null(v1) && !is.na(v1[rr]) && v1[rr] > 0
        on2 <- pvc2_on && !is.null(v2) && !is.na(v2[rr]) && v2[rr] > 0
        if (on1 || on2)
          add(sprintf("Phen %d | VC prior scale | Response %d", hh, rr))
      }
    }
  }

  .pl_fit(lab, target)
}

#' New-event layout variant (print_sumstats_0.r ordering). The new-event vector
#' leads with theta0 then calp; it has no common/emplacement/VC/FGSN model
#' blocks of its own. Falls back to make_param_labels for the shared tail if the
#' structure matches.
make_param_labels_0 <- function(pc, target = NULL) {
  ntheta0 <- .pl_int(pc, "ntheta0")
  ncalp   <- .pl_int(pc, "ncalp")
  npars   <- .pl_int(pc, "npars")
  nmpars  <- .pl_int(pc, "nmpars")
  if (is.null(target)) target <- if (npars > 0) npars else
                                   if (nmpars > 0) nmpars else ntheta0 + ncalp
  if (target <= 0) return(character(0))

  lab <- character(0)
  if (ntheta0 > 0) {
    nm <- if (!is.null(pc$theta_names) && length(pc$theta_names) == ntheta0)
            pc$theta_names else paste0("theta0_", seq_len(ntheta0))
    lab <- c(lab, paste0("NewEvent | ", nm))
  }
  if (ncalp > 0) {
    nm <- if (!is.null(pc$cal_par_names) && length(pc$cal_par_names) == ncalp)
            pc$cal_par_names else paste0("calp_", seq_len(ncalp))
    lab <- c(lab, paste0("CalParam | ", nm))
  }
  .pl_fit(lab, target)
}
