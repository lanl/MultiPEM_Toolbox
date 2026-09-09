########################################################################
# helpers.r
#
# Secondary code for the NUTS sampler (see nuts.r for the main entry
# point). This file:
#   * compiles and caches the C++ core (nuts_core.cpp), transparently
#     working around a macOS Command Line Tools issue where libc++ headers
#     (e.g. <cmath>) are not found;
#   * provides step-size dual averaging (Nesterov) updates,
#   * provides Welford online estimation of a diagonal mass matrix,
#   * wires up numerical gradients (numDeriv Richardson, or a vectorized
#     central difference) when the user supplies no gradient,
#   * exposes utilities for building compiled (XPtr) targets.
#
# Algorithms follow Hoffman & Gelman (2014) and Betancourt (2017); the
# warmup/adaptation schedule mirrors Stan's windowed adaptation.
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

# ---------------------------------------------------------------------------
# Compilation of the C++ core, with an automatic macOS toolchain workaround.
#
# Newer Command Line Tools ship a libc++ include directory that clang does not
# search by default, so Rcpp compilation fails with "'cmath' file not found".
# The SDK's own libc++ headers work; we add them via -isystem. We only do this
# when a probe actually fails, so the code stays portable.
# ---------------------------------------------------------------------------

# Directory containing this file, captured at source() time. Used to locate
# nuts_core.cpp regardless of the current working directory. When sourced,
# a source() frame on the call stack carries $ofile with this file's path.
# We scan the whole stack for it rather than assuming sys.frame(1): when the
# file is sourced with local=TRUE (e.g. from inside a function), sys.frame(1)
# is an unrelated outer frame with no $ofile. If no frame has $ofile (e.g. the
# code was pasted rather than sourced), fall back to the working directory.
.nuts_helpers_dir <- tryCatch({
  frames <- sys.frames()
  ofile <- NULL
  for (fr in rev(frames)) {
    of <- fr$ofile
    if (!is.null(of) && nzchar(of)) { ofile <- of; break }
  }
  if (is.null(ofile)) getwd()
  else normalizePath(dirname(ofile), mustWork = FALSE)
}, error = function(e) getwd())

# Parented to baseenv() (not emptyenv()) because Rcpp::sourceCpp(env=) evaluates
# its generated DLL-loader code inside this environment and needs `<-`, dyn.load,
# etc. to be in scope.
.nuts_env <- new.env(parent = baseenv())

# Run `R CMD config <what>` (optionally under a chosen DEVELOPER_DIR) and
# return the trimmed value, or NA if the query fails. A misconfigured active
# developer directory (e.g. xcode-select pointing at a broken Xcode.app) makes
# `make`/`xcrun` error out, so R CMD config emits *error text* on stdout with a
# nonzero exit status. We detect both -- a "status" attribute and tell-tale
# error strings -- and report NA so callers can fall back to another toolchain.
.r_cmd_config <- function(what, dev_dir = NULL) {
  old <- Sys.getenv("DEVELOPER_DIR", unset = NA)
  if (!is.null(dev_dir)) {
    Sys.setenv(DEVELOPER_DIR = dev_dir)
    on.exit(if (is.na(old)) Sys.unsetenv("DEVELOPER_DIR")
            else Sys.setenv(DEVELOPER_DIR = old), add = TRUE)
  }
  out <- tryCatch(
    suppressWarnings(system2("R", c("CMD", "config", what),
                             stdout = TRUE, stderr = TRUE)),
    error = function(e) NULL)
  if (is.null(out) || !is.null(attr(out, "status"))) return(NA_character_)
  out <- out[nzchar(out)]
  if (length(out) == 0) return(NA_character_)
  # Guard against error text leaking through as if it were the config value.
  if (any(grepl("error:|couldn't spawn|xcode-select|xcrun:",
                out, ignore.case = TRUE))) return(NA_character_)
  out[1]
}

.cmath_probe_ok <- function(dev_dir = NULL) {
  # Probe with a direct compile rather than Rcpp::sourceCpp: it is faster and,
  # crucially, gives us fully-controlled error handling (a failed sourceCpp can
  # emit an error that is awkward to trap). We ask the same compiler R uses
  # whether it can find <cmath> with R's own flags.
  cxx <- .r_cmd_config("CXX", dev_dir)
  if (is.na(cxx)) return(FALSE)        # toolchain query failed -> not usable
  src <- tempfile(fileext = ".cpp")
  writeLines(c("#include <cmath>", "int main(){ return 0; }"), src)
  on.exit(unlink(src), add = TRUE)
  old <- Sys.getenv("DEVELOPER_DIR", unset = NA)
  if (!is.null(dev_dir)) {
    Sys.setenv(DEVELOPER_DIR = dev_dir)
    on.exit(if (is.na(old)) Sys.unsetenv("DEVELOPER_DIR")
            else Sys.setenv(DEVELOPER_DIR = old), add = TRUE)
  }
  # Split "clang++ -arch x86_64 -std=gnu++17" into command + args.
  toks <- strsplit(cxx, "\\s+")[[1]]
  status <- tryCatch(
    suppressWarnings(system2(toks[1],
      c(toks[-1], "-fsyntax-only", "-x", "c++", src),
      stdout = FALSE, stderr = FALSE)),
    error = function(e) 1L)
  identical(as.integer(status), 0L)
}

.macos_sdk_cxx_include <- function(dev_dir = NULL) {
  old <- Sys.getenv("DEVELOPER_DIR", unset = NA)
  if (!is.null(dev_dir)) {
    Sys.setenv(DEVELOPER_DIR = dev_dir)
    on.exit(if (is.na(old)) Sys.unsetenv("DEVELOPER_DIR")
            else Sys.setenv(DEVELOPER_DIR = old), add = TRUE)
  }
  sdk <- tryCatch(
    suppressWarnings(system2("xcrun", "--show-sdk-path",
                             stdout = TRUE, stderr = TRUE)),
    error = function(e) NULL)
  if (is.null(sdk) || !is.null(attr(sdk, "status"))) return(NULL)
  sdk <- sdk[nzchar(sdk)]
  if (length(sdk) == 0 || any(grepl("error:", sdk, ignore.case = TRUE)))
    return(NULL)
  inc <- file.path(sdk[1], "usr", "include", "c++", "v1")
  if (dir.exists(inc)) inc else NULL
}

# Find a macOS developer directory whose toolchain actually resolves. Tried in
# preference order: the current DEVELOPER_DIR / xcode-select setting, the
# Command Line Tools, then a full Xcode install. Returns NULL if none work.
.macos_working_developer_dir <- function() {
  cur <- tryCatch(
    suppressWarnings(system2("xcode-select", "-p",
                             stdout = TRUE, stderr = FALSE)),
    error = function(e) character(0))
  cands <- unique(c(
    Sys.getenv("DEVELOPER_DIR", unset = ""),
    cur,
    "/Library/Developer/CommandLineTools",
    "/Applications/Xcode.app/Contents/Developer"))
  cands <- cands[nzchar(cands) & dir.exists(cands)]
  for (dd in cands) if (!is.na(.r_cmd_config("CXX", dd))) return(dd)
  NULL
}

# Locate a usable Rtools toolchain on Windows and return the bin directories to
# prepend to PATH (make/gcc/g++ live here). Rtools is R's official Windows
# build toolchain; when it is installed but not on PATH, Rcpp compilation fails
# with "make: command not found" or similar. Returns NULL if none is found.
.windows_rtools_paths <- function() {
  # Prefer pkgbuild if available -- it knows how to locate the matching Rtools.
  if (requireNamespace("pkgbuild", quietly = TRUE)) {
    p <- tryCatch(pkgbuild::rtools_path(), error = function(e) character(0))
    p <- p[nzchar(p) & dir.exists(p)]
    if (length(p)) return(p)
  }
  # Environment hints set by the Rtools installers, newest first.
  roots <- c(Sys.getenv("RTOOLS45_HOME"), Sys.getenv("RTOOLS44_HOME"),
             Sys.getenv("RTOOLS43_HOME"), Sys.getenv("RTOOLS42_HOME"),
             Sys.getenv("RTOOLS40_HOME"))
  roots <- roots[nzchar(roots)]
  if (!length(roots)) {
    # Common install locations as a last resort.
    guesses <- c("C:/rtools45", "C:/rtools44", "C:/rtools43", "C:/rtools42",
                 "C:/rtools40", "C:/Rtools")
    roots <- guesses[dir.exists(guesses)]
  }
  roots <- roots[dir.exists(roots)]
  if (!length(roots)) return(NULL)
  root <- roots[1]
  # coreutils/make live in usr/bin; the mingw compiler in the arch bin dir
  # (name varies across Rtools versions), so include every candidate present.
  cands <- c(file.path(root, "usr", "bin"),
             file.path(root, "x86_64-w64-mingw32.static.posix", "bin"),
             file.path(root, "mingw64", "bin"),
             file.path(root, "mingw32", "bin"),
             file.path(root, "bin"))
  cands <- cands[dir.exists(cands)]
  if (!length(cands)) return(NULL)
  cands
}

# Write a temporary Makevars that preserves the user's ~/.R/Makevars and appends
# the given line(s); returns the temp file path (caller unlinks it).
.append_makevars <- function(extra_lines) {
  mkfile <- tempfile(fileext = ".mk")
  base_mkv <- path.expand("~/.R/Makevars")
  lines <- if (file.exists(base_mkv)) readLines(base_mkv) else character(0)
  writeLines(c(lines, extra_lines), mkfile)
  mkfile
}

# Compile a .cpp file with Rcpp, transparently repairing a broken host
# toolchain first. macOS: if the active developer directory does not resolve,
# switch DEVELOPER_DIR to one that does, then apply the SDK libc++ (-isystem)
# fix so <cmath> et al. are found. Windows: if the compiler probe fails, prepend
# an installed Rtools to PATH. All environment changes are scoped to this call
# and restored on exit. Fails fast with an actionable message when no usable
# toolchain can be found, instead of letting sourceCpp emit a cryptic error.
.rcpp_source_fixed <- function(src, env, verbose = FALSE, rebuild = FALSE) {
  sysname <- Sys.info()[["sysname"]]

  saved <- list()
  set_env <- function(name, value) {
    if (!name %in% names(saved)) saved[[name]] <<- Sys.getenv(name, unset = NA)
    args <- list(value); names(args) <- name
    do.call(Sys.setenv, args)
  }
  mkfile <- NULL
  on.exit({
    for (nm in names(saved)) {
      v <- saved[[nm]]
      if (is.na(v)) Sys.unsetenv(nm)
      else { a <- list(v); names(a) <- nm; do.call(Sys.setenv, a) }
    }
    if (!is.null(mkfile)) unlink(mkfile)
  }, add = TRUE)

  if (sysname == "Darwin" && !.cmath_probe_ok()) {
    dev_dir <- NULL
    if (is.na(.r_cmd_config("CXX"))) {
      # The active developer directory itself is broken -- find a working one.
      dev_dir <- .macos_working_developer_dir()
      if (is.null(dev_dir))
        stop(".rcpp_source_fixed: no working C/C++ toolchain found. Run\n",
             "  sudo xcode-select --switch /Library/Developer/CommandLineTools\n",
             "or install the tools with `xcode-select --install`.",
             call. = FALSE)
      set_env("DEVELOPER_DIR", dev_dir)
      if (verbose) message(".rcpp_source_fixed: using DEVELOPER_DIR=", dev_dir)
    }
    # Add the SDK's libc++ headers (needs a working dev dir for xcrun).
    inc <- .macos_sdk_cxx_include(dev_dir)
    if (!is.null(inc)) {
      mkfile <- .append_makevars(sprintf("CPPFLAGS += -isystem %s", inc))
      set_env("R_MAKEVARS_USER", mkfile)
      if (verbose)
        message(".rcpp_source_fixed: applying macOS libc++ include fix: ", inc)
    } else if (is.null(dev_dir)) {
      warning(".rcpp_source_fixed: C++ probe failed and no SDK libc++ headers ",
              "found; compilation may fail. Consider `xcode-select --install`.")
    }
  } else if (sysname == "Windows" && !.cmath_probe_ok()) {
    bins <- .windows_rtools_paths()
    if (!is.null(bins)) {
      set_env("PATH", paste(c(bins, Sys.getenv("PATH")),
                            collapse = .Platform$path.sep))
      if (verbose)
        message(".rcpp_source_fixed: prepended Rtools to PATH: ",
                paste(bins, collapse = "; "))
    } else {
      warning(".rcpp_source_fixed: C++ probe failed and Rtools not found; ",
              "compilation may fail. Install Rtools from ",
              "https://cran.r-project.org/bin/windows/Rtools/ and ensure it is ",
              "on PATH.")
    }
  }

  Rcpp::sourceCpp(src, env = env, verbose = verbose, rebuild = rebuild)
}

#' Source an arbitrary C++ file with Rcpp, applying the macOS libc++ include
#' fix when needed. Use this to compile XPtr targets etc.
#' @param src Path to a .cpp file.
#' @param env Environment to expose exported functions in (default global).
source_cpp_fixed <- function(src, env = globalenv(), verbose = FALSE) {
  .rcpp_source_fixed(src, env = env, verbose = verbose)
}

#' Compile (once) and cache the NUTS C++ core.
#' @param src Path to nuts_core.cpp (defaults to alongside this file's dir).
#' @param force Recompile even if already cached.
#' @return Invisibly TRUE. Exported C++ functions become available in .nuts_env.
compile_nuts_core <- function(src = NULL, force = FALSE, verbose = FALSE) {
  if (!force && isTRUE(.nuts_env$compiled)) return(invisible(TRUE))

  if (is.null(src)) {
    # Look next to this file first (so the code works from any working
    # directory), then fall back to the working directory.
    cand <- c(file.path(.nuts_helpers_dir, "nuts_core.cpp"),
              "nuts_core.cpp", file.path(getwd(), "nuts_core.cpp"))
    src <- cand[file.exists(cand)][1]
    if (is.na(src)) stop("compile_nuts_core: cannot locate nuts_core.cpp; pass src=")
  }

  .rcpp_source_fixed(src, env = .nuts_env, verbose = verbose, rebuild = force)
  .nuts_env$compiled <- TRUE
  invisible(TRUE)
}

# ---------------------------------------------------------------------------
# Numerical gradients (used only when the user supplies no analytic gradient).
# ---------------------------------------------------------------------------

.warned_numgrad <- local({
  done <- FALSE
  function() { if (!done) { done <<- TRUE
    warning("NUTS: no analytic gradient supplied; using numerical gradients ",
            "(slower and less accurate). Supply grad_f for best performance.",
            call. = FALSE) } }
})

#' Build a gradient function from a log-density f.
#' @param f log-density function taking a numeric vector.
#' @param method "richardson" (numDeriv, accurate) or "simple" (central diff).
#' @param h step for the simple central-difference method.
make_numeric_grad <- function(f, method = c("richardson", "simple"), h = 1e-5) {
  method <- match.arg(method)
  .warned_numgrad()
  if (method == "richardson") {
    if (!requireNamespace("numDeriv", quietly = TRUE))
      stop("numeric_grad='richardson' needs the numDeriv package.")
    function(theta) numDeriv::grad(f, theta)
  } else {
    # Vectorized central difference: 2*d evaluations of f.
    function(theta) {
      d <- length(theta)
      g <- numeric(d)
      for (i in seq_len(d)) {
        tp <- theta; tm <- theta
        step <- h * max(1, abs(theta[i]))
        tp[i] <- tp[i] + step; tm[i] <- tm[i] - step
        g[i] <- (f(tp) - f(tm)) / (2 * step)
      }
      g
    }
  }
}

# ---------------------------------------------------------------------------
# Dual averaging for step-size adaptation (Nesterov; Hoffman & Gelman Alg. 5/6).
# ---------------------------------------------------------------------------

dual_avg_init <- function(eps, delta = 0.8, gamma = 0.05, t0 = 10, kappa = 0.75) {
  list(mu = log(10 * eps), eps = eps, eps_bar = 1, H_bar = 0,
       delta = delta, gamma = gamma, t0 = t0, kappa = kappa, m = 0)
}

#' One dual-averaging update given the iteration's mean Metropolis accept stat.
dual_avg_update <- function(state, accept_stat) {
  if (!is.finite(accept_stat)) accept_stat <- 0
  accept_stat <- min(accept_stat, 1)
  state$m <- state$m + 1
  m <- state$m
  w <- 1 / (m + state$t0)
  state$H_bar <- (1 - w) * state$H_bar + w * (state$delta - accept_stat)
  log_eps <- state$mu - sqrt(m) / state$gamma * state$H_bar
  eta <- m^(-state$kappa)
  state$eps_bar <- exp(eta * log_eps + (1 - eta) * log(state$eps_bar))
  state$eps <- exp(log_eps)
  state
}

# ---------------------------------------------------------------------------
# Welford online mean/variance for diagonal mass-matrix (metric) adaptation.
# We estimate Var(theta) over a warmup window; the metric is M = diag(1/var),
# so inv_M = var. Stan-style regularization shrinks toward 1 for small samples.
# ---------------------------------------------------------------------------

welford_init <- function(d) list(n = 0, mean = numeric(d), m2 = numeric(d))

welford_update <- function(w, x) {
  w$n <- w$n + 1
  delta <- x - w$mean
  w$mean <- w$mean + delta / w$n
  w$m2 <- w$m2 + delta * (x - w$mean)
  w
}

#' Finalize a Welford accumulator to a regularized variance vector.
#' Returns inv_M (the diagonal of M^{-1}), i.e. the estimated variances.
welford_variance <- function(w) {
  if (w$n < 2) return(rep(1, length(w$mean)))
  var <- w$m2 / (w$n - 1)
  n <- w$n
  # Stan's regularization toward a unit metric.
  reg <- (n / (n + 5)) * var + 1e-3 * (5 / (n + 5))
  reg
}

# ---------------------------------------------------------------------------
# Windowed warmup schedule (Stan-style): an initial fast interval for step-size
# only, a sequence of expanding "slow" windows for metric estimation, and a
# final fast interval to finalize the step size. Returns a logical/marker plan.
# ---------------------------------------------------------------------------

#' Build the warmup adaptation schedule.
#' @param warmup total warmup iterations.
#' @param init_buffer iterations of step-size-only adaptation at the start.
#' @param term_buffer iterations of step-size-only adaptation at the end.
#' @param base_window first slow-window size (doubles each window).
#' @return list with vectors: is_slow (metric accumulates), window_end (metric
#'         update boundaries, TRUE on the last iter of each slow window).
make_warmup_schedule <- function(warmup, init_buffer = 75, term_buffer = 50,
                                  base_window = 25) {
  is_slow <- rep(FALSE, warmup)
  window_end <- rep(FALSE, warmup)
  if (warmup == 0) return(list(is_slow = is_slow, window_end = window_end))

  # Shrink buffers if warmup is small (mirror Stan's fallback proportions).
  if (init_buffer + term_buffer + base_window > warmup) {
    init_buffer <- max(1, floor(0.15 * warmup))
    term_buffer <- max(1, floor(0.10 * warmup))
    base_window <- max(1, warmup - init_buffer - term_buffer)
  }

  slow_start <- init_buffer + 1
  slow_end <- warmup - term_buffer
  if (slow_end >= slow_start) {
    is_slow[slow_start:slow_end] <- TRUE
    # Expanding windows: base_window, then doubling, last window absorbs remainder.
    pos <- slow_start
    win <- base_window
    while (pos <= slow_end) {
      nxt <- min(pos + win - 1, slow_end)
      # If the following window would overshoot, extend current to the end.
      if (nxt + 2 * win > slow_end) nxt <- slow_end
      window_end[nxt] <- TRUE
      pos <- nxt + 1
      win <- win * 2
    }
  }
  list(is_slow = is_slow, window_end = window_end)
}

# ---------------------------------------------------------------------------
# find_reasonable_epsilon wrapper (delegates to the compiled core).
# ---------------------------------------------------------------------------
find_reasonable_epsilon <- function(theta, f, grad_f, inv_M, eps = 1) {
  compile_nuts_core()
  .nuts_env$find_reasonable_epsilon_cpp(theta, f, grad_f, inv_M, eps)
}
