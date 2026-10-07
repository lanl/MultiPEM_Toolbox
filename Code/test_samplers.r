local({
  required <- c("adaptMCMC", "FME")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing sampler package(s): ", paste(missing, collapse = ", "),
         call. = FALSE)
  }

  source(file.path(tdir, "helpers.r"), local = TRUE)
  source(file.path(tdir, "nuts.r"), local = TRUE)
  compile_nuts_core(file.path(tdir, "nuts_core.cpp"))

  target_mean <- c(0.75, -1.25)
  target_covariance <- matrix(c(1, 0.35, 0.35, 0.5), nrow = 2L)
  target_precision <- solve(target_covariance)
  log_target <- function(value) {
    delta <- value - target_mean
    -0.5 * drop(crossprod(delta, target_precision %*% delta))
  }
  gradient_target <- function(value) {
    -drop(target_precision %*% (value - target_mean))
  }
  fused_target <- function(value) {
    list(logp = log_target(value), grad = gradient_target(value))
  }

  assert_sample <- function(draws, label, mean_tolerance = 0.12,
                            covariance_tolerance = 0.15) {
    draws <- as.matrix(draws)
    if (nrow(draws) < 1000L || ncol(draws) != 2L || !all(is.finite(draws))) {
      stop(label, " did not return a finite two-parameter sample", call. = FALSE)
    }
    mean_error <- max(abs(colMeans(draws) - target_mean))
    covariance_error <- max(abs(stats::cov(draws) - target_covariance))
    if (mean_error > mean_tolerance || covariance_error > covariance_tolerance) {
      stop(sprintf(
        "%s posterior errors exceed tolerance (mean %.6g; covariance %.6g)",
        label, mean_error, covariance_error
      ), call. = FALSE)
    }
    message(sprintf(
      "%s posterior assertions passed (mean error %.6g; covariance error %.6g)",
      label, mean_error, covariance_error
    ))
  }

  run_ram <- function() {
    set.seed(414)
    adaptMCMC::MCMC(
      log_target, n = 10000, init = c(0, 0), scale = target_covariance,
      acc.rate = 0.234, showProgressBar = FALSE
    )
  }
  ram <- run_ram()
  ram_repeat <- run_ram()
  if (!identical(ram$samples, ram_repeat$samples) ||
      !identical(ram$acceptance.rate, ram_repeat$acceptance.rate)) {
    stop("RAM is not reproducible under a fixed seed", call. = FALSE)
  }
  if (!is.finite(ram$acceptance.rate) || ram$acceptance.rate <= 0 ||
      ram$acceptance.rate >= 1) {
    stop("RAM returned an invalid acceptance rate", call. = FALSE)
  }
  assert_sample(ram$samples[-seq_len(2000L), , drop = FALSE], "RAM")

  run_fme <- function() {
    set.seed(415)
    FME::modMCMC(
      function(value) -2 * log_target(value), p = c(0, 0),
      jump = target_covariance, niter = 10000, burninlength = 2000,
      updatecov = 100, ntrydr = 2, verbose = FALSE
    )
  }
  fme <- run_fme()
  fme_repeat <- run_fme()
  if (!identical(fme$pars, fme_repeat$pars) ||
      !identical(fme$naccepted, fme_repeat$naccepted)) {
    stop("FME is not reproducible under a fixed seed", call. = FALSE)
  }
  if (fme$naccepted <= 0L || fme$naccepted >= 10000L) {
    stop("FME returned an invalid acceptance count", call. = FALSE)
  }
  assert_sample(fme$pars, "FME")

  callback <- NUTS(
    c(0, 0), f = log_target, grad_f = gradient_target,
    n_iter = 1800, warmup = 600, seed = 416,
    diagnostics = TRUE, verbose = FALSE
  )
  fused <- NUTS(
    c(0, 0), fg = fused_target,
    n_iter = 1800, warmup = 600, seed = 416,
    diagnostics = TRUE, verbose = FALSE
  )
  if (!identical(callback$theta, fused$theta) ||
      !identical(callback$diagnostics, fused$diagnostics)) {
    stop("NUTS callback and fused backends diverged under a fixed seed",
         call. = FALSE)
  }
  assert_sample(callback$theta, "NUTS")
  post_warmup <- callback$diagnostics[
    !callback$diagnostics$warmup, , drop = FALSE
  ]
  if (any(post_warmup$divergent) ||
      !all(is.finite(post_warmup$accept_stat)) ||
      any(post_warmup$accept_stat < 0 | post_warmup$accept_stat > 1)) {
    stop("NUTS returned invalid post-warmup diagnostics", call. = FALSE)
  }
  message("Fixed-seed FME, RAM, and NUTS backend checks passed")
})
