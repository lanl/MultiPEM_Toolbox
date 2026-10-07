if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("The numDeriv package is required", call. = FALSE)
}

patch_runfile_setting <- function(path, name, value)
{
  lines <- readLines(path, warn = FALSE)
  pattern <- paste0(
    "^[[:space:]]*", name,
    "[[:space:]]*=[[:space:]]*[^,]+[[:space:]]*$"
  )
  hits <- grep(pattern, lines)
  if (length(hits) != 1L) {
    stop(sprintf("%s must assign %s exactly once", basename(path), name),
         call. = FALSE)
  }
  lines[hits] <- paste(name, "=", value)
  writeLines(lines, path, useBytes = TRUE)
}

run_batch_checked <- function(script, transcript, restore)
{
  arguments <- c("CMD", "BATCH", "--no-save")
  if (!restore) arguments <- c(arguments, "--no-restore")
  status <- system2(
    file.path(R.home("bin"), "R"), c(arguments, script, transcript)
  )
  if (!identical(status, 0L)) {
    stop(sprintf("%s failed with status %d; inspect %s",
                 script, status, normalizePath(transcript, mustWork = FALSE)),
         call. = FALSE)
  }
  output <- readLines(transcript, warn = FALSE)
  if (any(grepl("Execution halted|^Error( in|:)", output))) {
    stop(sprintf("%s reported an R error", script), call. = FALSE)
  }
}

load_workspace <- function(path)
{
  workspace <- new.env(parent = globalenv())
  load(path, envir = workspace)
  if (!exists("p_cal", envir = workspace, inherits = FALSE)) {
    stop(sprintf("%s does not contain p_cal", path), call. = FALSE)
  }
  workspace$p_cal
}

assert_finite_sample <- function(value, label, minimum_rows = 10L)
{
  value <- as.matrix(value)
  if (nrow(value) < minimum_rows || !ncol(value) || !all(is.finite(value))) {
    stop(sprintf("%s is not a finite posterior sample", label), call. = FALSE)
  }
  invisible(value)
}

assert_gradient <- function(label, objective, gradient, point, pc,
                            atol = 2e-5, rtol = 2e-5)
{
  point <- as.numeric(point)
  analytical <- as.numeric(gradient(point, pc))
  numerical <- as.numeric(numDeriv::grad(
    function(value) objective(value, pc), point, method.args = list(r = 6)
  ))
  if (length(analytical) != length(point) || !all(is.finite(analytical))) {
    stop(sprintf("%s returned an invalid gradient", label), call. = FALSE)
  }
  finite <- is.finite(numerical)
  if (any(!finite)) {
    boundary <- vapply(which(!finite), function(index) {
      step <- 1e-6 * max(1, abs(point[[index]]))
      lower <- upper <- point
      lower[[index]] <- lower[[index]] - step
      upper[[index]] <- upper[[index]] + step
      values <- c(objective(lower, pc), objective(upper, pc))
      length(values) != 2L || any(!is.finite(values))
    }, logical(1))
    if (!all(boundary)) {
      stop(sprintf("%s returned an invalid numerical gradient away from a boundary",
                   label), call. = FALSE)
    }
    minimum <- max(1L, ceiling(0.8 * length(point)))
    if (sum(finite) < minimum) {
      stop(sprintf("%s has only %d of %d numerically checkable coordinates",
                   label, sum(finite), length(point)), call. = FALSE)
    }
    message(sprintf(
      "%s: omitted boundary coordinate(s) %s from the numerical comparison",
      label, paste(which(!finite), collapse = ", ")
    ))
  }
  error <- max(abs(analytical[finite] - numerical[finite]))
  tolerance <- atol + rtol * max(1, max(abs(numerical[finite])))
  if (error > tolerance) {
    stop(sprintf("%s error %.8g exceeds %.8g", label, error, tolerance),
         call. = FALSE)
  }
  message(sprintf("%s passed (maximum absolute error %.8g)", label, error))
  invisible(error)
}

prepare_mle_scripts <- function(case_dir, bayes)
{
  scripts <- file.path(case_dir, c("runMPEM.r", "runMPEM_0.r"))
  if (any(!file.exists(scripts))) {
    stop(sprintf("runfile pair is incomplete under %s", case_dir), call. = FALSE)
  }
  if (!all(vapply(scripts, function(path) {
    any(grepl("^[[:space:]]*set.seed\\(", readLines(path, warn = FALSE)))
  }, logical(1)))) {
    stop("Bayesian smoke runfiles must set deterministic seeds", call. = FALSE)
  }
  for (script in scripts) {
    patch_runfile_setting(script, "nstart", "1")
    patch_runfile_setting(script, "ncores_mle", "1")
    patch_runfile_setting(script, "mle_grad_ck", "FALSE")
    patch_runfile_setting(script, "iBayes", if (bayes) "TRUE" else "FALSE")
  }
  patch_runfile_setting(scripts[[2L]], "nimpute", "1")
  if (bayes) {
    for (script in scripts) {
      patch_runfile_setting(script, "nburn", "20")
      patch_runfile_setting(script, "nmcmc", "40")
      patch_runfile_setting(script, "nthin", "2")
      patch_runfile_setting(script, "ncores_map", "1")
      patch_runfile_setting(script, "ncores_mc", "1")
      patch_runfile_setting(script, "prior_grad_ck", "FALSE")
      patch_runfile_setting(script, "post_grad_ck", "FALSE")
    }
  }
}

run_bayesian_smoke <- function(case_dir, label)
{
  prepare_mle_scripts(case_dir, bayes = TRUE)
  old <- setwd(case_dir)
  on.exit(setwd(old), add = TRUE)
  unlink(c(".RData", "opt.RData", "opt_nev.RData",
           "calibration-pristine.RData", "bayes-calibration.out",
           "bayes-event.out"), force = TRUE)

  run_batch_checked("runMPEM.r", "bayes-calibration.out", restore = FALSE)
  if (!file.exists(".RData") || !file.exists("opt.RData")) {
    stop(sprintf("%s Bayesian calibration outputs are incomplete", label),
         call. = FALSE)
  }
  calibration <- load_workspace(".RData")
  assert_finite_sample(calibration$mpi, paste(label, "calibration posterior"))
  assert_gradient(
    paste(label, "calibration likelihood gradient"),
    calibration$ll_cal, calibration$gll_cal, calibration$mle_cal, calibration
  )
  assert_gradient(
    paste(label, "calibration posterior gradient"),
    calibration$lpost_full, calibration$glpost_full,
    calibration$map_cal, calibration
  )
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve the Bayesian calibration checkpoint", call. = FALSE)
  }
  pristine_hash <- unname(tools::md5sum("calibration-pristine.RData"))

  run_batch_checked("runMPEM_0.r", "bayes-event.out", restore = TRUE)
  if (!file.exists("opt_nev.RData")) {
    stop(sprintf("%s Bayesian event output is incomplete", label), call. = FALSE)
  }
  event <- load_workspace(".RData")
  assert_finite_sample(event$tmpi_0, paste(label, "event posterior"))
  assert_gradient(
    paste(label, "event likelihood gradient"),
    event$ll_0, event$gll_0, event$mle, event
  )
  assert_gradient(
    paste(label, "event posterior gradient"),
    event$lpost_0, event$glpost_0, event$map, event
  )
  if (!identical(unname(tools::md5sum("calibration-pristine.RData")),
                 pristine_hash)) {
    stop(sprintf("%s changed the pristine Bayesian checkpoint", label),
         call. = FALSE)
  }
  if (identical(unname(tools::md5sum(".RData")), pristine_hash)) {
    stop(sprintf("%s event run did not update its working checkpoint", label),
         call. = FALSE)
  }
  message(sprintf("%s bounded Bayesian calibration/event workflow passed", label))
  invisible(TRUE)
}

run_eiv_derivative_smoke <- function(case_dir, label)
{
  prepare_mle_scripts(case_dir, bayes = FALSE)
  old <- setwd(case_dir)
  on.exit(setwd(old), add = TRUE)
  unlink(c(".RData", "opt.RData", "opt_nev.RData", "eiv-calibration.out",
           "eiv-event.out"), force = TRUE)

  run_batch_checked("runMPEM.r", "eiv-calibration.out", restore = FALSE)
  calibration <- load_workspace(".RData")
  if (!isTRUE(calibration$eiv)) {
    stop(sprintf("%s did not enable errors-in-variables", label), call. = FALSE)
  }
  assert_gradient(
    paste(label, "EIV calibration likelihood gradient"),
    calibration$ll_cal, calibration$gll_cal,
    calibration$mle_cal, calibration, atol = 5e-5, rtol = 5e-5
  )

  run_batch_checked("runMPEM_0.r", "eiv-event.out", restore = TRUE)
  event <- load_workspace(".RData")
  assert_gradient(
    paste(label, "EIV event likelihood gradient"),
    event$ll_0, event$gll_0, event$mle, event,
    atol = 5e-5, rtol = 5e-5
  )
  message(sprintf("%s EIV derivative workflow passed", label))
  invisible(TRUE)
}
