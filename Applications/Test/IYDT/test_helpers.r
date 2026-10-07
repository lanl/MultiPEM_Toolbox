if (!requireNamespace("numDeriv", quietly = TRUE)) {
  stop("The numDeriv package is required", call. = FALSE)
}

assert_jacobian_close <- function(actual, expected, label = "Jacobian",
                                  atol = 1e-6, rtol = 1e-6)
{
  actual <- as.matrix(actual)
  expected <- as.matrix(expected)

  if (!identical(dim(actual), dim(expected))) {
    stop(sprintf(
      "%s dimensions differ: analytical %s; numerical %s",
      label, paste(dim(actual), collapse = " x "),
      paste(dim(expected), collapse = " x ")
    ), call. = FALSE)
  }
  if (!all(is.finite(actual)) || !all(is.finite(expected))) {
    stop(sprintf("%s contains a non-finite value", label), call. = FALSE)
  }

  error <- max(abs(actual - expected))
  scale <- max(1, max(abs(expected)))
  tolerance <- atol + rtol * scale
  if (error > tolerance) {
    stop(sprintf(
      "%s mismatch: maximum absolute error %.8g exceeds %.8g",
      label, error, tolerance
    ), call. = FALSE)
  }
  message(sprintf("%s passed (maximum absolute error %.8g)", label, error))
  invisible(error)
}

combine_jacobian <- function(value)
{
  if (!is.list(value)) return(as.matrix(value))
  fields <- intersect(c("jbeta", "jcalp", "jtheta"), names(value))
  if (!length(fields)) stop("Jacobian list has no recognized blocks", call. = FALSE)
  do.call(cbind, value[fields])
}

assert_close <- function(actual, expected, label, atol = 1e-8, rtol = 1e-8)
{
  if (!identical(dim(actual), dim(expected)) || length(actual) != length(expected)) {
    stop(sprintf("%s dimensions differ", label), call. = FALSE)
  }
  if (!all(is.finite(actual)) || !all(is.finite(expected))) {
    stop(sprintf("%s contains a non-finite value", label), call. = FALSE)
  }
  error <- max(abs(actual - expected))
  tolerance <- atol + rtol * max(1, max(abs(expected)))
  if (error > tolerance) {
    stop(sprintf("%s: error %.8g exceeds %.8g", label, error, tolerance),
         call. = FALSE)
  }
  invisible(error)
}

check_jacobian <- function(label, forward, gradient, zeta, params,
                           atol = 1e-6, rtol = 1e-6)
{
  value <- forward(zeta, params)
  analytical <- combine_jacobian(gradient(zeta, params))
  numerical <- numDeriv::jacobian(
    forward, zeta, method.args = list(r = 6), params = params
  )
  if (length(value) != nrow(analytical) || ncol(analytical) != length(zeta)) {
    stop(sprintf("%s returned inconsistent dimensions", label), call. = FALSE)
  }
  if (!is.null(names(value))) {
    stop(sprintf("%s forward result unexpectedly has names", label), call. = FALSE)
  }
  assert_close(analytical, numerical, label, atol, rtol)
  message(sprintf("%s passed", label))
  invisible(list(value = value, jacobian = analytical))
}

check_checkpoint_reuse <- function(label, calibration_state, event_function)
{
  work <- tempfile("multipem-checkpoint-")
  dir.create(work)
  on.exit(unlink(work, recursive = TRUE, force = TRUE), add = TRUE)
  pristine <- file.path(work, "calibration-pristine.RData")
  event_copy <- file.path(work, "event-working.RData")
  save(calibration_state, file = pristine, version = 2)
  pristine_hash <- unname(tools::md5sum(pristine))
  if (!file.copy(pristine, event_copy)) stop("checkpoint copy failed", call. = FALSE)

  event_workspace <- new.env(parent = emptyenv())
  load(event_copy, envir = event_workspace)
  event_workspace$event_result <- event_function(event_workspace$calibration_state)
  if (!length(event_workspace$event_result) ||
      !all(is.finite(event_workspace$event_result))) {
    stop(sprintf("%s produced an invalid event result", label), call. = FALSE)
  }
  save(list = ls(event_workspace), envir = event_workspace,
       file = event_copy, version = 2)

  if (!identical(unname(tools::md5sum(pristine)), pristine_hash)) {
    stop(sprintf("%s changed the pristine checkpoint", label), call. = FALSE)
  }
  if (identical(unname(tools::md5sum(event_copy)), pristine_hash)) {
    stop(sprintf("%s did not update the event copy", label), call. = FALSE)
  }
  pristine_workspace <- new.env(parent = emptyenv())
  load(pristine, envir = pristine_workspace)
  if (exists("event_result", envir = pristine_workspace, inherits = FALSE)) {
    stop(sprintf("%s contaminated the pristine checkpoint", label), call. = FALSE)
  }
  message(sprintf("%s passed", label))
  invisible(event_workspace$event_result)
}
