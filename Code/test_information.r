local({
  if (!requireNamespace("Matrix", quietly = TRUE) ||
      !requireNamespace("numDeriv", quietly = TRUE)) {
    stop("Matrix and numDeriv are required for information-matrix tests",
         call. = FALSE)
  }
  suppressPackageStartupMessages(library(Matrix))
  source(file.path(tdir, "tryCatch.W.E.r"), local = TRUE)
  source(file.path(tdir, "transform.r"), local = TRUE)
  source(file.path(tdir, "observed_information_cal.r"), local = TRUE)
  source(file.path(tdir, "observed_information_0.r"), local = TRUE)

  target <- c(0.4, -1.1, 0.7)
  precision <- matrix(
    c(3.0, 0.4, -0.2,
      0.4, 2.2, 0.3,
     -0.2, 0.3, 1.6),
    nrow = 3L, byrow = TRUE
  )
  log_likelihood <- function(value) {
    delta <- value - target
    -0.5 * drop(crossprod(delta, precision %*% delta))
  }
  hessian <- numDeriv::hessian(log_likelihood, target)
  if (max(abs(hessian + precision)) > 1e-7) {
    stop("Numerical Hessian does not match the Gaussian information matrix",
         call. = FALSE)
  }
  optimizer <- list(par = target, hessian = hessian)
  expected_covariance <- solve(precision)

  calibration <- new.env(parent = emptyenv())
  calibration$ncalp <- 3L
  calibration$tryCatch.W.E <- tryCatch.W.E
  calibration_result <- obs_info_cal(calibration, optimizer, imle = TRUE)
  if (!identical(calibration_result$acov_cal, 1) ||
      max(abs(as.matrix(calibration_result$II_calp) -
              expected_covariance)) > 1e-7) {
    stop("Calibration observed-information covariance is incorrect",
         call. = FALSE)
  }
  calibration <- obs_info_cal(calibration, optimizer, imle = FALSE)
  if (max(abs(as.matrix(calibration$IHess) - expected_covariance)) > 1e-7) {
    stop("Calibration posterior proposal covariance is incorrect",
         call. = FALSE)
  }

  event <- new.env(parent = emptyenv())
  event$ntheta0 <- 2L
  event$ncalp <- 1L
  event$opt_B <- FALSE
  event$itheta0_bounds <- rep(list(integer()), 3L)
  event$tryCatch.W.E <- tryCatch.W.E
  event_result <- obs_info_0(event, optimizer, imle = TRUE)
  if (!identical(event_result$acov_0, 1) ||
      !identical(event_result$acov_cal, 1) ||
      max(abs(as.matrix(event_result$II_nev_it) -
              expected_covariance[1:2, 1:2])) > 1e-7 ||
      max(abs(as.matrix(event_result$II_nev) -
              expected_covariance[1:2, 1:2])) > 1e-7 ||
      max(abs(as.matrix(event_result$II_calp) -
              expected_covariance[3, 3, drop = FALSE])) > 1e-7) {
    stop("Event observed-information covariance blocks are incorrect",
         call. = FALSE)
  }
  event <- obs_info_0(event, optimizer, imle = FALSE)
  if (max(abs(as.matrix(event$IHess) - expected_covariance)) > 1e-7) {
    stop("Event posterior proposal covariance is incorrect", call. = FALSE)
  }

  bounded_precision <- matrix(
    c(3.4, 0.2, 0.1, -0.1,
      0.2, 2.8, -0.2, 0.3,
      0.1, -0.2, 2.5, 0.2,
     -0.1, 0.3, 0.2, 1.9),
    nrow = 4L, byrow = TRUE
  )
  bounded <- new.env(parent = emptyenv())
  bounded$ntheta0 <- 3L
  bounded$ncalp <- 1L
  bounded$opt_B <- TRUE
  bounded$itheta0_bounds <- list(1L, 2L, 3L)
  bounded$theta0_bounds <- matrix(
    c(1, Inf, -Inf, 5, 2, 10), nrow = 3L, byrow = TRUE
  )
  bounded$theta0_range <- 8
  bounded$notExp <- notExp
  bounded$dnotExp <- dnotExp
  bounded$notLog <- notLog
  bounded$inv_transform <- inv_transform
  bounded$tryCatch.W.E <- tryCatch.W.E
  bounded_coordinates <- c(0.2, -0.3, 0.4)
  bounded_physical <- c(
    bounded$theta0_bounds[1L, 1L] + notExp(bounded_coordinates[1L]),
    bounded$theta0_bounds[2L, 2L] - notExp(bounded_coordinates[2L]),
    bounded$theta0_bounds[3L, 1L] +
      bounded$theta0_range /
      (1 + 1 / notExp(bounded_coordinates[3L]))
  )
  bounded_optimizer <- list(
    par = c(bounded_physical, 0.7), hessian = -bounded_precision
  )
  bounded_tcal <- as.list(bounded)
  bounded_jacobian <- diag(4L)
  bounded_jacobian[1L, 1L] <- dnotExp(bounded_coordinates[1L])
  bounded_jacobian[2L, 2L] <- -dnotExp(bounded_coordinates[2L])
  bounded_tau <- notExp(bounded_coordinates[3L])
  bounded_jacobian[3L, 3L] <-
    dnotExp(bounded_coordinates[3L]) * bounded$theta0_range /
    (1 + bounded_tau)^2
  bounded_expected <- solve(
    t(bounded_jacobian) %*% bounded_precision %*% bounded_jacobian
  )
  bounded_result <- obs_info_0(
    bounded, bounded_optimizer, imle = FALSE, t_cal = bounded_tcal
  )
  if (max(abs(as.matrix(bounded_result$IHess) - bounded_expected)) > 1e-7) {
    stop("Bounded event proposal covariance missed a chain-rule branch",
         call. = FALSE)
  }

  transformed <- new.env(parent = emptyenv())
  transformed$ntheta0 <- 3L
  transformed$ncalp <- 1L
  transformed$opt_B <- FALSE
  transformed$itransform <- TRUE
  transformed$itheta0_bounds <- list(1L, 2L, 3L)
  transformed$theta0_range <- 8
  transformed$notExp <- notExp
  transformed$dnotExp <- dnotExp
  transformed$tryCatch.W.E <- tryCatch.W.E
  user_jacobian <- diag(c(1.2, 0.8, 1.5))
  transformed$j_tau <- function(x, pc) user_jacobian
  transformed$tau <- function(x, pc) drop(user_jacobian %*% x)
  transformed_optimizer <- list(
    par = c(0.15, -0.25, 0.35, -0.4), hessian = -bounded_precision
  )
  transformed_coordinates <- transformed$tau(
    transformed_optimizer$par[1:3], transformed
  )
  output_jacobian <- diag(3L)
  output_jacobian[1L, 1L] <- dnotExp(transformed_coordinates[1L])
  output_jacobian[2L, 2L] <- -dnotExp(transformed_coordinates[2L])
  transformed_tau <- notExp(transformed_coordinates[3L])
  output_jacobian[3L, 3L] <-
    dnotExp(transformed_coordinates[3L]) * transformed$theta0_range /
    (1 + transformed_tau)^2
  base_covariance <- solve(bounded_precision)
  transformed_expected <- output_jacobian %*% user_jacobian %*%
    base_covariance[1:3, 1:3] %*% t(user_jacobian) %*%
    t(output_jacobian)
  transformed_result <- obs_info_0(
    transformed, transformed_optimizer, imle = TRUE
  )
  if (max(abs(as.matrix(transformed_result$II_nev_it) -
              base_covariance[1:3, 1:3])) > 1e-7 ||
      max(abs(as.matrix(transformed_result$II_nev) -
              transformed_expected)) > 1e-7 ||
      max(abs(as.matrix(transformed_result$II_calp) -
              base_covariance[4, 4, drop = FALSE])) > 1e-7) {
    stop("Transformed event covariance propagation is incorrect",
         call. = FALSE)
  }
  message(paste(
    "Calibration, event, bounded, and transformed observed-information",
    "assertions passed"
  ))
})
