assert_contains <- function(output, pattern, label)
{
  if (!any(grepl(pattern, output, fixed = TRUE))) {
    stop(sprintf("%s did not contain %s", label, shQuote(pattern)),
         call. = FALSE)
  }
}

calibration_pc <- new.env(parent = emptyenv())
calibration_pc$ncalp <- 0L
calibration_pc$nsource <- 0L
calibration_pc$pbeta <- 0L
calibration_pc$ptbeta <- 0L
calibration_pc$pvc_1 <- 0L
calibration_pc$pvc_2 <- 0L
calibration_pc$p_A <- 0L
calibration_pc$iPrior <- FALSE
calibration_pc$H <- 1L
calibration_pc$h <- list(list(Rh = 1L))

calibration_output <- capture.output(
  calibration_result <- print_ss(log(2), calibration_pc)
)
assert_contains(calibration_output, "OBSERVATIONAL ERROR COVARIANCE PARAMETERS",
                "calibration summary")
assert_contains(calibration_output, "Phenomenology 1", "calibration summary")
assert_contains(calibration_output, "Variances", "calibration summary")
if (!identical(calibration_result, calibration_pc)) {
  stop("print_ss did not return its parameter-control environment", call. = FALSE)
}

event_pc <- new.env(parent = emptyenv())
event_pc$ntheta0 <- 1L
event_pc$theta_names <- "W"
event_pc$iPrior <- TRUE
event_pc$Sigma_mle_0 <- list(acov_0 = 0)
event_pc$ncalp <- 0L

event_output <- capture.output(event_result <- print_ss_0(log(10), event_pc))
assert_contains(event_output, "NEW EVENT INFERENCE PARAMETERS", "event summary")
assert_contains(event_output, "ESTIMATE:", "event summary")
assert_contains(event_output, "W", "event summary")
if (!identical(event_result$tmap_0, setNames(log(10), "W"))) {
  stop("print_ss_0 did not retain the named event estimate", call. = FALSE)
}
message("summary-output contracts passed")
