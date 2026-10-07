calibration_data <- read.csv(
  file.path(tdir, "data", "seismic_cal.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
true_beta <- c(1.1, -0.65, 0.4, 0.08, -0.3)
calibration_params <- list(
  pbeta = 5, iresp = TRUE, yield_scaling = 1 / 3,
  X = as.matrix(calibration_data[, c("lRange", "W", "C2N", "HOB")]),
  cal = FALSE, cal_par_names = character(), ncalp = 0,
  theta_names = character(), notExp = notExp, dnotExp = dnotExp
)
benchmark <- f_s(true_beta, calibration_params)
objective <- function(intercept) {
  candidate <- c(intercept, true_beta[-1])
  sum((f_s(candidate, calibration_params) - benchmark)^2)
}
fit <- optimize(objective, true_beta[1] + c(-2, 2), tol = 1e-12)
assert_close(fit$minimum, true_beta[1], "IYDT checkpoint calibration", 1e-7, 1e-7)

calibration_state <- list(
  beta = c(fit$minimum, true_beta[-1]),
  benchmark_prediction = benchmark,
  application = "IYDT"
)
event_data <- read.csv(
  file.path(tdir, "data", "seismic_new.csv"),
  check.names = FALSE, stringsAsFactors = FALSE
)
event_function <- function(state) {
  params <- list(
    beta = state$beta, theta_names = c("W", "HOB"), iresp = TRUE,
    yield_scaling = 1 / 3,
    X = as.matrix(event_data[, "lRange", drop = FALSE]), notExp = notExp
  )
  f0_s(unlist(event_data[1, c("W", "HOB")], use.names = FALSE), params)
}
check_checkpoint_reuse(
  "IYDT pristine calibration checkpoint reuse", calibration_state, event_function
)
