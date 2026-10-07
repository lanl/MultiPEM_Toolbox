check_prior_gradient <- function(label, log_prior, gradient, value, params)
{
  log_value <- log_prior(value, params)
  if (length(log_value) != 1L || !is.finite(log_value)) {
    stop(sprintf("%s returned an invalid log density", label), call. = FALSE)
  }
  analytical <- as.numeric(gradient(value, params))
  numerical <- as.numeric(numDeriv::grad(log_prior, value, p = params))
  assert_close(
    analytical, numerical, label, atol = 1e-7, rtol = 1e-7
  )
  message(sprintf("%s passed", label))
}

normal_params <- list(
  pi_w_mu = log(20), pi_w_sd = 1.3,
  pi_h_mu = 2.5, pi_h_sd = 18,
  pi_c_mu = log(1.2), pi_c_sd = 0.7
)

for (value in list(c(log(0.8), -10), c(log(25), 12))) {
  check_prior_gradient("IYDT event prior gradient", lp_0, lq_0, value, normal_params)
}
for (value in c(log(0.5), log(3))) {
  check_prior_gradient(
    "IYDT calibration prior gradient", lp_c, lq_c, value, normal_params
  )
  check_prior_gradient("IYDT yield prior gradient", lp_w, lq_w, value, normal_params)
}

transform_params <- list(
  notExp = notExp, dnotExp = dnotExp, d2notExp = d2notExp
)
for (beta3 in c(-2, -0.5, 0.5, 2)) {
  beta <- c(0.4, -0.2, beta3, 0.1, -0.3)
  check_prior_gradient(
    paste("IYDT seismic coefficient prior gradient", beta3),
    lp_s, lq_s, beta, transform_params
  )
}

for (pair in list(c(-2, 0.5), c(-0.5, 2), c(0.5, -2), c(2, -0.5))) {
  beta <- c(-1, 0.25, pair)
  check_prior_gradient(
    paste("IYDT optical coefficient prior gradient", paste(pair, collapse = ",")),
    lp_o, lq_o, beta, transform_params
  )
}
