tdir <- "../Test"
source(file.path(tdir, "test_bayes_smoke.r"), local = TRUE)
message("***** IYDT Bayesian and shared-inference tests *****")
run_bayesian_smoke(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater", "I-SUGAR-0"),
  "IYDT Crater I-SUGAR-0"
)
run_eiv_derivative_smoke(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater",
            "I-EIV-SUGAR-0"),
  "IYDT Crater I-EIV-SUGAR-0"
)
