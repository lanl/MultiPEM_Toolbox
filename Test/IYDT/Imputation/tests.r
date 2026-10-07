tdir <- "../Test"
source(file.path(tdir, "test_bayes_smoke.r"), local = TRUE)
source(file.path(tdir, "test_imputation.r"), local = TRUE)

message("***** IYDT two-imputation RAM/SMC test *****")
run_multiple_imputation_smoke(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater", "I-SUGAR-0"),
  "IYDT Crater I-SUGAR-0", c(10), c(18)
)
