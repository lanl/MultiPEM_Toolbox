tdir <- "../Test"
source(file.path(tdir, "test_bayes_smoke.r"), local = TRUE)
source(file.path(tdir, "test_advanced_integration.r"), local = TRUE)

message("***** IYDT reduced profile-likelihood test *****")
run_profile_likelihood_smoke(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater", "I-SUGAR-0"),
  "IYDT Crater I-SUGAR-0"
)
