tdir <- "../Test"
source(file.path(tdir, "test_bayes_smoke.r"), local = TRUE)
source(file.path(tdir, "test_advanced_integration.r"), local = TRUE)

message("***** IYDT staged multi-phenomenology integration tests *****")
run_multphen_staged_inputs(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp"),
  "IYDT optical/crater"
)
