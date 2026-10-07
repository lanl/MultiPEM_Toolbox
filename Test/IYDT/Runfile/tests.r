tdir <- "../Test"
source(file.path(tdir, "test_runfile_smoke.r"), local = TRUE)
message("***** IYDT real runfile smoke test *****")
run_runfile_smoke(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater", "I-SUGAR-0"),
  "IYDT Crater I-SUGAR-0"
)
