tdir <- "../Test"
source(file.path(tdir, "test_bayes_smoke.r"), local = TRUE)
source(file.path(tdir, "test_advanced_integration.r"), local = TRUE)

message("***** IYDT sequential/parallel equivalence tests *****")
run_parallel_equivalence(
  file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater", "I-SUGAR-0"),
  "IYDT Crater I-SUGAR-0"
)
