tdir <- "../Test"
source(file.path(tdir, "test_bayes_smoke.r"), local = TRUE)
source(file.path(tdir, "test_advanced_integration.r"), local = TRUE)

message("***** IYDT real-application RAM and NUTS tests *****")
case_dir <- file.path("..", "Smoke", "Runfiles", "IYDT-gsrp", "Crater",
                      "I-SUGAR-0")
run_application_backend(case_dir, "RAM", "IYDT Crater I-SUGAR-0")
run_application_backend(case_dir, "NUTS", "IYDT Crater I-SUGAR-0")
