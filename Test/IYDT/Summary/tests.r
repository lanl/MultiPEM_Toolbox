adir <- "../Code"
tdir <- "../Test"
source(file.path(adir, "print_sumstats.r"), local = TRUE)
source(file.path(adir, "print_sumstats_0.r"), local = TRUE)
message("***** IYDT summary-output tests *****")
source(file.path(tdir, "test_summary.r"), local = TRUE)
