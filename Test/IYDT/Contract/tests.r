adir <- "../Code"
tdir <- "../Test"
source(file.path(adir, "forward.r"), local = TRUE)
message("***** IYDT invalid-input contract tests *****")
source(file.path(tdir, "test_contract.r"), local = TRUE)
