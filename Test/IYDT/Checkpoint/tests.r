adir <- "../Code"
tdir <- "../Test"

source(file.path(adir, "forward.r"), local = TRUE)
source(file.path(adir, "forward_0.r"), local = TRUE)
source(file.path(adir, "jacobian.r"), local = TRUE)
source(file.path(tdir, "test_helpers.r"), local = TRUE)

message("***** IYDT checkpoint-reuse test *****")
source(file.path(tdir, "test_checkpoint_reuse.r"), local = TRUE)
