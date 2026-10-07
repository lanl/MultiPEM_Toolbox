adir <- "../Code"
tdir <- "../Test"

source(file.path(adir, "forward.r"), local = TRUE)
source(file.path(adir, "jacobian.r"), local = TRUE)
source(file.path(adir, "transform.r"), local = TRUE)
source(file.path(tdir, "test_helpers.r"), local = TRUE)

message("***** IYDT transform tests *****")
source(file.path(tdir, "test_transform.r"), local = TRUE)
