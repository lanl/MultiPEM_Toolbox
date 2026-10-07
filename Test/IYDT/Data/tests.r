adir <- "../Code"
tdir <- "../Test"

source(file.path(adir, "forward.r"), local = TRUE)
source(file.path(adir, "forward_0.r"), local = TRUE)
source(file.path(adir, "jacobian.r"), local = TRUE)
source(file.path(adir, "jacobian_0.r"), local = TRUE)
source(file.path(tdir, "test_helpers.r"), local = TRUE)

message("***** IYDT data/configuration tests *****")
source(file.path(tdir, "test_data_config.r"), local = TRUE)
