adir <- "../Code"
tdir <- "../Test"

source(file.path(adir, "forward.r"), local = TRUE)
source(file.path(adir, "jacobian.r"), local = TRUE)
source(file.path(adir, "lp_0.r"), local = TRUE)
source(file.path(adir, "glp_0.r"), local = TRUE)
source(file.path(adir, "lp_c.r"), local = TRUE)
source(file.path(adir, "glp_c.r"), local = TRUE)
source(file.path(adir, "Phenomenology", "lp_beta_s.r"), local = TRUE)
source(file.path(adir, "Phenomenology", "glp_beta_s.r"), local = TRUE)
source(file.path(adir, "Phenomenology", "lp_beta_o.r"), local = TRUE)
source(file.path(adir, "Phenomenology", "glp_beta_o.r"), local = TRUE)
source(file.path(adir, "Phenomenology", "lp_w.r"), local = TRUE)
source(file.path(adir, "Phenomenology", "glp_w.r"), local = TRUE)
source(file.path(tdir, "test_helpers.r"), local = TRUE)

message("***** IYDT prior tests *****")
source(file.path(tdir, "test_prior.r"), local = TRUE)
