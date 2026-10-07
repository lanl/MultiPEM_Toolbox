pc <- list(ntheta0 = 2, tpars = list(yield_scaling = 1 / 3))

for (x in list(c(log(0.5), -12), c(log(25), 14), c(log(1e6), 0))) {
  transformed <- tau(x, pc)
  assert_close(
    inv_tau(transformed, pc), x, "IYDT transform round trip",
    atol = 1e-12, rtol = 1e-12
  )
  analytical <- j_tau(x, pc)
  numerical <- numDeriv::jacobian(tau, x, pc = pc)
  assert_close(
    analytical, numerical, "IYDT transform Jacobian",
    atol = 1e-8, rtol = 1e-8
  )
  assert_close(
    log(abs(det(analytical))), log_absdet_j_tau(x, pc),
    "IYDT transform log determinant", atol = 1e-12, rtol = 1e-12
  )
  assert_close(
    dlog_absdet_j_tau(x, pc),
    numDeriv::grad(log_absdet_j_tau, x, pc = pc),
    "IYDT log-Jacobian gradient", atol = 1e-8, rtol = 1e-8
  )
}

points <- c(-2, -0.5, 0, 0.5, 2)
assert_close(
  dnotExp(points), numDeriv::grad(function(value) sum(notExp(value)), points),
  "IYDT notExp derivative", atol = 1e-7, rtol = 1e-7
)

epsilon <- 1e-7
for (boundary in c(-1, 1)) {
  offsets <- boundary + c(-epsilon, 0, epsilon)
  values <- notExp(offsets)
  derivatives <- dnotExp(offsets)
  if (max(abs(values - values[2])) > 1e-6) {
    stop(sprintf("notExp is discontinuous near %s", boundary), call. = FALSE)
  }
  if (max(abs(derivatives - derivatives[2])) > 1e-6) {
    stop(sprintf("dnotExp is discontinuous near %s", boundary), call. = FALSE)
  }
}

message("IYDT transform and boundary checks passed")
