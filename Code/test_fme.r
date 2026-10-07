########################################################################
#                                                                      #
# © 2023. Triad National Security, LLC. All rights reserved.           #
# This program was produced under U.S. Government contract             #
# 89233218CNA000001 for Los Alamos National Laboratory (LANL), which   #
# is operated by Triad National Security, LLC for the U.S. Department  #
# of Energy/National Nuclear Security Administration. All rights in    #
# the program are reserved by Triad National Security, LLC, and the    #
# U.S. Department of Energy/National Nuclear Security Administration.  #
# The Government is granted for itself and others acting on its behalf #
# a nonexclusive, paid-up, irrevocable worldwide license in this       #
# material to reproduce, prepare derivative works, distribute copies   #
# to the public, perform publicly and display publicly, and to permit  #
# others to do so.                                                     #
#                                                                      #
########################################################################

if (!requireNamespace("FME", quietly = TRUE)) {
  stop("The FME package is required", call. = FALSE)
}

set.seed(100)
mu_0 <- -3
sig_0 <- 4
sig <- 2
y <- rnorm(1, mu_0, sig)

ll <- function(x) dnorm(y, x, sig, log = TRUE)
lp <- function(x) dnorm(x, mu_0, sig_0, log = TRUE)

burnin <- 5000L
niter <- burnin + 20000L
samp_po <- FME::modMCMC(
  function(x) -2 * ll(x), 0, prior = function(x) -2 * lp(x),
  jump = 0.1, niter = niter, burninlength = burnin,
  updatecov = 100, ntrydr = 2, verbose = FALSE
)

sample <- as.numeric(samp_po$pars)
sig2_p <- sig^2 * sig_0^2 / (sig^2 + sig_0^2)
mu_p <- sig2_p * (y / sig^2 + mu_0 / sig_0^2)
sig_p <- sqrt(sig2_p)

if (length(sample) != niter - burnin || !all(is.finite(sample))) {
  stop("FME returned an invalid posterior sample", call. = FALSE)
}
mean_error <- abs(mean(sample) - mu_p)
sd_error <- abs(sd(sample) - sig_p)
if (mean_error > 0.1) {
  stop(sprintf("FME posterior mean error %.6f exceeds 0.1", mean_error),
       call. = FALSE)
}
if (sd_error > 0.1) {
  stop(sprintf("FME posterior SD error %.6f exceeds 0.1", sd_error),
       call. = FALSE)
}
message(sprintf(
  "FME posterior assertions passed (mean error %.6f; SD error %.6f)",
  mean_error, sd_error
))
