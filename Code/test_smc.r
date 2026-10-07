########################################################################
# Application-independent regression tests for the SMC implementation. #
########################################################################

required <- c("ramcmc", "doFuture", "future", "iterators")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("SMC tests require: ", paste(missing, collapse = ", "), call. = FALSE)
}

doFuture::registerDoFuture()
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
source(file.path(tdir, "pSMC-ram-mv.r"), local = TRUE)

expect_smc_error <- function(expression, pattern)
{
  error <- tryCatch({
    force(expression)
    NULL
  }, error = identity)
  if (is.null(error) || !grepl(pattern, conditionMessage(error))) {
    stop(sprintf("expected an error containing '%s'", pattern), call. = FALSE)
  }
  invisible(error)
}

particle_count <- 360L

assert_smc_fit <- function(fit, label, target_mean, target_cov,
                           expected_rows = particle_count,
                           check_moments = TRUE)
{
  sample <- as.matrix(fit$sample)
  if (!identical(dim(sample), c(as.integer(expected_rows), 2L)) ||
      any(!is.finite(sample))) {
    stop(sprintf("%s returned an invalid %d by 2 particle matrix",
                 label, expected_rows),
         call. = FALSE)
  }
  if (!length(fit$acr) || any(!is.finite(fit$acr)) ||
      any(fit$acr < 0 | fit$acr > 1) || !any(fit$acr > 0)) {
    stop(sprintf("%s returned invalid acceptance diagnostics", label),
         call. = FALSE)
  }
  if (check_moments && max(abs(colMeans(sample) - target_mean)) > 0.22) {
    stop(sprintf("%s posterior mean is outside tolerance", label), call. = FALSE)
  }
  if (check_moments && max(abs(cov(sample) - target_cov)) > 0.32) {
    stop(sprintf("%s posterior covariance is outside tolerance", label),
         call. = FALSE)
  }
  invisible(sample)
}

target_cov <- matrix(c(0.70, 0.24, 0.24, 0.50), 2L, 2L)
target_mean <- c(0.60, -0.80)
model <- list(
  mean = target_mean,
  precision = solve(target_cov),
  lprior_0 = function(x, pc) {
    if (length(x) == 2L && all(abs(x) <= 6)) 0 else -Inf
  },
  ll_0 = function(x, pc) {
    delta <- x - pc$mean
    -0.5 * drop(crossprod(delta, pc$precision %*% delta))
  }
)
box <- cbind(c(-5, -5), c(5, 5))

message("SMC Gaussian target, diagnostics, and reproducibility")
future::plan(future::sequential)
set.seed(48291)
sequential <- SMC(N = particle_count, M = 8, range = box, pc = model,
                  range_control = "off")
sequential_sample <- assert_smc_fit(
  sequential, "sequential SMC", target_mean, target_cov
)
set.seed(48291)
repeat_fit <- SMC(N = particle_count, M = 8, range = box, pc = model,
                  range_control = "off")
if (!identical(sequential, repeat_fit)) {
  stop("fixed-seed sequential SMC is not exactly reproducible", call. = FALSE)
}

message("SMC sequential/two-worker statistical equivalence")
future::plan(future::multisession, workers = 2L)
set.seed(48291)
parallel <- SMC(N = particle_count, M = 8, range = box, pc = model,
                range_control = "off")
parallel_sample <- assert_smc_fit(
  parallel, "two-worker SMC", target_mean, target_cov
)
if (max(abs(colMeans(sequential_sample) - colMeans(parallel_sample))) > 0.25 ||
    max(abs(cov(sequential_sample) - cov(parallel_sample))) > 0.38) {
  stop("sequential and two-worker SMC results are not statistically equivalent",
       call. = FALSE)
}
future::plan(future::sequential)

message("SMC adaptive/warn range controls")
narrow_box <- cbind(target_mean - 0.05, target_mean + 0.05)
set.seed(731)
adapted <- SMC(
  N = 240, M = 6, range = narrow_box, pc = model,
  range_control = "adapt", range_max_iter = 8, range_tol = 0.08
)
assert_smc_fit(adapted, "range-adapted SMC", target_mean, target_cov,
               expected_rows = 240L)

wide_box <- cbind(c(-20, -20), c(20, 20))
set.seed(913)
warning_messages <- capture.output(
  warned <- SMC(N = 80, M = 2, range = wide_box, pc = model,
                range_control = "warn"),
  type = "message"
)
assert_smc_fit(warned, "range-warning SMC", target_mean, target_cov,
               expected_rows = 80L, check_moments = FALSE)
if (!any(grepl("range.*WIDE|initial ESS", warning_messages,
               ignore.case = TRUE))) {
  stop("SMC warn mode did not report the deliberately wide range",
       call. = FALSE)
}

message("SMC particle-pool sampling")
pool <- rbind(seq(-4, 4, length.out = 401),
              seq(4, -4, length.out = 401))
set.seed(331)
without_replacement <- unrestricted(40, range = NULL, psamp = pool)
if (!identical(dim(without_replacement), c(2L, 40L)) ||
    anyDuplicated(as.data.frame(t(without_replacement)))) {
  stop("SMC particle pools should sample without replacement when N <= L",
       call. = FALSE)
}
set.seed(332)
with_replacement <- unrestricted(450, range = NULL, psamp = pool)
if (!identical(dim(with_replacement), c(2L, 450L)) ||
    !anyDuplicated(as.data.frame(t(with_replacement)))) {
  stop("SMC particle pools should sample with replacement when N > L",
       call. = FALSE)
}
set.seed(333)
prior_pool <- rbind(runif(500, -5, 5), runif(500, -5, 5))
pool_fit <- SMC(N = 120, M = 4, psamp = prior_pool, pc = model,
                range_control = "off")
assert_smc_fit(pool_fit, "particle-pool SMC", target_mean, target_cov,
               expected_rows = 120L, check_moments = FALSE)

message("SMC mixed finite/non-finite likelihood recovery")
partial_model <- model
partial_model$ll_0 <- function(x, pc) {
  if (x[[1L]] < -1) return(NaN)
  delta <- x - pc$mean
  -0.5 * drop(crossprod(delta, pc$precision %*% delta))
}
set.seed(811)
partial_messages <- capture.output(
  partial_fit <- SMC(N = 120, M = 4, range = box, pc = partial_model,
                     range_control = "off"),
  type = "message"
)
assert_smc_fit(partial_fit, "partially non-finite SMC", target_mean, target_cov,
               expected_rows = 120L, check_moments = FALSE)
if (!any(grepl("nonfinite", partial_messages, fixed = TRUE))) {
  stop("SMC did not report sanitizing the deliberately non-finite particles",
       call. = FALSE)
}

message("SMC effective-sample-size calculation")
log_likelihood <- c(-4, -2, -1, -0.5, 0)
base_weights <- rep(1 / length(log_likelihood), length(log_likelihood))
increment <- 0.65
weights <- base_weights * exp(increment * log_likelihood)
weights <- weights / sum(weights)
manual_ess <- 1 / sum(weights^2)
calculated_ess <- adapt_seq(
  increment, 0, log_likelihood, base_weights,
  length(log_likelihood), r_ess = 0
)
if (!isTRUE(all.equal(calculated_ess, manual_ess, tolerance = 1e-14)) ||
    abs(sum(weights) - 1) > 1e-14 || manual_ess < 1 ||
    manual_ess > length(log_likelihood)) {
  stop("SMC ESS or normalized-weight calculation is incorrect", call. = FALSE)
}

message("SMC input validation and zero-mass failures")
expect_smc_error(validate_range(c(-1, 1)), "matrix with two columns")
expect_smc_error(validate_range(matrix(c(-Inf, 1), 1L, 2L)), "non-finite")
expect_smc_error(validate_range(matrix(c(1, 1), 1L, 2L)), "lower bound")
expect_smc_error(validate_range(box, D_psamp = 3L), "psamp")
expect_smc_error(SMC(N = 20, M = 1, pc = model), "either `range` or `psamp`")

bad_model <- model
bad_model$ll_0 <- function(x, pc) NaN
expect_smc_error(
  SMC(N = 30, M = 1, range = box, pc = bad_model, range_control = "off"),
  "zero posterior mass"
)

sanitized <- sanitize_lpden(matrix(
  c(-Inf, NA, 0, Inf, NaN, -Inf), nrow = 2L
))
if (!identical(unname(attr(sanitized, "n_bad")), c(1L, 2L)) ||
    !identical(sanitized[1, 1], -Inf) || any(is.na(sanitized))) {
  stop("SMC non-finite density sanitization contract failed", call. = FALSE)
}

message("SMC regression tests passed")
