assert_rejected <- function(label, expression)
{
  rejected <- FALSE
  value <- tryCatch(
    withCallingHandlers(
      expression,
      warning = function(w) invokeRestart("muffleWarning")
    ),
    error = function(e) {
      rejected <<- TRUE
      NULL
    }
  )
  if (!rejected) {
    rejected <- !is.numeric(value) || !length(value) || !all(is.finite(value))
  }
  if (!rejected) {
    stop(sprintf("%s unexpectedly produced a finite model result", label),
         call. = FALSE)
  }
  message(sprintf("%s rejected", label))
  invisible(TRUE)
}

soft_exp <- function(x)
{
  value <- exp(pmin(x, 1))
  high <- x > 1
  value[high] <- exp(1) * (x[high]^2 + 1) / 2
  value
}

base_x <- matrix(
  c(log(120), log(8), 0.1, 50), nrow = 1,
  dimnames = list(NULL, c("lRange", "W", "C2N", "HOB"))
)
base_params <- list(
  pbeta = 5, cal = FALSE, cal_par_names = character(), ncalp = 0,
  theta_names = character(), iresp = TRUE, yield_scaling = 1 / 3,
  X = base_x, notExp = soft_exp
)
beta <- c(1.1, -0.65, 0.4, 0.08, -0.3)

missing_range <- base_params
missing_range$X <- missing_range$X[, setdiff(colnames(missing_range$X), "lRange"),
                                   drop = FALSE]
assert_rejected("missing required lRange covariate", f_s(beta, missing_range))
assert_rejected("short empirical parameter vector", f_s(beta[1], base_params))

nonfinite_yield <- base_params
nonfinite_yield$X[, "W"] <- Inf
assert_rejected("non-finite yield covariate", f_s(beta, nonfinite_yield))

invalid_scaling <- base_params
invalid_scaling$yield_scaling <- NA_real_
assert_rejected("non-finite yield scaling", f_s(beta, invalid_scaling))
