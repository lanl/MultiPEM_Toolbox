# Loaded before each test run. It limits only concurrent worker-pool arguments;
# the logical ncores* values remain unchanged so calculation structure is kept.
local({
detected <- suppressWarnings({
  if (requireNamespace("parallelly", quietly = TRUE)) {
    as.numeric(parallelly::availableCores())
  } else {
    as.numeric(parallel::detectCores(logical = TRUE))
  }
})
detected <- if (length(detected)) detected[[1L]] else NA_real_
if (!is.finite(detected) || detected < 1) detected <- 1

explicit <- suppressWarnings(as.numeric(Sys.getenv("TEST_MPEM_MAX_WORKERS", "")))
explicit <- if (length(explicit)) explicit[[1L]] else NA_real_
if (!is.finite(explicit) || explicit < 1) explicit <- Inf
maximum <- max(1L, floor(min(detected, explicit)))
log_file <- Sys.getenv("TEST_MPEM_CPU_LOG", "")
reported <- new.env(parent = emptyenv())

append_cpu_log <- function(text) {
  if (nzchar(log_file)) {
    cat(
      format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
      "\t", text, "\n", sep = "", file = log_file, append = TRUE
    )
  }
}

display_name <- function(label) {
  switch(label,
    ncor = "ncores_mle",
    ncor_map = "ncores_map",
    ncor_mc = "ncores_mc",
    ncor_smc = "ncores_smc",
    label
  )
}

limit_workers <- function(value, label, mode) {
  if (is.null(value) || length(value) != 1L || !is.numeric(value) ||
      !is.finite(value) || value <= maximum) return(value)

  requested <- as.numeric(value)
  assigned <- as.numeric(maximum)
  name <- display_name(label)
  key <- paste(name, requested, assigned, mode, sep = "|")
  if (!exists(key, envir = reported, inherits = FALSE)) {
    text <- sprintf(
      paste0(
        "available CPU count is %d; effective %s worker assignment reduced ",
        "from %g to %d; logical %s remains %g to preserve calculation structure"
      ),
      maximum, name, requested, maximum, name, requested
    )
    warning_text <- paste0("MultiPEM test CPU adjustment: ", text)
    warning(warning_text, call. = FALSE, immediate. = TRUE)
    append_cpu_log(warning_text)
    assign(key, TRUE, envir = reported)
  }
  as.integer(maximum)
}

options(
  testmpem.available.workers = maximum,
  testmpem.limit.workers = limit_workers
)

capacity_message <- sprintf(
  "MultiPEM test CPU capacity: %d concurrent worker(s)", maximum
)
message(capacity_message)
append_cpu_log(capacity_message)
})
