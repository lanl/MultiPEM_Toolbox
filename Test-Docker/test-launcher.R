#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("test-launcher.R must be run with Rscript")
script_path <- sub("^--file=", "", script_argument[[1L]])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
launcher <- file.path(script_dir, "testmpem.R")
rscript <- file.path(R.home("bin"), "Rscript")

run_launcher <- function(arguments, path = launcher, environment = character()) {
  output <- suppressWarnings(system2(
    rscript, c(shQuote(path), vapply(arguments, shQuote, character(1L))),
    stdout = TRUE, stderr = TRUE, env = environment
  ))
  status <- attr(output, "status", exact = TRUE)
  if (is.null(status)) status <- 0L
  list(status = as.integer(status), output = output)
}

expect_success <- function(arguments, pattern = NULL) {
  result <- run_launcher(arguments)
  if (result$status != 0L ||
      (!is.null(pattern) && !any(grepl(pattern, result$output, fixed = TRUE)))) {
    stop("Launcher command did not succeed as expected: ",
         paste(arguments, collapse = " "), call. = FALSE)
  }
  invisible(result)
}

expect_failure <- function(arguments, pattern, path = launcher,
                           environment = character()) {
  result <- run_launcher(arguments, path = path, environment = environment)
  if (result$status == 0L || !any(grepl(pattern, result$output, fixed = TRUE))) {
    stop("Launcher command did not fail as expected: ",
         paste(arguments, collapse = " "), "\n", paste(result$output, collapse = "\n"),
         call. = FALSE)
  }
  invisible(result)
}

listed <- expect_success("list")
if (!"global" %in% listed$output) stop("Launcher list omitted global")
known_suite <- listed$output[listed$output != "global"][[1L]]
expect_success("--help", "MultiPEM Docker test runner")
expect_failure(c("run", "../IYDT"), "invalid test path")
expect_failure(c("run", "IYDT/No-Such-Suite"), "unknown test")
expect_failure(c("run", known_suite, "--results"), "--results needs a value")
expect_failure(c("run", known_suite, "--cpus", "0"),
               "--cpus must be a positive number")
expect_failure("status", "status requires one job name")
expect_failure("not-a-command", "unknown command")
expect_failure("--help", "TEST_MPEM_VOLUME_LABEL must be empty, z, or Z",
               environment = "TEST_MPEM_VOLUME_LABEL=invalid")

fixture_root <- tempfile("testmpem-launcher-contract.")
dir.create(file.path(fixture_root, "Test-Docker"), recursive = TRUE)
dir.create(file.path(fixture_root, "Test", "IYDT", "Seismic"), recursive = TRUE)
invisible(file.copy(
  launcher, file.path(fixture_root, "Test-Docker", "testmpem.R")
))
writeLines("invisible(TRUE)",
           file.path(fixture_root, "Test", "IYDT", "Seismic", "tests.r"))
writeLines("tests\tcode\nIYDT\tIYDT-gsrp",
           file.path(fixture_root, "Test-Docker", "test-groups.tsv"))
expect_failure(
  c("run", "IYDT/Seismic"),
  paste(
    "test-groups.tsv must contain tests, code, fixtures, secondary_code,",
    "smoke_runfiles, smoke_case, and smoke_data columns"
  ),
  path = file.path(fixture_root, "Test-Docker", "testmpem.R")
)
unlink(fixture_root, recursive = TRUE, force = TRUE)

message("Docker test-launcher contract and negative-path checks passed")
