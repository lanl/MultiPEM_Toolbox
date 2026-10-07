#!/usr/bin/env Rscript

local({
  arguments <- commandArgs(trailingOnly = TRUE)
  if (length(arguments) != 2L) {
    stop("Usage: run-test.R TEST_SCRIPT OUTPUT.RData", call. = FALSE)
  }

  test_script <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
  output_file <- arguments[[2L]]

  # Run canonical test code in the global environment, just as an ordinary
  # Rscript invocation does, then preserve only objects created there.
  sys.source(test_script, envir = .GlobalEnv, keep.source = FALSE)
  save.image(file = output_file)
})
