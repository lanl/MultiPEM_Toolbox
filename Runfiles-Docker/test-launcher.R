#!/usr/bin/env Rscript

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("test-launcher.R must be run with Rscript")
script_path <- sub("^--file=", "", script_argument[[1L]])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- dirname(script_dir)
launcher <- file.path(script_dir, "mpem.R")
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
analyses <- listed$output[nzchar(listed$output)]
event_analysis <- analyses[vapply(analyses, function(analysis) {
  file.exists(file.path(repo_root, "Runfiles", analysis, "runMPEM_0.r"))
}, logical(1))][[1L]]

expect_success("--help", "MultiPEM Docker job runner")
expect_failure(c("run", "../Runfiles"), "invalid analysis path")
expect_failure(c("run", "No-Such-Application/No-Such-Analysis"),
               "unknown analysis")
expect_failure(c("run", event_analysis, "--stage", "invalid"), "invalid stage")
expect_failure(c("run", event_analysis, "--workspace", "missing.RData"),
               "--workspace is only valid with --stage event")
expect_failure(c("run", event_analysis, "--cpus", "0"),
               "--cpus must be a positive number")
expect_failure("status", "status requires one job name")
expect_failure("not-a-command", "unknown command")
expect_failure("--help", "MPEM_VOLUME_LABEL must be empty, z, or Z",
               environment = "MPEM_VOLUME_LABEL=invalid")

make_existing_job <- function(root, analysis, recorded = analysis) {
  metadata <- file.path(root, ".mpem")
  dir.create(file.path(metadata, "checkpoints"), recursive = TRUE)
  dir.create(file.path(root, "work", analysis), recursive = TRUE)
  writeLines(recorded, file.path(metadata, "analysis"))
  writeLines("file\tline\tvariable\tcontrol",
             file.path(metadata, "worker-controls.tsv"))
}

missing_root <- tempfile("mpem-missing-checkpoint.")
make_existing_job(missing_root, event_analysis)
expect_failure(
  c("run", event_analysis, "--stage", "event", "--results", missing_root),
  "event stage requires a pristine calibration checkpoint"
)

checkpoint_dir <- file.path(missing_root, ".mpem", "checkpoints")
writeBin(charToRaw("not an R workspace"),
         file.path(checkpoint_dir, "calibration.RData"))
writeLines(
  paste(strrep("0", 64L), " calibration.RData"),
  file.path(checkpoint_dir, "calibration.RData.sha256")
)
expect_failure(
  c("run", event_analysis, "--stage", "event", "--results", missing_root),
  "calibration checkpoint checksum does not match"
)

incompatible_root <- tempfile("mpem-incompatible-results.")
make_existing_job(incompatible_root, event_analysis, recorded = "another/analysis")
expect_failure(
  c("run", event_analysis, "--stage", "event", "--results", incompatible_root),
  "results belong to another/analysis"
)

fixture_root <- tempfile("mpem-launcher-contract.")
for (directory in c("Runfiles-Docker", "Runfiles", "Applications", "Code", "Test")) {
  dir.create(file.path(fixture_root, directory), recursive = TRUE)
}
invisible(file.copy(
  launcher, file.path(fixture_root, "Runfiles-Docker", "mpem.R")
))
writeLines("runfiles\tcode\nIYDT-gsrp\tIYDT-gsrp",
           file.path(fixture_root, "Runfiles-Docker", "applications.tsv"))
expect_failure(
  "doctor", "applications.tsv must contain runfiles, code, and data columns",
  path = file.path(fixture_root, "Runfiles-Docker", "mpem.R")
)

unlink(c(missing_root, incompatible_root, fixture_root),
       recursive = TRUE, force = TRUE)
message("Docker runfile-launcher contract and negative-path checks passed")
