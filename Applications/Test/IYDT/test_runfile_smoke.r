run_runfile_smoke <- function(case_dir, label)
{
  required <- c("runMPEM.r", "runMPEM_0.r")
  missing <- required[!file.exists(file.path(case_dir, required))]
  if (length(missing)) {
    stop(sprintf("%s is missing %s", label, paste(missing, collapse = ", ")),
         call. = FALSE)
  }

  replace_setting <- function(path, name, value)
  {
    lines <- readLines(path, warn = FALSE)
    pattern <- paste0("^[[:space:]]*", name, "[[:space:]]*=")
    hits <- grep(pattern, lines)
    if (length(hits) != 1L) {
      stop(sprintf("%s must assign %s exactly once", basename(path), name),
           call. = FALSE)
    }
    lines[hits] <- paste(name, "=", value)
    writeLines(lines, path, useBytes = TRUE)
  }

  for (script in file.path(case_dir, required)) {
    replace_setting(script, "nstart", "1")
    replace_setting(script, "ncores_mle", "1")
    replace_setting(script, "mle_grad_ck", "FALSE")
    replace_setting(script, "iBayes", "FALSE")
  }
  replace_setting(file.path(case_dir, "runMPEM_0.r"), "nimpute", "1")

  old <- setwd(case_dir)
  on.exit(setwd(old), add = TRUE)
  unlink(c(".RData", "opt.RData", "opt_nev.RData",
           "calibration-pristine.RData", "calibration.out", "event.out"),
         force = TRUE)

  run_batch <- function(script, transcript, restore)
  {
    arguments <- c("CMD", "BATCH", "--no-save")
    if (!restore) arguments <- c(arguments, "--no-restore")
    arguments <- c(arguments, script, transcript)
    status <- system2(file.path(R.home("bin"), "R"), arguments)
    if (!identical(status, 0L)) {
      stop(sprintf("%s failed with status %d; inspect %s", script, status,
                   file.path(case_dir, transcript)), call. = FALSE)
    }
    output <- readLines(transcript, warn = FALSE)
    if (any(grepl("Execution halted|^Error( in|:)", output))) {
      stop(sprintf("%s reported an R error", script), call. = FALSE)
    }
  }

  run_batch("runMPEM.r", "calibration.out", restore = FALSE)
  calibration_outputs <- c(".RData", "opt.RData")
  if (any(!file.exists(calibration_outputs))) {
    stop(sprintf("%s calibration did not create %s", label,
                 paste(calibration_outputs[!file.exists(calibration_outputs)],
                       collapse = ", ")), call. = FALSE)
  }
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve the calibration checkpoint", call. = FALSE)
  }
  pristine_hash <- unname(tools::md5sum("calibration-pristine.RData"))

  run_batch("runMPEM_0.r", "event.out", restore = TRUE)
  if (!file.exists("opt_nev.RData")) {
    stop(sprintf("%s event analysis did not create opt_nev.RData", label),
         call. = FALSE)
  }
  if (!identical(unname(tools::md5sum("calibration-pristine.RData")),
                 pristine_hash)) {
    stop(sprintf("%s changed the pristine calibration checkpoint", label),
         call. = FALSE)
  }
  if (identical(unname(tools::md5sum(".RData")), pristine_hash)) {
    stop(sprintf("%s event analysis did not update its working checkpoint", label),
         call. = FALSE)
  }
  message(sprintf("%s real calibration/event runfiles passed", label))
  invisible(TRUE)
}
