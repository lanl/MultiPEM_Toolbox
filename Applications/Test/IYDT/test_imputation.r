copy_imputation_case <- function(source, destination)
{
  if (dir.exists(destination)) unlink(destination, recursive = TRUE, force = TRUE)
  if (!dir.create(destination, recursive = TRUE)) {
    stop(sprintf("could not create imputation case %s", destination),
         call. = FALSE)
  }
  entries <- list.files(source, all.files = TRUE, no.. = TRUE,
                        full.names = TRUE)
  if (length(entries) &&
      !all(file.copy(entries, destination, recursive = TRUE,
                     copy.mode = TRUE, copy.date = TRUE))) {
    stop(sprintf("could not copy imputation case %s", source), call. = FALSE)
  }
  invisible(destination)
}

run_multiple_imputation_smoke <- function(case_dir, label,
                                          smc_lower, smc_upper)
{
  patch_vector_setting <- function(path, name, value) {
    lines <- readLines(path, warn = FALSE)
    hits <- grep(paste0("^[[:space:]]*", name, "[[:space:]]+="), lines)
    if (length(hits) != 1L) {
      stop(sprintf("%s must assign %s exactly once", basename(path), name),
           call. = FALSE)
    }
    lines[hits] <- paste(name, "=", value)
    writeLines(lines, path, useBytes = TRUE)
  }
  work_dir <- file.path(dirname(case_dir),
                        paste0(basename(case_dir), "-imputation"))
  copy_imputation_case(case_dir, work_dir)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)
  prepare_mle_scripts(work_dir, bayes = TRUE)
  scripts <- file.path(work_dir, c("runMPEM.r", "runMPEM_0.r"))
  scripts <- normalizePath(scripts, winslash = "/", mustWork = TRUE)
  for (script in scripts) {
    patch_runfile_setting(script, "parallel_plan", '"multisession"')
  }
  event_script <- scripts[[2L]]
  patch_runfile_setting(event_script, "nimpute", "2")
  patch_runfile_setting(event_script, "ncores_mc", "2")
  patch_runfile_setting(event_script, "ncores_smc", "1")
  patch_vector_setting(event_script, "lb_smc", deparse(smc_lower))
  patch_vector_setting(event_script, "ub_smc", deparse(smc_upper))

  old <- setwd(work_dir)
  on.exit(setwd(old), add = TRUE)
  unlink(c(".RData", "opt.RData", "opt_nev.RData",
           "calibration-pristine.RData", "imputation-calibration.out",
           "imputation-ram.out", "imputation-smc.out"), force = TRUE)

  run_batch_checked("runMPEM.r", "imputation-calibration.out", restore = FALSE)
  calibration <- load_workspace(".RData")
  calibration_sample <- assert_finite_sample(
    calibration$mpi, paste(label, "calibration posterior"), minimum_rows = 2L
  )
  if (nrow(calibration_sample) < 2L) {
    stop("multiple-imputation test needs two calibration draws", call. = FALSE)
  }
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve the multiple-imputation checkpoint", call. = FALSE)
  }
  pristine_hash <- unname(tools::md5sum("calibration-pristine.RData"))

  run_backend <- function(backend) {
    if (!file.copy("calibration-pristine.RData", ".RData", overwrite = TRUE)) {
      stop("could not restore the multiple-imputation checkpoint", call. = FALSE)
    }
    unlink("opt_nev.RData", force = TRUE)
    patch_runfile_setting(event_script, "iMCMC", sprintf('"%s"', backend))
    if (identical(backend, "RAM")) {
      patch_runfile_setting(event_script, "nburn", "10")
      patch_runfile_setting(event_script, "nmcmc", "20")
      patch_runfile_setting(event_script, "nthin", "2")
    } else {
      patch_runfile_setting(event_script, "nburn", "1")
      patch_runfile_setting(event_script, "nmcmc", "40")
      patch_runfile_setting(event_script, "nthin", "4")
    }
    transcript <- sprintf("imputation-%s.out", tolower(backend))
    run_batch_checked("runMPEM_0.r", transcript, restore = TRUE)
    event <- load_workspace(".RData")
    posterior <- assert_finite_sample(
      event$tmpi_0, paste(label, backend, "two-imputation posterior"),
      minimum_rows = 20L
    )
    if (nrow(posterior) != 20L || ncol(posterior) != event$ntheta0) {
      stop(sprintf("%s %s posterior was not pooled as 2 x 10 draws",
                   label, backend), call. = FALSE)
    }
    output <- readLines(transcript, warn = FALSE)
    for (index in 1:2) {
      marker <- paste0("Imputation ", index, ":")
      if (!any(grepl(marker, output, fixed = TRUE))) {
        stop(sprintf("%s %s transcript is missing %s",
                     label, backend, marker), call. = FALSE)
      }
    }
    invisible(posterior)
  }

  run_backend("RAM")
  run_backend("SMC")
  if (!identical(unname(tools::md5sum("calibration-pristine.RData")),
                 pristine_hash)) {
    stop(sprintf("%s changed the pristine multiple-imputation checkpoint",
                 label), call. = FALSE)
  }
  message(sprintf("%s two-imputation RAM/SMC workflow passed", label))
  invisible(TRUE)
}
