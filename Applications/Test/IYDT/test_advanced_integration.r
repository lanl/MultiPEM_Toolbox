copy_case_tree <- function(source, destination)
{
  if (!dir.exists(source)) {
    stop(sprintf("case directory does not exist: %s", source), call. = FALSE)
  }
  if (dir.exists(destination)) unlink(destination, recursive = TRUE, force = TRUE)
  if (!dir.create(destination, recursive = TRUE)) {
    stop(sprintf("could not create case directory: %s", destination),
         call. = FALSE)
  }
  entries <- list.files(source, all.files = TRUE, no.. = TRUE,
                        full.names = TRUE)
  if (length(entries) &&
      !all(file.copy(entries, destination, recursive = TRUE,
                     copy.mode = TRUE, copy.date = TRUE))) {
    stop(sprintf("could not copy case directory: %s", source), call. = FALSE)
  }
  invisible(destination)
}

patch_advanced_runfiles <- function(case_dir, bayes = FALSE, backend = NULL,
                                    workers = 1L, starts = 1L)
{
  prepare_mle_scripts(case_dir, bayes = bayes)
  scripts <- file.path(case_dir, c("runMPEM.r", "runMPEM_0.r"))
  for (script in scripts) {
    patch_runfile_setting(script, "nstart", as.character(starts))
    patch_runfile_setting(script, "ncores_mle", as.character(workers))
    patch_runfile_setting(script, "parallel_plan", '"multisession"')
    if (bayes) {
      nuts <- identical(backend, "NUTS")
      patch_runfile_setting(script, "nburn", if (nuts) "3" else "10")
      patch_runfile_setting(script, "nmcmc", if (nuts) "6" else "20")
      patch_runfile_setting(script, "nthin", if (nuts) "1" else "2")
      patch_runfile_setting(script, "ncores_map", as.character(workers))
      patch_runfile_setting(script, "ncores_mc", as.character(workers))
      patch_runfile_setting(script, "iMCMC", sprintf('"%s"', backend))
    }
  }
  invisible(scripts)
}

assert_numeric_equal <- function(left, right, label, tolerance = 1e-8)
{
  left <- as.numeric(left)
  right <- as.numeric(right)
  if (!length(left) || length(left) != length(right) ||
      any(!is.finite(left)) || any(!is.finite(right))) {
    stop(sprintf("%s does not contain comparable finite values", label),
         call. = FALSE)
  }
  error <- max(abs(left - right))
  scale <- max(1, abs(left), abs(right))
  if (error > tolerance * scale) {
    stop(sprintf("%s differs by %.8g (tolerance %.8g)",
                 label, error, tolerance * scale), call. = FALSE)
  }
  message(sprintf("%s agreed (maximum absolute difference %.8g)",
                  label, error))
  invisible(error)
}

run_application_backend <- function(source_case, backend, label)
{
  backend <- match.arg(backend, c("RAM", "NUTS"))
  case_dir <- file.path(dirname(source_case),
                        paste0(basename(source_case), "-", tolower(backend)))
  copy_case_tree(source_case, case_dir)
  on.exit(unlink(case_dir, recursive = TRUE, force = TRUE), add = TRUE)
  patch_advanced_runfiles(case_dir, bayes = TRUE, backend = backend)

  old <- setwd(case_dir)
  on.exit(setwd(old), add = TRUE)
  run_batch_checked("runMPEM.r", paste0(tolower(backend), "-calibration.out"),
                    restore = FALSE)
  calibration <- load_workspace(".RData")
  assert_finite_sample(calibration$mpi,
                       paste(label, backend, "calibration posterior"),
                       minimum_rows = if (backend == "NUTS") 6L else 10L)
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve backend calibration checkpoint", call. = FALSE)
  }
  pristine_hash <- unname(tools::md5sum("calibration-pristine.RData"))

  run_batch_checked("runMPEM_0.r", paste0(tolower(backend), "-event.out"),
                    restore = TRUE)
  event <- load_workspace(".RData")
  assert_finite_sample(event$tmpi_0,
                       paste(label, backend, "event posterior"),
                       minimum_rows = if (backend == "NUTS") 6L else 10L)
  if (!identical(unname(tools::md5sum("calibration-pristine.RData")),
                 pristine_hash)) {
    stop(sprintf("%s changed the pristine backend checkpoint", label),
         call. = FALSE)
  }
  transcript_paths <- c(paste0(tolower(backend), "-calibration.out"),
                        paste0(tolower(backend), "-event.out"))
  transcripts <- unlist(lapply(
    transcript_paths, readLines, warn = FALSE
  ), use.names = FALSE)
  required <- if (backend == "NUTS") "DIVERGENCES:" else "ACCEPTANCE RATES:"
  if (!any(grepl(required, transcripts, fixed = TRUE))) {
    stop(sprintf("%s transcripts do not demonstrate the %s backend",
                 label, backend), call. = FALSE)
  }
  message(sprintf("%s real-application %s calibration/event workflow passed",
                  label, backend))
  invisible(TRUE)
}

run_parallel_equivalence <- function(source_case, label)
{
  parent <- dirname(source_case)
  sequential_dir <- file.path(parent, paste0(basename(source_case), "-worker-1"))
  parallel_dir <- file.path(parent, paste0(basename(source_case), "-workers-2"))
  copy_case_tree(source_case, sequential_dir)
  copy_case_tree(source_case, parallel_dir)
  on.exit(unlink(c(sequential_dir, parallel_dir), recursive = TRUE, force = TRUE),
          add = TRUE)

  run_case <- function(case_dir, workers) {
    patch_advanced_runfiles(case_dir, bayes = FALSE, workers = workers,
                            starts = 2L)
    old <- setwd(case_dir)
    on.exit(setwd(old), add = TRUE)
    run_batch_checked("runMPEM.r", "parallel-calibration.out", restore = FALSE)
    calibration <- load_workspace(".RData")
    if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
      stop("could not preserve parallel calibration checkpoint", call. = FALSE)
    }
    run_batch_checked("runMPEM_0.r", "parallel-event.out", restore = TRUE)
    event <- load_workspace(".RData")
    list(
      calibration_mle = calibration$mle_cal,
      calibration_loglik = calibration$ll_cal(calibration$mle_cal, calibration),
      event_mle = event$mle,
      event_loglik = event$ll_0(event$mle, event)
    )
  }

  sequential <- run_case(sequential_dir, 1L)
  parallel <- run_case(parallel_dir, 2L)
  for (name in names(sequential)) {
    assert_numeric_equal(sequential[[name]], parallel[[name]],
                         paste(label, name), tolerance = 2e-8)
  }
  message(sprintf("%s one-worker/two-worker MLE workflow passed", label))
  invisible(TRUE)
}

run_multphen_staged_inputs <- function(app_root, label)
{
  optical <- file.path(app_root, "Optical", "I-SUGAR-hob-0")
  crater <- file.path(app_root, "Crater", "I-SUGAR-0")
  combined <- file.path(app_root, "2-Phen-oc", "I-SUGAR-hob-0")
  for (case_dir in c(optical, crater, combined)) {
    if (!dir.exists(case_dir)) {
      stop(sprintf("multi-phenomenology fixture is missing %s", case_dir),
           call. = FALSE)
    }
  }

  run_single <- function(case_dir, transcript) {
    prepare_mle_scripts(case_dir, bayes = FALSE)
    patch_runfile_setting(file.path(case_dir, "runMPEM.r"), "parallel_plan",
                          '"sequential"')
    old <- setwd(case_dir)
    on.exit(setwd(old), add = TRUE)
    run_batch_checked("runMPEM.r", transcript, restore = FALSE)
    workspace <- load_workspace(".RData")
    if (!file.exists("opt.RData") || !all(is.finite(workspace$mle_cal))) {
      stop(sprintf("single-phenomenology stage failed in %s", case_dir),
           call. = FALSE)
    }
    normalizePath("opt.RData", winslash = "/", mustWork = TRUE)
  }

  optical_opt <- run_single(optical, "multphen-optical.out")
  crater_opt <- run_single(crater, "multphen-crater.out")
  opt_dir <- file.path(app_root, "2-Phen-oc", "Opt")
  dir.create(opt_dir, recursive = TRUE, showWarnings = FALSE)
  staged <- file.path(opt_dir, c("opt_1_0.RData", "opt_2_0.RData"))
  if (!all(file.copy(c(optical_opt, crater_opt), staged, overwrite = TRUE))) {
    stop("could not stage single-phenomenology opt.RData inputs", call. = FALSE)
  }
  staged <- normalizePath(staged, winslash = "/", mustWork = TRUE)
  staged_hash <- unname(tools::md5sum(staged))

  prepare_mle_scripts(combined, bayes = FALSE)
  for (script in file.path(combined, c("runMPEM.r", "runMPEM_0.r"))) {
    patch_runfile_setting(script, "parallel_plan", '"sequential"')
  }
  old <- setwd(combined)
  on.exit(setwd(old), add = TRUE)
  run_batch_checked("runMPEM.r", "multphen-calibration.out", restore = FALSE)
  calibration <- load_workspace(".RData")
  if (!identical(calibration$H, 2L) || length(calibration$h) != 2L ||
      !all(is.finite(calibration$mle_cal)) || !file.exists("opt.RData")) {
    stop("combined calibration did not integrate both phenomenologies",
         call. = FALSE)
  }
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve multi-phenomenology checkpoint", call. = FALSE)
  }
  run_batch_checked("runMPEM_0.r", "multphen-event.out", restore = TRUE)
  event <- load_workspace(".RData")
  if (!isTRUE(event$nev) || !all(is.finite(event$mle)) ||
      !file.exists("opt_nev.RData")) {
    stop("combined event analysis did not produce finite results",
         call. = FALSE)
  }
  if (!identical(unname(tools::md5sum(staged)), staged_hash)) {
    stop("combined workflow modified a staged opt.RData input", call. = FALSE)
  }
  message(sprintf("%s staged two-phenomenology calibration/event workflow passed",
                  label))
  invisible(TRUE)
}

direct_common_parameter_fits <- function(parameters, pc)
{
  if (isTRUE(pc$eiv) || isTRUE(pc$ptbeta > 0)) {
    stop("direct prediction check requires a non-EIV common-parameter case",
         call. = FALSE)
  }
  parameters <- as.numeric(parameters)
  theta0 <- numeric()
  if (isTRUE(pc$nev)) {
    theta0 <- parameters[seq_len(pc$ntheta0)]
    if (isTRUE(pc$itransform)) theta0 <- pc$tau(theta0, pc = pc)
    if (exists("itheta0_bounds", where = pc, inherits = FALSE)) {
      theta0 <- pc$transform(theta0, pc = pc)
    }
    parameters <- parameters[-seq_len(pc$ntheta0)]
  }
  calp <- numeric()
  if (pc$ncalp > 0) {
    calp <- parameters[seq_len(pc$ncalp)]
    parameters <- parameters[-seq_len(pc$ncalp)]
  }
  beta_all <- if (pc$pbeta > 0) parameters[seq_len(pc$pbeta)] else numeric()
  result <- vector("list", pc$H)

  for (hh in seq_len(pc$H)) {
    h <- pc$h[[hh]]
    h_names <- names(h)
    beta_count <- sum(h$pbeta)
    beta <- if (beta_count) beta_all[seq_len(beta_count)] else numeric()
    if (beta_count) beta_all <- beta_all[-seq_len(beta_count)]
    result[[hh]] <- vector("list", h$nsource)
    for (ii in seq_len(h$nsource)) {
      pm <- new.env(hash = TRUE)
      if (exists("notExp", where = pc, inherits = FALSE)) pm$notExp <- pc$notExp
      if ("llpars" %in% h_names) {
        for (name in names(h$llpars)) pm[[name]] <- h$llpars[[name]]
      }
      event_source <- "nev" %in% h_names && isTRUE(h$nev[[ii]])
      cp <- numeric()
      if ("cal_par_names" %in% h_names && !event_source) {
        indices <- which(pc$cal_par_names %in% h$cal_par_names)
        cp <- calp[indices]
        pm$cal <- length(cp) > 0L
        if (pm$cal) {
          pm$cal_par_names <- h$cal_par_names
          pm$ncalp <- length(pm$cal_par_names)
        }
      } else {
        pm$cal <- FALSE
      }
      if ("theta_names" %in% h_names && !is.null(h$theta_names[[ii]])) {
        pm$theta_names <- h$theta_names[[ii]]
      }
      result[[hh]][[ii]] <- vector("list", h$Rh)
      for (rr in seq_len(h$Rh)) {
        if (h$n[[ii]][[rr]] <= 0) next
        pm$X <- h$X[[ii]][[rr]]
        start <- if (rr == 1L) 0L else sum(h$pbeta[seq_len(rr - 1L)])
        beta_r <- beta[start + seq_len(h$pbeta[[rr]])]
        pm$pbeta <- length(beta_r)
        if ("iResponse" %in% h_names) pm$iresp <- h$iResponse[[rr]]
        zeta <- c(beta_r, cp)
        if (event_source) {
          event_theta <- if ("itheta0" %in% h_names) {
            theta0[h$itheta0]
          } else theta0
          zeta <- c(zeta, event_theta)
        }
        result[[hh]][[ii]][[rr]] <-
          pc$ffm[[h$f[[rr]]]](zeta, pm)
      }
    }
  }
  result
}

assert_prediction_contract <- function(prediction, pc, parameters, label)
{
  if (!is.environment(prediction) || !is.list(prediction$h) ||
      length(prediction$h) != pc$H) {
    stop(sprintf("%s prediction has an invalid phenomenology structure", label),
         call. = FALSE)
  }
  numeric_leaves <- 0L
  inspect <- function(value, path) {
    if (is.null(value)) return(invisible())
    if (is.environment(value)) value <- as.list.environment(value, all.names = TRUE)
    if (is.list(value)) {
      for (name in names(value)) inspect(value[[name]], paste(path, name, sep = "/"))
      if (is.null(names(value))) {
        for (index in seq_along(value)) inspect(value[[index]],
                                                paste0(path, "/", index))
      }
      return(invisible())
    }
    if (is.numeric(value) || inherits(value, "Matrix")) {
      values <- as.numeric(value)
      if (length(values) && any(!is.finite(values))) {
        stop(sprintf("%s contains non-finite values at %s", label, path),
             call. = FALSE)
      }
      numeric_leaves <<- numeric_leaves + 1L
    }
    invisible()
  }
  inspect(prediction$h, "h")
  for (hh in seq_len(pc$H)) {
    section <- prediction$h[[hh]]
    if (!is.list(section$observed) || !is.list(section$fitted) ||
        length(section$observed) != length(section$fitted)) {
      stop(sprintf("%s observed/fitted structure is invalid", label),
           call. = FALSE)
    }
    for (ii in seq_along(section$observed)) {
      observed <- section$observed[[ii]]
      fitted <- section$fitted[[ii]]
      if (length(observed) != length(fitted)) {
        stop(sprintf("%s source %d observed/fitted responses differ", label, ii),
             call. = FALSE)
      }
      for (rr in seq_along(observed)) {
        if (length(observed[[rr]]) != length(fitted[[rr]])) {
          stop(sprintf("%s source %d response %d lengths differ", label, ii, rr),
               call. = FALSE)
        }
        if (!isTRUE(all.equal(as.numeric(observed[[rr]]),
                              as.numeric(pc$h[[hh]]$Y[[ii]][[rr]]),
                              tolerance = 0))) {
          stop(sprintf("%s source %d response %d observations changed",
                       label, ii, rr), call. = FALSE)
        }
      }
    }
    for (gg in seq_along(section$Omega)) {
      omega <- section$Omega[[gg]]
      if (!is.null(omega)) {
        omega <- as.matrix(omega)
        if (!isTRUE(all.equal(omega, t(omega), tolerance = 1e-9)) ||
            any(!is.finite(omega))) {
          stop(sprintf("%s produced an invalid covariance matrix", label),
               call. = FALSE)
        }
        inverse <- as.matrix(section$IOmega[[gg]])
        identity_error <- max(abs(omega %*% inverse - diag(nrow(omega))))
        if (!is.finite(identity_error) || identity_error > 1e-7) {
          stop(sprintf("%s covariance inverse error %.8g exceeds tolerance",
                       label, identity_error), call. = FALSE)
        }
        grouped_rows <- sum(vapply(
          pc$h[[hh]]$Source_Groups[[gg]],
          function(ii) sum(pc$h[[hh]]$n[[ii]]), numeric(1L)
        ))
        if (nrow(omega) != grouped_rows) {
          stop(sprintf("%s covariance dimension does not match grouped residuals",
                       label), call. = FALSE)
        }
      }
    }
  }
  direct <- direct_common_parameter_fits(parameters, pc)
  for (hh in seq_len(pc$H)) {
    for (ii in seq_along(direct[[hh]])) {
      for (rr in seq_along(direct[[hh]][[ii]])) {
        if (!is.null(direct[[hh]][[ii]][[rr]]) &&
            !isTRUE(all.equal(
              as.numeric(prediction$h[[hh]]$fitted[[ii]][[rr]]),
              as.numeric(direct[[hh]][[ii]][[rr]]), tolerance = 1e-12
            ))) {
          stop(sprintf(
            "%s source %d response %d differs from its direct forward model",
            label, ii, rr
          ), call. = FALSE)
        }
      }
    }
  }
  if (!numeric_leaves) stop(sprintf("%s prediction is empty", label), call. = FALSE)
  invisible(TRUE)
}

run_prediction_regression <- function(source_case, label)
{
  case_dir <- file.path(dirname(source_case), paste0(basename(source_case),
                                                    "-prediction"))
  copy_case_tree(source_case, case_dir)
  on.exit(unlink(case_dir, recursive = TRUE, force = TRUE), add = TRUE)
  patch_advanced_runfiles(case_dir, bayes = FALSE)
  old <- setwd(case_dir)
  on.exit(setwd(old), add = TRUE)
  run_batch_checked("runMPEM.r", "prediction-calibration.out", restore = FALSE)
  calibration <- load_workspace(".RData")
  source(file.path("..", "..", "..", "Code", "predict.r"), local = TRUE)
  calibration_prediction <- predict(calibration$mle_cal, calibration)
  assert_prediction_contract(calibration_prediction, calibration,
                             calibration$mle_cal,
                             paste(label, "calibration"))
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve prediction calibration checkpoint", call. = FALSE)
  }
  pristine_hash <- unname(tools::md5sum("calibration-pristine.RData"))

  run_batch_checked("runMPEM_0.r", "prediction-event.out", restore = TRUE)
  event <- load_workspace(".RData")
  event_prediction <- predict(c(event$mle, event$mle_cal), event)
  repeated_prediction <- predict(c(event$mle, event$mle_cal), event)
  assert_prediction_contract(event_prediction, event,
                             c(event$mle, event$mle_cal),
                             paste(label, "event"))
  if (!isTRUE(all.equal(as.list.environment(event_prediction, all.names = TRUE),
                        as.list.environment(repeated_prediction, all.names = TRUE),
                        tolerance = 0))) {
    stop(sprintf("%s prediction is not deterministic", label), call. = FALSE)
  }
  if (!identical(unname(tools::md5sum("calibration-pristine.RData")),
                 pristine_hash)) {
    stop(sprintf("%s changed the pristine prediction checkpoint", label),
         call. = FALSE)
  }
  message(sprintf("%s calibration/event prediction contract passed", label))
  invisible(TRUE)
}

patch_profile_script <- function(path)
{
  lines <- readLines(path, warn = FALSE)
  replace <- function(pattern, value, expected) {
    hits <- grep(pattern, lines)
    if (length(hits) != expected) {
      stop(sprintf("%s expected %d matches for %s", basename(path), expected,
                   pattern), call. = FALSE)
    }
    lines[hits] <<- value
  }
  replace("^[[:space:]]*ngrid[[:space:]]*<-", "ngrid <- 5", 1L)
  replace("^[[:space:]]*w[[:space:]]*<-", "w <- seq(12,16,length=ngrid)", 1L)
  replace("^[[:space:]]*pl[[:space:]]*=", 'pl = "sequential"', 2L)
  replace("^[[:space:]]*ncor[[:space:]]*=", "ncor = 1", 2L)
  writeLines(lines, path, useBytes = TRUE)
}

run_profile_likelihood_smoke <- function(source_case, label)
{
  case_dir <- file.path(dirname(source_case), paste0(basename(source_case),
                                                    "-profile"))
  copy_case_tree(source_case, case_dir)
  on.exit(unlink(case_dir, recursive = TRUE, force = TRUE), add = TRUE)
  patch_advanced_runfiles(case_dir, bayes = FALSE)
  patch_profile_script(file.path(case_dir, "profile_ll.r"))
  old <- setwd(case_dir)
  on.exit(setwd(old), add = TRUE)
  run_batch_checked("runMPEM.r", "profile-calibration.out", restore = FALSE)
  if (!file.copy(".RData", "calibration-pristine.RData", overwrite = TRUE)) {
    stop("could not preserve profile calibration checkpoint", call. = FALSE)
  }
  pristine_hash <- unname(tools::md5sum("calibration-pristine.RData"))
  run_batch_checked("runMPEM_0.r", "profile-event.out", restore = TRUE)
  run_batch_checked("profile_ll.r", "profile-likelihood.out", restore = TRUE)
  workspace <- new.env(parent = globalenv())
  load(".RData", envir = workspace)
  for (name in c("pll", "pll_0")) {
    value <- get(name, envir = workspace, inherits = FALSE)
    if (length(value) != 5L || any(!is.finite(value))) {
      stop(sprintf("%s did not produce five finite %s values", label, name),
           call. = FALSE)
    }
  }
  if (length(workspace$pll_conv) != 5L || any(workspace$pll_conv != 0)) {
    stop(sprintf("%s profile optimizations did not all converge", label),
         call. = FALSE)
  }
  pdfs <- c("profile_ll.pdf", "profile_ll_w.pdf", "profile_ll_0.pdf")
  if (any(!file.exists(pdfs)) || any(file.info(pdfs)$size <= 0)) {
    stop(sprintf("%s profile plots are incomplete", label), call. = FALSE)
  }
  if (!identical(unname(tools::md5sum("calibration-pristine.RData")),
                 pristine_hash)) {
    stop(sprintf("%s changed the pristine profile checkpoint", label),
         call. = FALSE)
  }
  message(sprintf("%s reduced profile-likelihood workflow passed", label))
  invisible(TRUE)
}
