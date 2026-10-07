#!/usr/bin/env Rscript

options(warn = 1)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) {
  stop("mpem.R must be run with Rscript", call. = FALSE)
}
script_path <- sub("^--file=", "", script_argument[[1L]])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- dirname(script_dir)
image_name <- Sys.getenv("MPEM_IMAGE", "mpem/runtime:4.4.3")
results_root <- Sys.getenv("MPEM_RESULTS_ROOT", file.path(script_dir, "results"))
mapping_file <- file.path(script_dir, "applications.tsv")

usage <- function() {
  cat(paste0(
    "MultiPEM Docker job runner\n\n",
    "Usage:\n",
    "  Rscript mpem.R doctor\n",
    "  Rscript mpem.R build [--no-cache]\n",
    "  Rscript mpem.R refresh-lock\n",
    "  Rscript mpem.R list\n",
    "  Rscript mpem.R run ANALYSIS [options]\n",
    "  Rscript mpem.R jobs\n",
    "  Rscript mpem.R status JOB\n",
    "  Rscript mpem.R logs [-f] JOB\n",
    "  Rscript mpem.R wait JOB\n",
    "  Rscript mpem.R stop JOB\n",
    "  Rscript mpem.R clean JOB\n",
    "  Rscript mpem.R compare LEFT.RData RIGHT.RData [TOLERANCE]\n",
    "  Rscript mpem.R verify [--global-only|--apps-only] [REFERENCE_IMAGE]\n\n",
    "Run options:\n",
    "  --stage full|calibration|event  Stage to run (default: full)\n",
    "  --results DIR                  New or existing job directory\n",
    "  --workspace FILE               Pristine calibration for a new event job\n",
    "  --cpus N                       Docker CPU and worker ceiling\n",
    "  --memory SIZE                  Optional Docker memory ceiling\n",
    "  --foreground                   Attach instead of running detached\n",
    "  --name NAME                    Explicit container name\n\n",
    "ANALYSIS is relative to Runfiles, for example:\n",
    "  IYDT-gsrp/Seismic/I-SUGAR-hob\n"
  ))
}

fail <- function(...) stop(paste0(...), call. = FALSE)

quote_arguments <- function(arguments) {
  vapply(as.character(arguments), shQuote, character(1L), USE.NAMES = FALSE)
}

run_external <- function(command, arguments = character(), capture = FALSE,
                         check = TRUE) {
  quoted <- quote_arguments(arguments)
  if (capture) {
    output <- suppressWarnings(system2(command, quoted, stdout = TRUE, stderr = TRUE))
    status <- attr(output, "status", exact = TRUE)
    if (is.null(status)) status <- 0L
  } else {
    status <- suppressWarnings(system2(command, quoted, stdout = "", stderr = ""))
    output <- character()
  }
  status <- as.integer(status)
  if (check && status != 0L) {
    if (length(output)) cat(paste0(output, "\n"), sep = "", file = stderr())
    fail(command, " exited with status ", status)
  }
  list(status = status, output = output)
}

docker_call <- function(arguments, capture = FALSE, check = TRUE) {
  run_external("docker", arguments, capture = capture, check = check)
}

require_docker <- function() {
  if (!nzchar(Sys.which("docker"))) fail("docker is not installed or not on PATH")
  result <- docker_call("info", capture = TRUE, check = FALSE)
  if (result$status != 0L) fail("the Docker engine is not available")
}

image_id <- function(required = TRUE, name = image_name) {
  result <- docker_call(
    c("image", "inspect", "--format", "{{.Id}}", name),
    capture = TRUE, check = FALSE
  )
  if (result$status != 0L || !length(result$output)) {
    if (required) fail("image ", name, " is not built; run Rscript mpem.R build")
    return(NA_character_)
  }
  trimws(result$output[[1L]])
}

absolute_directory <- function(path) {
  if (!dir.exists(path) && !dir.create(path, recursive = TRUE)) {
    fail("cannot create directory: ", path)
  }
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

absolute_file <- function(path) {
  if (!file.exists(path) || dir.exists(path)) fail("file does not exist: ", path)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

read_first_line <- function(path) {
  value <- readLines(path, n = 1L, warn = FALSE)
  if (length(value)) value[[1L]] else ""
}

write_value <- function(value, path) writeLines(as.character(value), path, useBytes = TRUE)

volume_label <- Sys.getenv("MPEM_VOLUME_LABEL", "")
if (nzchar(volume_label) && !volume_label %in% c("z", "Z")) {
  fail("MPEM_VOLUME_LABEL must be empty, z, or Z")
}

volume_spec <- function(host_path, container_path, read_only = FALSE) {
  host_path <- normalizePath(host_path, winslash = "/", mustWork = TRUE)
  options <- character()
  if (read_only) options <- c(options, "ro")
  if (nzchar(volume_label)) options <- c(options, volume_label)
  suffix <- if (length(options)) paste0(":", paste(options, collapse = ",")) else ""
  paste0(host_path, ":", container_path, suffix)
}

read_application_mapping <- function(application) {
  if (!file.exists(mapping_file)) fail("missing application mapping: ", mapping_file)
  mappings <- read.delim(
    mapping_file, sep = "\t", quote = "", comment.char = "#",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  if (!identical(names(mappings), c("runfiles", "code", "data"))) {
    fail("applications.tsv must contain runfiles, code, and data columns")
  }
  match <- mappings[mappings$runfiles == application, , drop = FALSE]
  if (nrow(match) != 1L) fail("unsupported application: ", application)
  list(code = match$code[[1L]], data = match$data[[1L]])
}

copy_contents <- function(source, destination) {
  if (!dir.exists(source)) fail("source directory does not exist: ", source)
  if (!dir.exists(destination) && !dir.create(destination, recursive = TRUE)) {
    fail("cannot create directory: ", destination)
  }
  entries <- list.files(source, all.files = TRUE, no.. = TRUE, full.names = TRUE)
  if (length(entries)) {
    copied <- file.copy(
      entries, destination, overwrite = FALSE, recursive = TRUE,
      copy.mode = TRUE, copy.date = TRUE
    )
    if (!all(copied)) fail("failed to copy directory contents from ", source)
  }
}

copy_one_file <- function(source, destination, overwrite = TRUE) {
  if (!file.copy(
    source, destination, overwrite = overwrite, copy.mode = TRUE, copy.date = TRUE
  )) fail("failed to copy ", source, " to ", destination)
}

hash_tree <- function(path) {
  result <- docker_call(c(
    "run", "--rm", "--entrypoint", "sh",
    "--volume", volume_spec(path, "/snapshot", TRUE), image_name,
    "-c", "cd /snapshot && find . -type f -print0 | sort -z | xargs -0 -r sha256sum"
  ), capture = TRUE)
  result$output
}

hash_file <- function(path, displayed_name) {
  result <- docker_call(c(
    "run", "--rm", "--entrypoint", "sha256sum",
    "--volume", volume_spec(path, "/input", TRUE), image_name, "/input"
  ), capture = TRUE)
  sub("  /input$", paste0("  ", displayed_name), result$output)
}

verify_checkpoint <- function(checkpoint_dir) {
  checkpoint <- file.path(checkpoint_dir, "calibration.RData")
  checksum <- file.path(checkpoint_dir, "calibration.RData.sha256")
  if (!file.exists(checkpoint)) {
    fail("event stage requires a pristine calibration checkpoint; use --workspace when creating the job")
  }
  if (!file.exists(checksum)) {
    fail("event stage requires the calibration checkpoint checksum")
  }
  expected <- trimws(readLines(checksum, n = 1L, warn = FALSE))
  actual <- trimws(hash_file(checkpoint, "calibration.RData")[[1L]])
  if (length(expected) != 1L || !nzchar(expected) ||
      !identical(expected, actual)) {
    fail("calibration checkpoint checksum does not match; restore the pristine checkpoint")
  }
  invisible(TRUE)
}

analysis_components <- function(analysis) {
  analysis <- gsub("\\\\", "/", analysis)
  if (!nzchar(analysis) || startsWith(analysis, "/") || grepl("^[A-Za-z]:", analysis)) {
    fail("invalid analysis path: ", analysis)
  }
  components <- strsplit(analysis, "/", fixed = TRUE)[[1L]]
  if (any(!nzchar(components)) || any(components %in% c(".", ".."))) {
    fail("invalid analysis path: ", analysis)
  }
  components
}

path_from_components <- function(...) do.call(file.path, as.list(c(...)))

prepare_new_job <- function(job_root, analysis, workspace_file) {
  components <- analysis_components(analysis)
  application <- components[[1L]]
  relative_analysis <- components[-1L]
  mapping <- read_application_mapping(application)
  work_root <- file.path(job_root, "work")
  application_root <- file.path(work_root, application)
  target_dir <- path_from_components(application_root, relative_analysis)
  checkpoint_dir <- file.path(job_root, ".mpem", "checkpoints")

  dir.create(application_root, recursive = TRUE, showWarnings = FALSE)
  copy_contents(file.path(repo_root, "Code"), file.path(work_root, "Code"))
  copy_contents(file.path(repo_root, "Runfiles", application), application_root)
  copy_contents(
    file.path(repo_root, "Applications", "Code", mapping$code),
    file.path(application_root, "Code")
  )
  copy_contents(
    file.path(repo_root, "Applications", "Data", mapping$data),
    file.path(application_root, "Data")
  )

  phenomenology_code <- file.path(application_root, "Code", "Phenomenology")
  if (dir.exists(phenomenology_code)) {
    children <- list.dirs(application_root, full.names = TRUE, recursive = FALSE)
    children <- children[!basename(children) %in% c("Code", "Data", "Plots")]
    for (child in children) {
      destination <- file.path(child, "Code")
      if (!file.exists(destination)) copy_contents(phenomenology_code, destination)
    }
  }

  if (!dir.exists(target_dir)) fail("analysis disappeared while preparing workspace: ", analysis)
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  if (nzchar(workspace_file)) {
    workspace_file <- absolute_file(workspace_file)
    working_workspace <- file.path(target_dir, ".RData")
    checkpoint <- file.path(checkpoint_dir, "calibration.RData")
    copy_one_file(workspace_file, working_workspace)
    Sys.chmod(working_workspace, mode = "0644")
    copy_one_file(workspace_file, checkpoint)
    Sys.chmod(checkpoint, mode = "0444")
  }

  metadata <- file.path(job_root, ".mpem")
  instrumentation <- run_external(
    file.path(
      R.home("bin"),
      if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"
    ),
    c(
      file.path(script_dir, "docker", "instrument-workers.R"),
      work_root,
      file.path(metadata, "worker-controls.tsv")
    ),
    capture = TRUE
  )
  if (length(instrumentation$output)) {
    cat(paste0(instrumentation$output, "\n"), sep = "")
  }
  write_value(analysis, file.path(metadata, "analysis"))
  write_value(
    format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
    file.path(metadata, "created")
  )
  writeLines(hash_tree(work_root), file.path(metadata, "initial-files.sha256"), useBytes = TRUE)
  checkpoint <- file.path(checkpoint_dir, "calibration.RData")
  if (file.exists(checkpoint)) {
    writeLines(
      hash_file(checkpoint, "calibration.RData"),
      file.path(checkpoint_dir, "calibration.RData.sha256"), useBytes = TRUE
    )
  }
}

list_analyses <- function() {
  root <- file.path(repo_root, "Runfiles")
  files <- list.files(
    root, pattern = "^runMPEM(_0)?\\.r$", recursive = TRUE, full.names = FALSE
  )
  analyses <- sort(unique(gsub("\\\\", "/", dirname(files))))
  cat(paste0(analyses, "\n"), sep = "")
  0L
}

parse_run_options <- function(arguments) {
  result <- list(
    stage = "full", results = "", workspace = "", cpus = "", memory = "",
    foreground = FALSE, name = ""
  )
  index <- 1L
  while (index <= length(arguments)) {
    option <- arguments[[index]]
    valued <- c("--stage", "--results", "--workspace", "--cpus", "--memory", "--name")
    if (option %in% valued) {
      if (index == length(arguments)) fail(option, " needs a value")
      key <- switch(option,
        "--stage" = "stage", "--results" = "results", "--workspace" = "workspace",
        "--cpus" = "cpus", "--memory" = "memory", "--name" = "name"
      )
      result[[key]] <- arguments[[index + 1L]]
      index <- index + 2L
    } else if (option == "--foreground") {
      result$foreground <- TRUE
      index <- index + 1L
    } else fail("unknown run option: ", option)
  }
  result
}

run_job <- function(arguments) {
  if (!length(arguments)) fail("run requires an analysis")
  analysis <- gsub("\\\\", "/", arguments[[1L]])
  components <- analysis_components(analysis)
  options <- parse_run_options(arguments[-1L])
  if (!options$stage %in% c("full", "calibration", "event")) {
    fail("invalid stage: ", options$stage)
  }
  if (nzchar(options$workspace) && options$stage != "event") {
    fail("--workspace is only valid with --stage event")
  }
  cpu_limit <- NA_real_
  if (nzchar(options$cpus)) {
    cpu_limit <- suppressWarnings(as.numeric(options$cpus))
    if (length(cpu_limit) != 1L || !is.finite(cpu_limit) || cpu_limit <= 0) {
      fail("--cpus must be a positive number")
    }
  }

  source_dir <- path_from_components(repo_root, "Runfiles", components)
  if (!dir.exists(source_dir)) fail("unknown analysis: ", analysis)
  calibration_deck <- file.exists(file.path(source_dir, "runMPEM.r"))
  event_deck <- file.exists(file.path(source_dir, "runMPEM_0.r"))
  if (!calibration_deck && !event_deck) fail("analysis has no runMPEM input deck: ", analysis)
  if (options$stage == "event" && !event_deck) fail("analysis has no runMPEM_0.r event deck")

  require_docker()
  runtime_id <- image_id()
  timestamp <- format(Sys.time(), tz = "UTC", format = "%Y%m%dT%H%M%SZ")
  slug <- gsub("[^a-z0-9.-]", "", gsub("[/_]", "-", tolower(analysis)))
  slug <- substr(if (nzchar(slug)) slug else "analysis", 1L, 48L)

  if (nzchar(options$results)) {
    job_root <- absolute_directory(options$results)
  } else {
    dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
    job_root <- absolute_directory(file.path(
      results_root, paste(slug, timestamp, Sys.getpid(), sep = "-")
    ))
  }

  metadata <- file.path(job_root, ".mpem")
  analysis_record <- file.path(metadata, "analysis")
  if (file.exists(analysis_record)) {
    recorded_analysis <- read_first_line(analysis_record)
    if (!identical(recorded_analysis, analysis)) {
      fail("results belong to ", recorded_analysis, ", not ", analysis)
    }
    if (nzchar(options$workspace)) fail("--workspace is only valid when creating a job")
    image_record <- file.path(metadata, "image-id")
    if (file.exists(image_record) && !identical(read_first_line(image_record), runtime_id)) {
      fail("job was created with a different runtime image; use a new results directory")
    }
    if (!file.exists(file.path(metadata, "worker-controls.tsv"))) {
      fail("job is missing CPU-control metadata; use a new results directory")
    }
  } else {
    existing <- list.files(job_root, all.files = TRUE, no.. = TRUE)
    if (length(existing)) fail("new results directory is not empty: ", job_root)
    prepare_new_job(job_root, analysis, options$workspace)
  }

  work_root <- file.path(job_root, "work")
  target_dir <- path_from_components(work_root, components)
  if (!dir.exists(target_dir)) fail("job workspace is incomplete")
  checkpoint_dir <- file.path(metadata, "checkpoints")
  dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  checkpoint <- file.path(checkpoint_dir, "calibration.RData")
  if (options$stage == "event") verify_checkpoint(checkpoint_dir)

  container_name <- options$name
  if (!nzchar(container_name)) {
    container_name <- paste("multipem", slug, timestamp, Sys.getpid(), sep = "-")
  }
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]+$", container_name)) {
    fail("invalid Docker container name: ", container_name)
  }

  checkpoint_mount <- volume_spec(checkpoint_dir, "/checkpoints", options$stage == "event")
  docker_arguments <- c(
    "run", "--name", container_name, "--init",
    "--label", "org.multipem.managed=true",
    "--label", paste0("org.multipem.analysis=", analysis),
    "--label", paste0("org.multipem.results=", job_root),
    "--label", paste0("org.multipem.stage=", options$stage),
    "--label", paste0("org.multipem.image-id=", runtime_id)
  )
  owner <- file.info(work_root)
  if (.Platform$OS.type == "unix" && is.finite(owner$uid) && is.finite(owner$gid)) {
    docker_arguments <- c(docker_arguments, "--user", paste0(owner$uid, ":", owner$gid))
  }
  docker_arguments <- c(
    docker_arguments, "--env", "HOME=/tmp",
    "--volume", volume_spec(work_root, "/work/Runfiles"),
    "--volume", checkpoint_mount
  )
  if (nzchar(options$cpus)) {
    effective_workers <- max(1L, floor(cpu_limit))
    docker_arguments <- c(
      docker_arguments,
      "--cpus", options$cpus,
      "--env", paste0("MPEM_MAX_WORKERS=", effective_workers)
    )
  }
  if (nzchar(options$memory)) docker_arguments <- c(docker_arguments, "--memory", options$memory)
  if (!options$foreground) docker_arguments <- c(docker_arguments, "-d")
  docker_arguments <- c(docker_arguments, image_name, "run", analysis, options$stage)

  write_value(container_name, file.path(metadata, "last-container"))
  write_value(runtime_id, file.path(metadata, "image-id"))
  write_value(options$stage, file.path(metadata, "last-stage"))

  if (options$foreground) return(docker_call(docker_arguments, check = FALSE)$status)
  result <- docker_call(docker_arguments, capture = TRUE)
  container_id <- trimws(tail(result$output, 1L))
  cat("Started ", container_name, " (", substr(container_id, 1L, 12L), ")\n", sep = "")
  cat("Results: ", job_root, "\n", sep = "")
  cat("Follow:  Rscript mpem.R logs -f ", container_name, "\n", sep = "")
  cat("Wait:    Rscript mpem.R wait ", container_name, "\n", sep = "")
  0L
}

require_job <- function(name) {
  if (!nzchar(name)) fail("a job name is required")
  exists <- docker_call(c("inspect", name), capture = TRUE, check = FALSE)
  if (exists$status != 0L) fail("unknown Docker job: ", name)
  managed <- docker_call(c(
    "inspect", "--format", "{{index .Config.Labels \"org.multipem.managed\"}}", name
  ), capture = TRUE)
  if (!identical(trimws(managed$output[[1L]]), "true")) {
    fail("container is not managed by MultiPEM: ", name)
  }
}

show_status <- function(name) {
  require_docker()
  require_job(name)
  format <- paste0(
    "job={{.Name}} state={{.State.Status}} exit={{.State.ExitCode}} ",
    "started={{.State.StartedAt}} finished={{.State.FinishedAt}}{{println}}",
    "analysis={{index .Config.Labels \"org.multipem.analysis\"}}{{println}}",
    "stage={{index .Config.Labels \"org.multipem.stage\"}}{{println}}",
    "results={{index .Config.Labels \"org.multipem.results\"}}"
  )
  result <- docker_call(c("inspect", "--format", format, name), capture = TRUE)
  cat(paste0(sub("^job=/", "job=", result$output), "\n"), sep = "")
  0L
}

list_jobs <- function() {
  require_docker()
  docker_call(c(
    "ps", "-a", "--filter", "label=org.multipem.managed=true",
    "--format", paste0(
      "table {{.Names}}\\t{{.Status}}\\t{{.Label \"org.multipem.stage\"}}\\t",
      "{{.Label \"org.multipem.analysis\"}}"
    )
  ))
  0L
}

show_logs <- function(arguments) {
  follow <- length(arguments) && arguments[[1L]] %in% c("-f", "--follow")
  if (follow) arguments <- arguments[-1L]
  if (length(arguments) != 1L) fail("logs requires one job name")
  require_docker()
  require_job(arguments[[1L]])
  docker_call(c("logs", if (follow) "-f", arguments[[1L]]))
  0L
}

wait_job <- function(name) {
  require_docker()
  require_job(name)
  result <- docker_call(c("wait", name), capture = TRUE)
  exit_code <- suppressWarnings(as.integer(trimws(tail(result$output, 1L))))
  show_status(name)
  if (is.na(exit_code) || exit_code > 255L) 1L else exit_code
}

stop_job <- function(name) {
  require_docker(); require_job(name); docker_call(c("stop", name)); 0L
}

clean_job <- function(name) {
  require_docker(); require_job(name); docker_call(c("rm", name))
  cat("Removed container ", name, "; result files were retained.\n", sep = "")
  0L
}

compare_workspaces <- function(arguments) {
  if (length(arguments) < 2L || length(arguments) > 3L) {
    fail("compare requires two .RData files and an optional tolerance")
  }
  left <- absolute_file(arguments[[1L]])
  right <- absolute_file(arguments[[2L]])
  tolerance <- if (length(arguments) == 3L) arguments[[3L]] else "1e-12"
  require_docker(); image_id()
  docker_call(c(
    "run", "--rm", "--entrypoint", "Rscript",
    "--volume", volume_spec(left, "/compare/left.RData", TRUE),
    "--volume", volume_spec(right, "/compare/right.RData", TRUE),
    image_name, "/opt/multipem/docker/compare-workspaces.R",
    "/compare/left.RData", "/compare/right.RData", tolerance
  ), check = FALSE)$status
}

run_to_file <- function(arguments, path) {
  result <- docker_call(arguments, capture = TRUE, check = FALSE)
  writeLines(result$output, path, useBytes = TRUE)
  if (result$status != 0L) fail("verification command failed; see ", path)
}

verification_user_arguments <- function(path) {
  owner <- file.info(path)
  if (.Platform$OS.type == "unix" && is.finite(owner$uid) && is.finite(owner$gid)) {
    c("--user", paste0(owner$uid, ":", owner$gid))
  } else character()
}

run_verification_test <- function(image, root, workdir, source, destination, log) {
  script <- paste(
    "set -eu",
    paste0("mkdir -p ", destination),
    paste0("cp -R ", source, "/. ", destination, "/"),
    paste0("cd ", workdir),
    "Rscript /run-test.R tests.r .RData",
    sep = "; "
  )
  run_to_file(c(
    "run", "--rm", verification_user_arguments(root),
    "--env", "HOME=/tmp", "--entrypoint", "sh",
    "--volume", volume_spec(root, "/verification"),
    "--volume", volume_spec(
      file.path(script_dir, "docker", "run-test.R"), "/run-test.R", TRUE
    ),
    image, "-c", script
  ), log)
}

verify_equivalence <- function(arguments) {
  mode <- "all"
  if (length(arguments) && arguments[[1L]] == "--global-only") {
    mode <- "global"; arguments <- arguments[-1L]
  } else if (length(arguments) && arguments[[1L]] == "--apps-only") {
    mode <- "apps"; arguments <- arguments[-1L]
  }
  if (length(arguments) > 1L) fail("too many verify arguments")
  reference_image <- if (length(arguments)) arguments[[1L]] else "mpem/accepted-baseline"
  require_docker(); image_id(); image_id(name = reference_image)

  verification_root <- tempfile("multipem-equivalence.")
  dir.create(verification_root)
  on.exit(unlink(verification_root, recursive = TRUE, force = TRUE), add = TRUE)
  reference_dir <- file.path(verification_root, "reference")
  runtime_dir <- file.path(verification_root, "runtime")
  dir.create(reference_dir); dir.create(runtime_dir)

  if (mode != "apps") {
    for (suite in c("reference", "runtime")) {
      global_root <- file.path(verification_root, suite, "Global")
      dir.create(global_root, recursive = TRUE)
      copy_one_file(
        file.path(repo_root, "Test", "tests.r"),
        file.path(global_root, "tests.r")
      )
    }
    cat("Running global verification in reference image...\n")
    run_verification_test(
      reference_image, file.path(reference_dir, "Global"),
      "/verification", "/opt/multipem/Code", "/verification/Code",
      file.path(reference_dir, "Global", "tests.out")
    )
    cat("Running global verification in runtime image...\n")
    run_verification_test(
      image_name, file.path(runtime_dir, "Global"),
      "/verification", "/opt/multipem/Code", "/verification/Code",
      file.path(runtime_dir, "Global", "tests.out")
    )
    status <- compare_workspaces(c(
      file.path(reference_dir, "Global", ".RData"),
      file.path(runtime_dir, "Global", ".RData"), "1e-12"
    ))
    if (status != 0L) fail("global workspace comparison failed")
  }

  if (mode != "global") {
    phenomena <- c(
      "Seismic", "Acoustic", "Crater", "Optical", "Data", "Transform", "Prior"
    )
    for (suite in c("reference", "runtime")) {
      test_root <- file.path(verification_root, suite, "IYDT")
      dir.create(test_root, recursive = TRUE, showWarnings = FALSE)
      copy_contents(
        file.path(repo_root, "Applications", "Test", "IYDT"),
        file.path(test_root, "Test")
      )
      for (phenomenon in phenomena) {
        destination <- file.path(test_root, phenomenon)
        dir.create(destination, recursive = TRUE, showWarnings = FALSE)
        copy_one_file(
          file.path(repo_root, "Test", "IYDT", phenomenon, "tests.r"),
          file.path(destination, "tests.r")
        )
      }
    }
    for (phenomenon in phenomena) {
      cat("Running IYDT ", phenomenon, " verification in reference image...\n", sep = "")
      run_verification_test(
        reference_image, reference_dir,
        paste0("/verification/IYDT/", phenomenon),
        "/opt/multipem/Applications/Code/IYDT-gsrp",
        "/verification/IYDT/Code",
        file.path(reference_dir, "IYDT", phenomenon, "tests.out")
      )
      cat("Running IYDT ", phenomenon, " verification in runtime image...\n", sep = "")
      run_verification_test(
        image_name, runtime_dir,
        paste0("/verification/IYDT/", phenomenon),
        "/opt/multipem/Applications/Code/IYDT-gsrp",
        "/verification/IYDT/Code",
        file.path(runtime_dir, "IYDT", phenomenon, "tests.out")
      )
      status <- compare_workspaces(c(
        file.path(reference_dir, "IYDT", phenomenon, ".RData"),
        file.path(runtime_dir, "IYDT", phenomenon, ".RData"), "1e-12"
      ))
      if (status != 0L) fail(phenomenon, " workspace comparison failed")
    }
  }
  cat("Reference and runtime ", mode, " calculations are numerically equivalent.\n", sep = "")
  0L
}

doctor <- function() {
  cat("R: ", R.version.string, "\n", sep = "")
  if (getRversion() < "4.0.0") fail("R 4.0 or newer is required")
  required <- file.path(repo_root, c("Runfiles", "Applications", "Code", "Test"))
  if (!all(dir.exists(required))) {
    fail("Runfiles-Docker must remain beside Runfiles, Applications, Code, and Test")
  }
  read_application_mapping("IYDT-gsrp")
  require_docker()
  info <- docker_call(c(
    "info", "--format", "{{.ServerVersion}}|{{.OSType}}|{{.Architecture}}"
  ), capture = TRUE)
  fields <- strsplit(trimws(info$output[[1L]]), "|", fixed = TRUE)[[1L]]
  if (length(fields) >= 2L && fields[[2L]] != "linux") {
    fail("Docker must be configured for Linux containers")
  }
  cat("Docker: ", paste(fields, collapse = " "), "\n", sep = "")
  id <- image_id(required = FALSE)
  if (is.na(id)) cat("Runtime image: not built (run Rscript mpem.R build)\n")
  else cat("Runtime image: ", id, "\n", sep = "")
  if (nzchar(volume_label)) cat("SELinux volume label: ", volume_label, "\n", sep = "")
  cat("Repository layout and required host commands are available.\n")
  0L
}

build_image <- function(arguments) {
  if (length(arguments) > 1L || (length(arguments) && arguments[[1L]] != "--no-cache")) {
    fail("build accepts only --no-cache")
  }
  require_docker()
  docker_arguments <- c("build", "--progress=plain")
  if (length(arguments)) docker_arguments <- c(docker_arguments, "--no-cache")
  docker_arguments <- c(
    docker_arguments, "-f", file.path(script_dir, "Dockerfile.runtime"),
    "-t", image_name, repo_root
  )
  docker_call(docker_arguments)
  0L
}

refresh_lock <- function() {
  require_docker(); image_id()
  docker_call(c(
    "run", "--rm", "--entrypoint", "Rscript",
    "--volume", volume_spec(script_dir, "/out"),
    "--volume", volume_spec(
      file.path(script_dir, "docker", "create-package-lock.R"),
      "/create-package-lock.R", TRUE
    ), image_name, "/create-package-lock.R"
  ))
  0L
}

main <- function() {
  arguments <- commandArgs(trailingOnly = TRUE)
  command <- if (length(arguments)) arguments[[1L]] else "help"
  rest <- if (length(arguments) > 1L) arguments[-1L] else character()
  switch(command,
    "doctor" = doctor(),
    "build" = build_image(rest),
    "refresh-lock" = refresh_lock(),
    "list" = list_analyses(),
    "run" = run_job(rest),
    "jobs" = list_jobs(),
    "status" = {
      if (length(rest) != 1L) fail("status requires one job name")
      show_status(rest[[1L]])
    },
    "logs" = show_logs(rest),
    "wait" = {
      if (length(rest) != 1L) fail("wait requires one job name")
      wait_job(rest[[1L]])
    },
    "stop" = {
      if (length(rest) != 1L) fail("stop requires one job name")
      stop_job(rest[[1L]])
    },
    "clean" = {
      if (length(rest) != 1L) fail("clean requires one job name")
      clean_job(rest[[1L]])
    },
    "compare" = compare_workspaces(rest),
    "verify" = verify_equivalence(rest),
    "help" = { usage(); 0L },
    "-h" = { usage(); 0L },
    "--help" = { usage(); 0L },
    { usage(); fail("unknown command: ", command) }
  )
}

status <- tryCatch(
  main(),
  error = function(error) {
    cat("Error: ", conditionMessage(error), "\n", sep = "", file = stderr())
    1L
  }
)
quit(save = "no", status = as.integer(status), runLast = FALSE)
