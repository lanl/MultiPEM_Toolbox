#!/usr/bin/env Rscript

options(warn = 1)

script_argument <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (!length(script_argument)) stop("testmpem.R must be run with Rscript", call. = FALSE)
script_path <- sub("^--file=", "", script_argument[[1L]])
script_dir <- dirname(normalizePath(script_path, winslash = "/", mustWork = TRUE))
repo_root <- dirname(script_dir)
image_name <- Sys.getenv("TEST_MPEM_IMAGE", "mpem/test-runtime:4.4.3")
results_root <- Sys.getenv("TEST_MPEM_RESULTS_ROOT", file.path(script_dir, "results"))
mapping_file <- file.path(script_dir, "test-groups.tsv")

usage <- function() {
  cat(paste0(
    "MultiPEM Docker test runner\n\n",
    "Usage:\n",
    "  Rscript testmpem.R doctor\n",
    "  Rscript testmpem.R build [--no-cache]\n",
    "  Rscript testmpem.R refresh-lock\n",
    "  Rscript testmpem.R list\n",
    "  Rscript testmpem.R run TEST [options]\n",
    "  Rscript testmpem.R jobs\n",
    "  Rscript testmpem.R status JOB\n",
    "  Rscript testmpem.R logs [-f] JOB\n",
    "  Rscript testmpem.R wait JOB\n",
    "  Rscript testmpem.R stop JOB\n",
    "  Rscript testmpem.R clean JOB\n",
    "  Rscript testmpem.R compare LEFT.RData RIGHT.RData [TOLERANCE]\n\n",
    "Run options:\n",
    "  --results DIR     New or existing job directory\n",
    "  --cpus N          Docker CPU and worker ceiling\n",
    "  --memory SIZE     Optional Docker memory ceiling\n",
    "  --foreground      Attach instead of running detached\n",
    "  --name NAME       Explicit container name\n\n",
    "TEST is global or relative to Test, for example IYDT/Seismic.\n"
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

image_id <- function(required = TRUE) {
  result <- docker_call(
    c("image", "inspect", "--format", "{{.Id}}", image_name),
    capture = TRUE, check = FALSE
  )
  if (result$status != 0L || !length(result$output)) {
    if (required) fail("image ", image_name, " is not built; run Rscript testmpem.R build")
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

volume_label <- Sys.getenv("TEST_MPEM_VOLUME_LABEL", "")
if (nzchar(volume_label) && !volume_label %in% c("z", "Z")) {
  fail("TEST_MPEM_VOLUME_LABEL must be empty, z, or Z")
}

volume_spec <- function(host_path, container_path, read_only = FALSE) {
  host_path <- normalizePath(host_path, winslash = "/", mustWork = TRUE)
  options <- character()
  if (read_only) options <- c(options, "ro")
  if (nzchar(volume_label)) options <- c(options, volume_label)
  suffix <- if (length(options)) paste0(":", paste(options, collapse = ",")) else ""
  paste0(host_path, ":", container_path, suffix)
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

test_components <- function(test) {
  test <- gsub("\\\\", "/", test)
  if (!nzchar(test) || startsWith(test, "/") || grepl("^[A-Za-z]:", test)) {
    fail("invalid test path: ", test)
  }
  if (identical(test, "global")) return(character())
  components <- strsplit(test, "/", fixed = TRUE)[[1L]]
  if (any(!nzchar(components)) || any(components %in% c(".", ".."))) {
    fail("invalid test path: ", test)
  }
  components
}

path_from_components <- function(...) do.call(file.path, as.list(c(...)))

read_mappings <- function() {
  if (!file.exists(mapping_file)) fail("missing test mapping: ", mapping_file)
  mappings <- read.delim(
    mapping_file, sep = "\t", quote = "", comment.char = "#",
    stringsAsFactors = FALSE, check.names = FALSE
  )
  expected <- c(
    "tests", "code", "fixtures", "secondary_code",
    "smoke_runfiles", "smoke_case", "smoke_data"
  )
  if (!identical(names(mappings), expected)) {
    fail(paste(
      "test-groups.tsv must contain tests, code, fixtures, secondary_code,",
      "smoke_runfiles, smoke_case, and smoke_data columns"
    ))
  }
  mappings
}

read_test_mapping <- function(group) {
  mappings <- read_mappings()
  match <- mappings[mappings$tests == group, , drop = FALSE]
  if (nrow(match) != 1L) fail("unsupported test group: ", group)
  list(
    code = match$code[[1L]], fixtures = match$fixtures[[1L]],
    secondary_code = match$secondary_code[[1L]],
    smoke_runfiles = match$smoke_runfiles[[1L]],
    smoke_case = match$smoke_case[[1L]],
    smoke_data = match$smoke_data[[1L]]
  )
}

list_test_names <- function() {
  root <- file.path(repo_root, "Test")
  files <- list.files(root, pattern = "^tests\\.r$", recursive = TRUE, full.names = FALSE)
  names <- gsub("\\\\", "/", dirname(files))
  names[names == "."] <- "global"
  sort(unique(names))
}

hash_tree <- function(path) {
  result <- docker_call(c(
    "run", "--rm", "--entrypoint", "sh",
    "--volume", volume_spec(path, "/snapshot", TRUE), image_name,
    "-c", "cd /snapshot && find . -type f -print0 | sort -z | xargs -0 -r sha256sum"
  ), capture = TRUE)
  result$output
}

prepare_new_job <- function(job_root, test) {
  components <- test_components(test)
  work_root <- file.path(job_root, "work")
  test_root <- file.path(work_root, "Test")
  metadata <- file.path(job_root, ".testmpem")
  dir.create(test_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(metadata, recursive = TRUE, showWarnings = FALSE)

  copy_contents(file.path(repo_root, "Test"), test_root)
  copy_contents(file.path(repo_root, "Code"), file.path(test_root, "Code"))

  if (length(components)) {
    group <- components[[1L]]
    mapping <- read_test_mapping(group)
    group_root <- file.path(test_root, group)
    copy_contents(
      file.path(repo_root, "Applications", "Code", mapping$code),
      file.path(group_root, "Code")
    )
    copy_contents(
      file.path(repo_root, "Applications", "Test", mapping$fixtures),
      file.path(group_root, "Test")
    )
    if (!identical(mapping$secondary_code, "-")) {
      copy_contents(
        file.path(repo_root, "Applications", "Code", mapping$secondary_code),
        file.path(group_root, paste0("Code-", mapping$secondary_code))
      )
    }
    if (length(components) >= 2L &&
        components[[2L]] %in% c(
          "Runfile", "Bayes", "Bayes-Backends", "MultiPhenomenology",
          "Parallel", "Prediction", "Profile", "Imputation"
        )) {
      smoke_root <- file.path(group_root, "Smoke", "Runfiles")
      copy_contents(file.path(repo_root, "Code"), file.path(smoke_root, "Code"))
      app_root <- file.path(smoke_root, mapping$smoke_runfiles)
      copy_contents(
        file.path(repo_root, "Applications", "Code", mapping$code),
        file.path(app_root, "Code")
      )
      copy_contents(
        file.path(repo_root, "Applications", "Data", mapping$smoke_data),
        file.path(app_root, "Data")
      )
      phenomenology <- strsplit(mapping$smoke_case, "/", fixed = TRUE)[[1L]][[1L]]
      phenomenon_code <- file.path(
        repo_root, "Applications", "Code", mapping$code, "Phenomenology"
      )
      if (dir.exists(phenomenon_code)) {
        copy_contents(
          phenomenon_code, file.path(app_root, phenomenology, "Code")
        )
      }
      copy_contents(
        file.path(repo_root, "Runfiles", mapping$smoke_runfiles,
                  mapping$smoke_case),
        file.path(app_root, mapping$smoke_case)
      )
      if (identical(components[[2L]], "MultiPhenomenology")) {
        if (!identical(group, "IYDT")) {
          fail("MultiPhenomenology staging is currently defined only for IYDT")
        }
        extra_cases <- c(
          "Optical/I-SUGAR-hob-0",
          "2-Phen-oc/I-SUGAR-hob-0"
        )
        for (case in extra_cases) {
          copy_contents(
            file.path(repo_root, "Runfiles", mapping$smoke_runfiles, case),
            file.path(app_root, case)
          )
        }
      }
    }
  }

  target_dir <- if (length(components)) {
    path_from_components(test_root, components)
  } else test_root
  if (!file.exists(file.path(target_dir, "tests.r"))) {
    fail("test disappeared while preparing workspace: ", test)
  }

  instrumentation <- run_external(
    file.path(R.home("bin"), if (.Platform$OS.type == "windows") "Rscript.exe" else "Rscript"),
    c(
      file.path(script_dir, "docker", "instrument-workers.R"),
      work_root,
      file.path(metadata, "worker-controls.tsv")
    ),
    capture = TRUE
  )
  if (length(instrumentation$output)) cat(paste0(instrumentation$output, "\n"), sep = "")
  write_value(test, file.path(metadata, "test"))
  write_value(
    format(Sys.time(), tz = "UTC", format = "%Y-%m-%dT%H:%M:%SZ"),
    file.path(metadata, "created")
  )
  writeLines(hash_tree(work_root), file.path(metadata, "initial-files.sha256"), useBytes = TRUE)
}

parse_run_options <- function(arguments) {
  result <- list(results = "", cpus = "", memory = "", foreground = FALSE, name = "")
  index <- 1L
  while (index <= length(arguments)) {
    option <- arguments[[index]]
    valued <- c("--results", "--cpus", "--memory", "--name")
    if (option %in% valued) {
      if (index == length(arguments)) fail(option, " needs a value")
      key <- switch(option,
        "--results" = "results", "--cpus" = "cpus",
        "--memory" = "memory", "--name" = "name"
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
  if (!length(arguments)) fail("run requires a test name")
  test <- gsub("\\\\", "/", arguments[[1L]])
  components <- test_components(test)
  options <- parse_run_options(arguments[-1L])

  source_dir <- if (length(components)) {
    path_from_components(repo_root, "Test", components)
  } else file.path(repo_root, "Test")
  if (!file.exists(file.path(source_dir, "tests.r"))) fail("unknown test: ", test)
  if (length(components)) read_test_mapping(components[[1L]])

  cpu_limit <- NA_real_
  if (nzchar(options$cpus)) {
    cpu_limit <- suppressWarnings(as.numeric(options$cpus))
    if (length(cpu_limit) != 1L || !is.finite(cpu_limit) || cpu_limit <= 0) {
      fail("--cpus must be a positive number")
    }
  }

  require_docker()
  runtime_id <- image_id()
  timestamp <- format(Sys.time(), tz = "UTC", format = "%Y%m%dT%H%M%SZ")
  slug <- gsub("[^a-z0-9.-]", "", gsub("[/_]", "-", tolower(test)))
  slug <- substr(if (nzchar(slug)) slug else "test", 1L, 48L)

  if (nzchar(options$results)) {
    job_root <- absolute_directory(options$results)
  } else {
    dir.create(results_root, recursive = TRUE, showWarnings = FALSE)
    job_root <- absolute_directory(file.path(
      results_root, paste(slug, timestamp, Sys.getpid(), sep = "-")
    ))
  }

  metadata <- file.path(job_root, ".testmpem")
  test_record <- file.path(metadata, "test")
  if (file.exists(test_record)) {
    recorded_test <- read_first_line(test_record)
    if (!identical(recorded_test, test)) {
      fail("results belong to ", recorded_test, ", not ", test)
    }
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
    prepare_new_job(job_root, test)
  }

  work_root <- file.path(job_root, "work")
  target_dir <- if (length(components)) {
    path_from_components(work_root, "Test", components)
  } else file.path(work_root, "Test")
  if (!file.exists(file.path(target_dir, "tests.r"))) fail("job workspace is incomplete")

  container_name <- options$name
  if (!nzchar(container_name)) {
    container_name <- paste("multipem-test", slug, timestamp, Sys.getpid(), sep = "-")
  }
  if (!grepl("^[A-Za-z0-9][A-Za-z0-9_.-]+$", container_name)) {
    fail("invalid Docker container name: ", container_name)
  }

  docker_arguments <- c(
    "run", "--name", container_name, "--init",
    "--label", "org.multipem.tests.managed=true",
    "--label", paste0("org.multipem.tests.test=", test),
    "--label", paste0("org.multipem.tests.results=", job_root),
    "--label", paste0("org.multipem.tests.image-id=", runtime_id)
  )
  owner <- file.info(work_root)
  if (.Platform$OS.type == "unix" && is.finite(owner$uid) && is.finite(owner$gid)) {
    docker_arguments <- c(docker_arguments, "--user", paste0(owner$uid, ":", owner$gid))
  }
  docker_arguments <- c(
    docker_arguments, "--env", "HOME=/tmp",
    "--volume", volume_spec(work_root, "/work")
  )
  if (nzchar(options$cpus)) {
    effective_workers <- max(1L, floor(cpu_limit))
    docker_arguments <- c(
      docker_arguments, "--cpus", options$cpus,
      "--env", paste0("TEST_MPEM_MAX_WORKERS=", effective_workers)
    )
  }
  if (nzchar(options$memory)) docker_arguments <- c(docker_arguments, "--memory", options$memory)
  if (!options$foreground) docker_arguments <- c(docker_arguments, "-d")
  docker_arguments <- c(docker_arguments, image_name, "run", test)

  write_value(container_name, file.path(metadata, "last-container"))
  write_value(runtime_id, file.path(metadata, "image-id"))
  cat("Results: ", job_root, "\n", sep = "")

  if (options$foreground) return(docker_call(docker_arguments, check = FALSE)$status)
  result <- docker_call(docker_arguments, capture = TRUE)
  container_id <- trimws(tail(result$output, 1L))
  cat("Started ", container_name, " (", substr(container_id, 1L, 12L), ")\n", sep = "")
  cat("Follow:  Rscript testmpem.R logs -f ", container_name, "\n", sep = "")
  cat("Wait:    Rscript testmpem.R wait ", container_name, "\n", sep = "")
  0L
}

require_job <- function(name) {
  if (!nzchar(name)) fail("a job name is required")
  exists <- docker_call(c("inspect", name), capture = TRUE, check = FALSE)
  if (exists$status != 0L) fail("unknown Docker test job: ", name)
  managed <- docker_call(c(
    "inspect", "--format", "{{index .Config.Labels \"org.multipem.tests.managed\"}}", name
  ), capture = TRUE)
  if (!identical(trimws(managed$output[[1L]]), "true")) {
    fail("container is not managed by the MultiPEM test runner: ", name)
  }
}

show_status <- function(name) {
  require_docker(); require_job(name)
  format <- paste0(
    "job={{.Name}} state={{.State.Status}} exit={{.State.ExitCode}} ",
    "started={{.State.StartedAt}} finished={{.State.FinishedAt}}{{println}}",
    "test={{index .Config.Labels \"org.multipem.tests.test\"}}{{println}}",
    "results={{index .Config.Labels \"org.multipem.tests.results\"}}"
  )
  result <- docker_call(c("inspect", "--format", format, name), capture = TRUE)
  cat(paste0(sub("^job=/", "job=", result$output), "\n"), sep = "")
  0L
}

list_jobs <- function() {
  require_docker()
  docker_call(c(
    "ps", "-a", "--filter", "label=org.multipem.tests.managed=true",
    "--format", paste0(
      "table {{.Names}}\\t{{.Status}}\\t",
      "{{.Label \"org.multipem.tests.test\"}}"
    )
  ))
  0L
}

show_logs <- function(arguments) {
  follow <- length(arguments) && arguments[[1L]] %in% c("-f", "--follow")
  if (follow) arguments <- arguments[-1L]
  if (length(arguments) != 1L) fail("logs requires one job name")
  require_docker(); require_job(arguments[[1L]])
  docker_call(c("logs", if (follow) "-f", arguments[[1L]]))
  0L
}

wait_job <- function(name) {
  require_docker(); require_job(name)
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
    image_name, "/opt/testmpem/docker/compare-workspaces.R",
    "/compare/left.RData", "/compare/right.RData", tolerance
  ), check = FALSE)$status
}

doctor <- function() {
  cat("R: ", R.version.string, "\n", sep = "")
  if (getRversion() < "4.0.0") fail("R 4.0 or newer is required")
  required <- file.path(repo_root, c("Test", "Applications", "Code"))
  if (!all(dir.exists(required))) {
    fail("Test-Docker must remain beside Test, Applications, and Code")
  }

  mappings <- read_mappings()
  if (anyDuplicated(mappings$tests)) fail("test-groups.tsv contains duplicate test groups")
  suites <- list_test_names()
  nested <- suites[suites != "global"]
  for (suite in nested) {
    group <- strsplit(suite, "/", fixed = TRUE)[[1L]][[1L]]
    mapping <- read_test_mapping(group)
    if (!dir.exists(file.path(repo_root, "Applications", "Code", mapping$code))) {
      fail("missing application code for test group ", group, ": ", mapping$code)
    }
    if (!dir.exists(file.path(repo_root, "Applications", "Test", mapping$fixtures))) {
      fail("missing fixtures for test group ", group, ": ", mapping$fixtures)
    }
    if (!identical(mapping$secondary_code, "-") &&
        !dir.exists(file.path(repo_root, "Applications", "Code", mapping$secondary_code))) {
      fail("missing secondary application code for test group ", group, ": ",
           mapping$secondary_code)
    }
    smoke_paths <- c(
      file.path(repo_root, "Runfiles", mapping$smoke_runfiles,
                mapping$smoke_case),
      file.path(repo_root, "Applications", "Data", mapping$smoke_data)
    )
    if (!all(dir.exists(smoke_paths))) {
      fail("missing runfile smoke inputs for test group ", group)
    }
  }

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
  if (is.na(id)) cat("Runtime image: not built (run Rscript testmpem.R build)\n")
  else cat("Runtime image: ", id, "\n", sep = "")
  if (nzchar(volume_label)) cat("SELinux volume label: ", volume_label, "\n", sep = "")
  cat("Registered tests: ", length(suites), "\n", sep = "")
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
    "list" = { cat(paste0(list_test_names(), "\n"), sep = ""); 0L },
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
