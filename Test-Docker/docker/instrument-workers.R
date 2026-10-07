#!/usr/bin/env Rscript

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  stop("Usage: instrument-workers.R WORK_ROOT REPORT_FILE", call. = FALSE)
}

work_root <- normalizePath(arguments[[1L]], winslash = "/", mustWork = TRUE)
report_file <- arguments[[2L]]
files <- list.files(
  work_root, pattern = "\\.[Rr]$", recursive = TRUE,
  full.names = TRUE, all.files = TRUE
)

pool_variables <- c(
  "ncor", "ncor_map", "ncor_mc", "ncor_smc",
  "ncores_mle", "ncores_map", "ncores_mc", "ncores_smc"
)
records <- list()
record_index <- 0L

record_change <- function(file, line_number, variable, control) {
  record_index <<- record_index + 1L
  records[[record_index]] <<- data.frame(
    file = substring(normalizePath(file, winslash = "/"), nchar(work_root) + 2L),
    line = line_number,
    variable = variable,
    control = control,
    stringsAsFactors = FALSE
  )
}

for (file in files) {
  lines <- readLines(file, warn = FALSE)
  output <- character()
  changed <- FALSE
  for (line_number in seq_along(lines)) {
    line <- lines[[line_number]]
    for (variable in pool_variables) {
      pattern <- paste0("workers[[:space:]]*=[[:space:]]*", variable, "\\b")
      if (grepl(pattern, line, perl = TRUE)) {
        replacement <- paste0(
          "workers=getOption(\"testmpem.limit.workers\", ",
          "function(value, label, mode) value)(",
          variable, ", \"", variable, "\", \"pool\")"
        )
        line <- gsub(pattern, replacement, line, perl = TRUE)
        record_change(file, line_number, variable, "worker-pool")
        changed <- TRUE
      }
    }
    output <- c(output, line)
  }
  if (changed) writeLines(output, file, useBytes = TRUE)
}

dir.create(dirname(report_file), recursive = TRUE, showWarnings = FALSE)
if (length(records)) {
  report <- do.call(rbind, records)
} else {
  report <- data.frame(
    file = character(), line = integer(), variable = character(),
    control = character(), stringsAsFactors = FALSE
  )
}
write.table(report, report_file, sep = "\t", row.names = FALSE, quote = FALSE)
cat(sprintf(
  "Instrumented %d worker controls in %d files.\n",
  nrow(report), length(unique(report$file))
))
