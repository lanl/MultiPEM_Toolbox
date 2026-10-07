#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Usage: compare-workspaces.R LEFT.RData RIGHT.RData [tolerance]", call. = FALSE)
}

left_file <- normalizePath(args[[1L]], mustWork = TRUE)
right_file <- normalizePath(args[[2L]], mustWork = TRUE)
tolerance <- if (length(args) == 3L) as.numeric(args[[3L]]) else 1e-12
if (!is.finite(tolerance) || tolerance < 0) {
  stop("Tolerance must be a non-negative finite number", call. = FALSE)
}

left <- new.env(parent = emptyenv())
right <- new.env(parent = emptyenv())
load(left_file, envir = left)
load(right_file, envir = right)

left_names <- sort(ls(left, all.names = TRUE))
right_names <- sort(ls(right, all.names = TRUE))
if (!identical(left_names, right_names)) {
  only_left <- setdiff(left_names, right_names)
  only_right <- setdiff(right_names, left_names)
  if (length(only_left)) cat("Only in left: ", paste(only_left, collapse = ", "), "\n", sep = "")
  if (length(only_right)) cat("Only in right: ", paste(only_right, collapse = ", "), "\n", sep = "")
  quit(status = 1L)
}

failures <- character()
for (name in left_names) {
  result <- all.equal(
    get(name, envir = left, inherits = FALSE),
    get(name, envir = right, inherits = FALSE),
    tolerance = tolerance,
    check.environment = FALSE
  )
  if (!isTRUE(result)) {
    failures <- c(failures, paste0(name, ": ", paste(result, collapse = "; ")))
  }
}

if (length(failures)) {
  cat(paste(failures, collapse = "\n"), "\n")
  quit(status = 1L)
}

cat(
  sprintf(
    "Equivalent: %d objects agree within tolerance %.3g\n",
    length(left_names), tolerance
  )
)
