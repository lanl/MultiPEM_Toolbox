#!/usr/bin/env Rscript

# Maintainer utility: run this only when intentionally refreshing the runtime
# dependency set.  Normal image builds consume renv.lock and never update it.
options(repos = c(CRAN = "http://cran.rstudio.com"))

runtime_packages <- c(
  "Matrix",
  "numDeriv",
  "doFuture",
  "future",
  "iterators",
  "adaptMCMC",
  "FME",
  "ramcmc",
  "Rcpp",
  "RcppEigen"
)

lock_packages <- unique(c("renv", runtime_packages))
missing_packages <- lock_packages[
  !vapply(lock_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages)) {
  install.packages(
    missing_packages,
    dependencies = c("Depends", "Imports", "LinkingTo")
  )
}

unavailable <- runtime_packages[
  !vapply(runtime_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(unavailable)) {
  stop("Required packages unavailable: ", paste(unavailable, collapse = ", "))
}

project <- tempfile("multipem-lock-")
dir.create(project)
renv::snapshot(
  project = project,
  lockfile = "/out/renv.lock",
  packages = runtime_packages,
  prompt = FALSE,
  force = TRUE
)

cat("Wrote /out/renv.lock\n")
