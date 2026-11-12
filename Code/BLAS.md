# Using an optimized BLAS/LAPACK with MultiPEM

MultiPEM's MLE optimization and NUTS sampling are dominated by dense linear
algebra (Cholesky factorizations, triangular solves, matrix products) in the
likelihood, gradient, and mass-matrix code. R ships with an unoptimized
*reference* BLAS by default. Linking R against a tuned BLAS can speed up these
operations substantially at no change to the MultiPEM code.

> **Scope.** The gain depends on the problem. Dense operations (the MLE stage,
> and the dense NUTS metric) benefit most. For a target whose per-evaluation
> cost is dominated by interpreted-R work rather than large matrix operations,
> the BLAS speedup can be modest — measure before and after on your own case.

## Verify which BLAS R is using

```r
sessionInfo()$BLAS        # path to the active BLAS shared library
extSoftVersion()[["BLAS"]]
```

If the path contains `libRblas` (reference) you are on the unoptimized library.

A quick before/after microbenchmark:

```r
set.seed(1); n <- 1000; A <- crossprod(matrix(rnorm(n*n), n))
system.time(for (i in 1:20) chol(A))   # compare under each BLAS
```

**Always re-run the test suite after switching BLAS.** An optimized BLAS can
differ from the reference in the last few floating-point digits, which is
expected; confirm results still match to tolerance
(your application tests).

---

## macOS — Apple Accelerate (vecLib)

The CRAN R build for macOS already ships an Accelerate-backed BLAS alongside the
reference one; switching is a single symlink, no build or download required.

```sh
cd /Library/Frameworks/R.framework/Resources/lib
ls libRblas*                       # expect libRblas.0.dylib and libRblas.vecLib.dylib
sudo ln -sf libRblas.vecLib.dylib libRblas.dylib   # use Accelerate
```

Revert at any time:

```sh
cd /Library/Frameworks/R.framework/Resources/lib
sudo ln -sf libRblas.0.dylib libRblas.dylib        # back to reference
```

Restart R, then confirm with `sessionInfo()$BLAS` (the path should now reference
`vecLib`/`Accelerate`). The change only affects newly started R sessions; it
does not disturb a job already running.

> If your R was not installed from the CRAN package and lacks
> `libRblas.vecLib.dylib`, build OpenBLAS (below) or reinstall the CRAN R.

---

## Windows — replace Rblas.dll

Windows has no Accelerate equivalent. The standard route is to replace the
reference `Rblas.dll` with an optimized drop-in build (e.g. an OpenBLAS build
packaged as `Rblas.dll`).

1. Locate `Rblas.dll` under your R installation, typically
   `C:\Program Files\R\R-x.y.z\bin\x64\Rblas.dll`.
2. Back it up: copy to `Rblas.dll.reference`.
3. Obtain an optimized `Rblas.dll` matching your R major version and
   architecture (x64) — for example an OpenBLAS build distributed for R — and
   copy it over the original.
4. Restart R and verify with `sessionInfo()`.

To revert, restore `Rblas.dll` from your backup.

> Match the DLL to your R version/architecture; a mismatched build will fail to
> load or crash R. Keep the backup so you can always roll back.

---

## OpenBLAS (any platform; the portable option)

OpenBLAS is a tuned, open-source BLAS/LAPACK that works across macOS, Linux, and
Windows and is a good choice when you want the same optimized library on every
machine.

**Install the library**

- macOS: `brew install openblas`, or build from the official OpenBLAS source
  release.
- Linux (Debian/Ubuntu): `sudo apt-get install libopenblas-dev`.
- Linux (RHEL/Fedora): `sudo dnf install openblas-devel`.
- Any platform: build from the OpenBLAS source release.

**Point R at it**

- **Linux**, the clean system-wide way, is the alternatives mechanism:

  ```sh
  sudo update-alternatives --config libblas.so.3     # pick the OpenBLAS entry
  sudo update-alternatives --config liblapack.so.3
  ```

  Or, without changing system defaults, preload for one session:

  ```sh
  LD_PRELOAD=/path/to/libopenblas.so Rscript your_run.r
  ```

- **macOS**: replace the symlink target as in the Accelerate section, pointing
  `libRblas.dylib` at the OpenBLAS dylib (e.g. the `brew` path under
  `$(brew --prefix openblas)/lib`), or build R with
  `--with-blas="-lopenblas" --with-lapack`.

- **Windows**: supply the OpenBLAS build as `Rblas.dll` per the Windows section.

**Threading note.** OpenBLAS is multithreaded. When MultiPEM already runs many
parallel MCMC chains / optimization starts (one per physical core), let BLAS use
a single thread to avoid oversubscription:

```sh
export OPENBLAS_NUM_THREADS=1       # also OMP_NUM_THREADS=1
```

Reserve multithreaded BLAS for runs that use *fewer* chains than cores, so the
freed cores can accelerate individual matrix operations.

---

## Which should I pick?

- **macOS, quickest win:** Accelerate (one symlink, already present).
- **Reproducibility across mixed machines:** OpenBLAS everywhere.
- **Windows:** optimized `Rblas.dll` (OpenBLAS build).

Whichever you choose, re-run the tests and benchmark your own case — the payoff
is problem-dependent.
