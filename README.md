# Multi-Phenomenology Explosion Monitoring (MultiPEM) Toolbox

MultiPEM combines observed signatures from one or more phenomenologies to
estimate characteristics of a new explosive event with quantified uncertainty.
The supplied IYDT application models seismic, acoustic, optical, and
surface-effects signatures using empirical or physics-based forward models
together with source, path, and observational-error terms.

The repository supports two assessment modes:

- A **rapid assessment** first calibrates forward and error models with benchmark
  data (`runMPEM.r`), then analyzes new-event data (`runMPEM_0.r`). A saved
  calibration can be reused for multiple new events.
- A **complete assessment** uses benchmark and new-event data in one joint
  calculation (`runMPEM.r`).

Analyses can run natively in R or through the versioned Docker workflow. Docker
is the recommended route when reproducibility across machines matters.

## Documentation and references

- [`multipem-um-la-ur-26-24464.pdf`](multipem-um-la-ur-26-24464.pdf) is the user
  manual. It documents the public calculation functions and walks through the
  IYDT-GSRP rapid and complete analyses.
- [`multipem-la-ur-23-21950.pdf`](multipem-la-ur-23-21950.pdf) is the detailed
  methodology report (LA-UR-23-21950, revision 4).
- [`ggag125.pdf`](ggag125.pdf) is the associated *Geophysical Journal
  International* article.
- [`MultiPEM-GSRP-030425.pdf`](MultiPEM-GSRP-030425.pdf) discusses crossed path
  effects and related extensions represented in the supplied run files.
- [`Runfiles-Docker/README`](Runfiles-Docker/README) is the operational reference
  for versioned analysis jobs.
- [`Test-Docker/README`](Test-Docker/README) is the operational reference for
  isolated verification-test jobs.

## Quick start with Docker

The host needs R 4.0 or newer with `Rscript`, a Docker client connected to a
running Linux-container engine, and permission to bind-mount the repository and
results locations. Host R packages are not required.

From `Runfiles-Docker`:

```text
Rscript mpem.R doctor
Rscript mpem.R build
Rscript mpem.R list
Rscript mpem.R run IYDT-gsrp/Seismic/I-SUGAR-hob
```

The last command creates an auditable private snapshot and starts a detached
container. It prints the job name and results directory. Use the job name with:

```text
Rscript mpem.R status JOB
Rscript mpem.R logs -f JOB
Rscript mpem.R wait JOB
Rscript mpem.R clean JOB
```

Read [`Runfiles-Docker/README`](Runfiles-Docker/README) before staging a
multi-phenomenology job: those decks can require `opt.RData` files produced by
earlier single-phenomenology jobs.

Run the verification suites independently from `Test-Docker`:

```text
Rscript testmpem.R doctor
Rscript testmpem.R build
Rscript testmpem.R run global
Rscript testmpem.R run IYDT/Seismic
```

## Native R requirements

The repository does not impose a native R version or package lock. R 4.4.3 is
recommended when a native run is intended to resemble the pinned Docker
runtime. Install the packages required by the selected algorithms:

- `Matrix`, `numDeriv`, `doFuture`, and `future` for shared calculations;
- `adaptMCMC` for RAM sampling;
- `FME` for DRAM/FME sampling;
- `Rcpp` and `RcppEigen`, plus a working C/C++ toolchain, for NUTS; and
- `ramcmc` and `iterators` for SMC.

The supplied IYDT application uses this shared package set. A future
application may impose additional packages through its forward-model code;
document and lock those dependencies before relying on them.

For example:

```r
install.packages(c(
  "Matrix", "numDeriv", "doFuture", "future", "adaptMCMC", "FME"
))
```

## Native repository setup

The canonical run files use relative paths. From `Runfiles`, create the shared
code link once:

```text
ln -s ../Code Code
```

Each application directory also needs `Data` and `Code` links. For example:

```text
cd IYDT-gsrp
ln -s ../../Applications/Data/IYDT-gsrp Data
ln -s ../../Applications/Code/IYDT-gsrp Code
```

Some phenomenology directories require an additional `Code` link for prior
functions. Their local README states when and where to create it. These links
are native-run conveniences; Docker prepares equivalent directories inside
each job snapshot and does not require canonical links.

From an analysis directory, run a complete deck or a rapid calibration with:

```text
R CMD BATCH runMPEM.r runMPEM.out &
```

The trailing `&` returns the shell prompt while R runs. For a rapid assessment,
monitor `runMPEM.out` and wait for successful calibration. Before starting any
event deck, preserve the resulting workspace:

```text
cp .RData .RData-calibration
```

Then start the event stage in the background:

```text
R CMD BATCH runMPEM_0.r runMPEM_0.out &
```

The event deck reads and then modifies `.RData`. Keep `.RData-calibration` as a
pristine benchmark/calibration result for all future event analyses, copying it
back to `.RData` or into a new event directory as needed. See
[`Runfiles/README`](Runfiles/README) and the application-specific README before
running a deck.

## Repository layout

| Path | Purpose |
| --- | --- |
| `Applications/Code/` | Application forward models, Jacobians, transforms, priors, and result formatting |
| `Applications/Data/` | Canonical CSV inputs used by application decks |
| `Applications/Test/` | Application-specific verification fixtures |
| `Code/` | Application-independent preprocessing, likelihood, posterior, optimization, and sampler code |
| `Runfiles/` | Canonical native run decks grouped by application and analysis |
| `Runfiles-Docker/` | Versioned Docker analysis runner, dependency lock, and job management |
| `Test/` | Canonical global and IYDT verification-suite entry points |
| `Test-Docker/` | Versioned Docker verification runner and job management |

The supplied run-file group is:

- `IYDT-gsrp`: the four-phenomenology example used by the manual and reports.

`Runfiles-Docker/applications.tsv` is the authoritative mapping from each
run-file group to its application code and data directories.

## Run decks and output

`runMPEM.r` performs either rapid calibration or a complete joint analysis.
`runMPEM_0.r`, when present, performs rapid new-event inference using the saved
calibration state. Common durable native outputs are:

- `.RData`, the saved calculation workspace;
- `runMPEM.out` and `runMPEM_0.out`, the R batch transcripts; and
- `opt.RData` or `opt_nev.RData`, saved optimization results when enabled by the
  deck.

Multi-phenomenology calibration decks may hard-code single-phenomenology
optimization files under a sibling `Opt` directory. Create that directory and
copy the exact filenames listed in the local README before starting the deck.

The supplied decks request multicore execution and often enable Bayesian
analysis. Review `ncores_*`, `iBayes`, `iMCMC`, sample counts, and random seeds
before a production run. Native `future` backend availability is
platform-dependent. The Docker runner caps effective worker pools to the CPU
capacity visible to the container while preserving the decks' logical counts.

## Verification and reproducibility

Run all supplied isolated verification suites as described in
[`Test-Docker/README`](Test-Docker/README). The analysis runner also provides
runtime comparison commands:

```text
Rscript mpem.R verify
Rscript mpem.R compare LEFT.RData RIGHT.RData 1e-12
```

The Docker runtime pins its base image and R package versions in
`Runfiles-Docker/Dockerfile.runtime` and `Runfiles-Docker/renv.lock`. Each job
records input hashes and writes only to its private results snapshot.

## Extending the toolbox

An application normally contributes:

1. code under `Applications/Code/NAME`;
2. immutable input data under `Applications/Data/NAME`;
3. one or more decks under `Runfiles/NAME`; and
4. a row in `Runfiles-Docker/applications.tsv`.

Add application tests under `Applications/Test` and register their group in
`Test-Docker/test-groups.tsv` when applicable. Follow the interface and
registration checklists in [`Applications/README`](Applications/README),
[`Runfiles/README`](Runfiles/README), and the two Docker workflow READMEs.

## License

This toolbox is open source under the BSD 3-Clause License. Its LANL-internal
software identifier is O4673.
