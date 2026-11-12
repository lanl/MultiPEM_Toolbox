# Note: observed-information conditioning and MCMC seeding

**Status: informational only — no code change recommended or made.**

## Finding

On the IYDT 4-Phen `I-EIV-SUGAR-hob` posterior (d = 104), the inverse observed
information at the MAP (`p_cal$IHess`, from `observed_information_cal.r`) is very
ill-conditioned: condition number ~3e13, with ~14 eigenvalues below 1e-10 of the
largest. The matrix is, however, **positive-definite** (0 negative eigenvalues),
so the Cholesky-with-jitter loop in `observed_information_cal.r` /
`observed_information_0.r` is working correctly — it is not producing a broken
or indefinite result.

The ill-conditioning is a **genuine property of the posterior, not a numerical
defect**: several model-block directions (notably the emplacement-condition
forward-model coefficients around coords 47-72, posterior variances ~1e5) are
nearly unconstrained by the calibration data. The FGSN block is well-behaved
(variances ~1e-3 to ~1). No regularization in the observed-information code can
manufacture information the data do not contain.

## Why no change was made

1. **The observed-information files compute the right thing.** Their job is
   uncertainty quantification (confidence intervals, `II_nev`, `II_calp`). A
   near-singular `IHess` is the correct, honest answer there — it reports that
   some directions are poorly determined. Forcing a better condition number
   would falsify those uncertainty estimates.

2. **The seeding fragility is real but pre-existing and non-blocking.**
   `calc_bayes_cal.r` seeds RAM/FME with `lambda * IHess` as the proposal
   covariance (`lambda = 2.38^2 / npars`). A condition-number-1e13 proposal
   covariance is inefficient along the flat directions, but RAM/FME adapt around
   it, and the toolbox has always seeded this way. Nothing we changed affects
   this path.

3. **NUTS does not use `IHess`.** It adapts a dense metric during warmup
   (`nuts_metric="dense"`), so this conditioning issue is irrelevant to the NUTS
   path. (It is also why `IHess` was a poor *fixed* seed in a seeded probe.)

## If this is ever revisited

The correct fix would be an **optional** eigenvalue floor applied **only at the
MCMC seeding site** (where `lambda * IHess` becomes a proposal covariance in
`calc_bayes_cal.r`), leaving the uncertainty-quantification path
(`II_nev`, `II_calp`, confidence intervals) untouched. This separates the two
concerns: honest uncertainties from `obs_info_*`, a well-conditioned seed for
RAM/FME. It is a nice-to-have hardening, not a bug fix, and would require its own
verification.
