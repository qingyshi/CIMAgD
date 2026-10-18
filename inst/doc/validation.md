# Validation of CIMAgD 0.1.0

Date: 2026-10-18. Platform: R 4.5.2, macOS arm64. This report validates the 2026-10-18 update in [version.md](version.md), relative to the 2026-10-16 development snapshot. The package version remains 0.1.0.

## Results

`devtools::test()` passes. The installed-package test run reports **40 test blocks / 250 passing expectations**, with zero failures, warnings or skips. `R CMD build .` succeeds, and `R CMD check --no-manual CIMAgD_0.1.0.tar.gz` completes with **Status: OK**. The check covers installation, namespace loading, R code, documentation, examples, data and tests. Remote repository indexes were unavailable in this restricted environment; installed dependencies sufficed. PDF-manual compilation is excluded by `--no-manual`.

## Checks against the version history

| Version entry | Validation and result |
|---|---|
| 1. Optional Jacobian diagnostics | Defaults resolve to none with diagnostics enabled and minimum-norm otherwise. Disabling diagnostics returns NULL reports, preserves full-rank coefficients/covariance, and still respects explicit none. |
| 1. Moore-Penrose solutions | Redundant CATE coefficients agree with an independent raw-coordinate SVD solution. Redundant H_s satisfies Moore-Penrose identities. Constant and duplicated tilting functions preserve calibration and fitted effects. Conflicting redundant targets fail under both regularization methods. |
| 1. Scaled ridge | GMM matches an independently constructed normalized penalized least-squares solution. Rescaling duplicated predictors by 10^8 preserves fitted effects. H_s ridge inverses transform correctly under diagonal scales of 10^6 and 10^-3. |
| 1. Inference after singularity | A regularized singular fit returns point estimates, NA ordinary covariance/intervals, an explanatory warning and separately labeled finite numerical sensitivity covariance. |
| 2. Optional tilting-span diagnostic | Skipped by default; can be enabled independently of Jacobian diagnostics. Known sufficient and insufficient spans, redundant diagnostic bases, residual Gram matrices and the non-rejection interpretation remain tested. |
| 3. Copula inputs and margins | Named lists, normalized tables and native versus moment parameters generate equivalent draws. Normal, Bernoulli, Gamma, Beta, lognormal and exponential marginal moments/support are checked. Invalid Beta/Gamma parameters fail. |
| 3. Copula uncertainty | Plugin uses the original n. Bootstrap refits marginal parameters and latent dependence. Mixed continuous/Bernoulli and Bernoulli/Bernoulli dependence recovery is checked against latent correlations. Gamma covariates run through both downstream uncertainty methods. Independent target/trial variance terms remain, without double counting plugin Q variance. |
| 4. Two-step GMM | Second-stage weights agree with the inverse pilot moment covariance in a full-rank example. Coefficients satisfy independently computed weighted normal equations. Final covariance agrees with the complete sandwich evaluated at the final coefficients. |
| 4. Singular moment covariance | Exact redundant covariance directions are removed without ridge and retain identifiable coefficients. Informative zero-variance directions fail unless ridge is explicit. Two-step weighting also works with final variance reporting disabled. |
| 4. Two-step copula bootstrap | Pilot weights use refitted-copula moment covariance plus trial covariance. Subsequent fits hold the estimated weight fixed; bootstrap Q covariance replaces the analytic integration-sample term. All 12 seeded test replicates succeed. |
| 5. Computation | Full-size before/after coefficient and covariance comparisons agree to numerical precision. Profiling and returned-object memory are reported below. |
| 6. Documentation | Roxygen help regenerated; package check finds no missing help, usage mismatches or invalid examples. The version record describes implemented changes. Package source/help no longer identify the paper by a local manuscript filename. |
| 7. Data independence | For each of dataHF and sumstatsPULSE, the RDA object is identical to both development and installed RDS copies. Examples/tests pass after removal of the raw CSV and its creation script. Runtime source, examples and tests contain no references to those removed files. |

The focused tests are in `test-regularization.R`, `test-two-step.R` and `test-copula.R`. The existing tests continue to cover input completion, exact linear coefficient recovery, normalized-ratio derivatives, joint trial covariance reconstruction, Q/trial/target cross terms, time models and the HF/PULSE example.

In particular, the A.6.2 covariance test independently enumerates every (X,A,Y) state in a finite binary model. Other tests cover contrast/summary cross-covariances, reported-variance anchoring, fixed normalization variance, zero-event reports, bounded working control risks, complete supplied Sigma and cross-time covariance. These checks protect the variance estimator while the solver and GMM weight change.

The simulation regression test retains 60 seeded binary replicates. Its tolerances are maximum coefficient bias below .02, target-effect bias below .012, empirical-to-estimated target variance ratio between .5 and 1.8, and pointwise coverage at least .83. This is a regression test, not a large simulation study establishing universal coverage or efficiency gains from two-step GMM.

## Full-size example and computational assessment

The additional benchmark uses 50,000 PULSE copula draws, four trials, the three available effect modifiers, target=Q, and seed 20261016. Detailed Jacobian diagnostics remain enabled. The previous snapshot computed the span diagnostic by default; the updated default skips it.

| Quantity | Previous snapshot | Updated default |
|---|---:|---:|
| One-step target risk difference | -0.03373277 | -0.03373277 |
| One-step target standard error | 0.00556189 | 0.00556189 |
| Returned object size, decimal MB | 132.33 | 67.53 |
| Elapsed time with profiling, seconds | 1.211 | 1.105 |

The maximum absolute coefficient difference is 7.63e-17; the maximum covariance-entry difference is 1.14e-18. Returned-object memory decreases by about 49%. Timings are single profiled runs on the same machine, not a reliable general speedup estimate. They support modest caching and avoiding an unrequested diagnostic, without changing convergence tolerances or statistical calculations.

The full-size two-step example also completes: moment Jacobian rank 4/4, moment covariance rank 26, ordinary covariance inverse with no ridge, target risk difference -0.03336176, SE 0.00533732 and pointwise 95% interval [-0.04382272, -0.02290080]. Different one-step and two-step estimates are expected when the reported moments are not exactly matched by the linear model. The paper's full-rank inference still relies on its valid-moment assumptions; this example does not establish those assumptions clinically.

Time=1 and the illustrative PULSE n=50000 remain documented example conventions. All 42 HF reports and PULSE numerical summaries are retained.

## Scope and limitations

| Paper component | Implementation |
|---|---|
| Main normalized tilting and linear-basis GMM | R/tilting.R, R/gmm.R, R/CIMA.R |
| A.1 subgroup and A.3 time-indexed target effects | R/target.R and basis construction in R/gmm.R |
| C2 / C3 Jacobian rank and conditioning | R/diagnostics.R, R/numerics.R |
| A.4.1 sufficient-span projection under Q | diagnose_tilting() |
| A.6.1 combined Q and joint trial influences | R/vcov.R |
| A.6.2 bounded working control regression and covariance anchoring | R/control.R |
| Fixed-copula plugin and refitted parametric bootstrap | R/copula.R, R/bootstrap.R |

Regularization selects a numerical solution without establishing identification; ordinary Wald inference is withheld for singular systems. The sensitivity covariance is not a substitute. Bootstrap uncertainty assumes the Gaussian copula model and the documented marginal/dependence refitting rules; it cannot recover the sampling behavior of an unknown original summary estimator. Fixed second-stage bootstrap weights use first-order equivalence under valid population GMM moments, not misspecification-robust inference.

Arm-report covariance reconstruction is a working approximation for complete-outcome randomized-trial reports, not exact inference for arbitrary adjusted or censored effects. Continuous reconstruction assumes constant conditional residual variance within each arm. Repeated horizons require complete joint Sigma; partial Q/target overlap is outside the exact-reuse/independence cases. Indirect comparisons, relative effects, A.5 misspecification inference and A.7 arm-mean transport estimation remain excluded.

The paper was not edited. Its SHA-256 remains:

```
4783b732165e2a35e5b84ccb83110bd0afed200a6811ac794d6233e749cb50a5
```
