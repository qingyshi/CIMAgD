# Validation of CIMAgD 0.1.1

Date: 2026-11-16. Platform: R 4.5.2, macOS arm64. This report validates version 0.1.1 in [version.md](version.md) against the November 2026 version of the paper and the committed 0.1.0 code.

## Results

`devtools::test()` passes. The installed-package test run reports **54 test blocks / 377 passing expectations**, with zero failures, warnings or skips. `R CMD build` succeeds, and `R CMD check --no-manual CIMAgD_0.1.1.tar.gz` completes with **Status: OK**. The check covers installation, namespace loading, R code, documentation, examples, data and tests. Remote repository indexes were not queried; installed dependencies sufficed. PDF-manual compilation is excluded by `--no-manual`.

## Checks against the version history

| Version entry | Validation and result |
|---|---|
| 1. Theory check | Every 0.1.0 test passes unchanged: exact linear coefficient recovery, numerical normalized-ratio and summary derivatives, the A.6.3 covariance oracle that enumerates every (X,A,Y) state, Q/trial/target cross terms, two-step sandwich covariance, regularization and copula propagation. On the full-size HF example, the 0.1.0 and 0.1.1 common-method point estimates agree exactly; the covariance changes only through the working control models of section 3. |
| 2. MVLP equals separate fits | For each of two trials, MVLP coefficients, transported overall and subgroup effects, and their covariance blocks agree with a separate single-trial fit to 1e-10. Trial influence terms are block diagonal across trials; the shared Q terms couple them. |
| 2. A.7 covariance and pooling | The trial-specific covariance equals the A.7 expression assembled independently from the target, Q and trial influence functions when Q reuses the target, and the sum of the target, Q and trial components when they are independent. For each estimand, weights equal `V^{-1} 1 / (1' V^{-1} 1)` and sum to one, the pooled estimate equals the weighted trial-specific estimates, the pooled variance equals `1 / (1' V^{-1} 1)`, and the joint pooled covariance equals `L' V L`. |
| 2. Special cases | Single-trial MVLP reproduces the common fit with unit weights. Two-step MVLP weights invert each trial's own moment-covariance block, are zero across trials, and solve each trial's weighted normal equations. Point-only MVLP returns the same pooled estimates and weights with NA covariance. Continuous and time-indexed fits (with complete Sigma) pool each horizon. Copula bootstrap MVLP preserves the component and pooling identities. |
| 2. Per-trial identification | When one trial reports subgroups for a single covariate, MVLP stops with a `cima_singularity` error naming that trial and its rank, while the common fit remains identified. Per-trial ranks are reported. Minimum-norm regularizes only the singular trial and withholds Wald inference. On the HF data, MVLP with LVEF in the CATE stops at DAPA-HF (rank 3 of 4). |
| 3. Automatic working control | Stable trial-specific fits are retained. A trial whose own working mean cannot identify a term switches to the common fit, which equals an independent inverse-variance weighted least-squares oracle over both trials' control reports; the other trial's working fit equals its stand-alone fit. Active risk bounds switch to an interior common fit; a failing common fit leaves the trial-specific fit in place; infeasible risks stop with an explicit message; a single trial never attempts a common fit. With MVLP, the borrowed working mean satisfies the trial's own CATE bounds and leaves its CATE estimate unchanged. |
| 4. README correction | The HF test asserts that the binary summaries in `mu` equal the `dataHF` subgroup proportions to three decimals. |
| Documentation | Roxygen help regenerated; package check finds no missing help, usage mismatches or invalid examples. The new MVLP example runs. |

The focused tests are in `test-mvlp.R` and `test-control.R`. The other tests continue to cover input completion, exact linear coefficient recovery, normalized-ratio derivatives, joint trial covariance reconstruction, Q/trial/target cross terms, time models, regularization, two-step GMM, copula propagation and the HF/PULSE example.

The simulation regression test retains 60 seeded binary replicates. Its tolerances are maximum coefficient bias below .02, target-effect bias below .012, empirical-to-estimated target variance ratio between .5 and 1.8, and pointwise coverage at least .83. This is a regression test, not a large simulation study establishing coverage or efficiency.

## Full-size HF example

The benchmark uses 50,000 PULSE copula draws, four trials, the three available effect modifiers, target=Q, seed 20261016 and the corrected summaries. The 0.1.0 column was computed from the committed 0.1.0 code with the same corrected summaries.

| Common CATE, `~ LVEF + preHHF + diabetes` | 0.1.0 code | 0.1.1 |
|---|---:|---:|
| One-step target risk difference | -0.03438889 | -0.03438889 |
| One-step standard error | 0.00543859 | 0.00542982 |
| One-step 95% interval | [-0.04504834, -0.02372944] | [-0.04503114, -0.02374665] |
| Two-step target risk difference | -0.03515536 | -0.03489483 |
| Two-step standard error | 0.00528687 | 0.00533477 |

The one-step point estimate does not depend on the covariance. The standard error changes because DAPA-HF and EMPEROR-Reduced report no LVEF subgroups: their trial-specific working means dropped the diabetes term and needed active risk bounds, so 0.1.1 uses the common working mean, which is full rank and interior, for these two trials only. EMPEROR-Preserved and DELIVER keep their trial-specific working means. Two-step weights depend on the moment covariance, so the two-step estimate also changes. The moment Jacobian has rank 4/4 and the two-step moment covariance has rank 26, inverted without ridge. Both one-step intervals round to the paper's Table 3 (-0.034; -0.045, -0.024). With the README's previous summaries (.405 and .474), 0.1.0 gave -0.03373277 with standard error 0.00556189, whose upper limit -0.0228 does not round to Table 3.

With MVLP and `~ preHHF + diabetes`, every trial has rank 3/3 and keeps its trial-specific working mean. The trial-specific overall effects (standard error; weight) are EMPEROR-Preserved -0.0292 (0.0095; 0.325), DELIVER -0.0279 (0.0095; 0.323), DAPA-HF -0.0389 (0.0117; 0.214) and EMPEROR-Reduced -0.0510 (0.0146; 0.138). The pooled effect is -0.03385666 with standard error 0.00541168 and 95% interval [-0.04446337, -0.02324996]; two-step MVLP gives -0.03419119 (standard error 0.00538811). The common CATE with the same basis gives -0.03426961 (standard error 0.00542357). Cross-trial covariances are of order 1e-8 because the 50,000 shared draws contribute little. A2p is doubtful here because the HFrEF and HFpEF trials enroll different LVEF ranges, so this is an interface check, not a clinical analysis.

Single runs took about 1.2 seconds for the common fit and 0.8 seconds for MVLP; the returned objects are about 67 and 72 decimal MB. Time=1 and the illustrative PULSE n=50000 remain documented example conventions. All 42 HF reports and PULSE numerical summaries are retained.

## Scope and limitations

| Paper component | Implementation |
|---|---|
| Main normalized tilting and linear-basis GMM | R/tilting.R, R/gmm.R, R/CIMA.R |
| A.1 subgroup and A.3 time-indexed target effects | R/target.R and basis construction in R/gmm.R |
| C2 / C3 Jacobian rank and conditioning, per trial for MVLP | R/diagnostics.R, R/numerics.R, R/gmm.R |
| A.4.1 sufficient-span projection under Q | diagnose_tilting() |
| A.6.1 combined Q and joint trial influences | R/vcov.R |
| A.6.3 trial-specific or common working control regression and covariance anchoring | R/control.R |
| A.7 trial-specific CATEs and minimum-variance linear pooling | R/gmm.R, R/target.R, R/numerics.R |
| Fixed-copula plugin and refitted parametric bootstrap | R/copula.R, R/bootstrap.R |

Regularization selects a numerical solution without establishing identification; ordinary Wald inference is withheld for singular systems. The sensitivity covariance is not a substitute. Bootstrap uncertainty assumes the Gaussian copula model and the documented marginal/dependence refitting rules; it cannot recover the sampling behavior of an unknown original summary estimator. Fixed second-stage bootstrap weights use first-order equivalence under valid population GMM moments, not misspecification-robust inference.

Arm-report covariance reconstruction is a working approximation for complete-outcome randomized-trial reports, not exact inference for arbitrary adjusted or censored effects. The common working control mean borrows control-arm information across trials as a working model only; it does not change either CATE estimator, and the variance remains an approximation. Continuous reconstruction assumes constant conditional residual variance within each arm. Repeated horizons require complete joint Sigma; partial Q/target overlap is outside the exact-reuse/independence cases.

MVLP relies on A1p and A2p, which aggregate data cannot verify, and on identification within every trial. Pooled subgroup effects additionally require the subgroup version of A1p. Minimum-variance weights can be negative when trial-specific effects are strongly correlated.

Indirect treatment comparisons (Section 2.4), relative effects (A.2), the misspecification-robust variance (A.6.2) and the conditional arm-mean estimator (A.8) are not implemented; A.8-type working control means are used only for covariance reconstruction.
