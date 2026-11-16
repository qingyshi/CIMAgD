# CIMAgD

CIMAgD estimates an additive conditional average treatment effect (CATE) from aggregate randomized-trial reports, then averages it in a target population. Version 0.1.1 implements the paper’s method: normalized exponential tilting, linear-basis GMM, target subgroups and optional time models, identification/span diagnostics, and joint influence-function variance estimation. Two estimation methods are available: one CATE common to all trials (`method = "common"`, the default) or trial-specific CATEs whose transported effects are combined by minimum-variance linear pooling (`method = "MVLP"`).

## Quick example

```r
library(CIMAgD)
data(dataHF)
data(sumstatsPULSE)

# PULSE reports these three covariates, but not NYHA or BMI.
reports <- subset(dataHF, !subgroup %in% c("NYHA", "BMI"))
mu <- list(
  "EMPEROR-Preserved" = c(54.3, 54.3^2 + 8.8^2, .229, .491),
  "DELIVER" = c(54.3, 54.3^2 + 8.9^2, .261, .448),
  "DAPA-HF" = c(31.1, 31.1^2 + 6.8^2, .274, .451),
  "EMPEROR-Reduced" = c(27.2, 27.2^2 + 6.1^2, .309, .498)
)
fit <- CIMA(
  reports, mu,
  summaryX = sumstatsPULSE,
  formula_cate = ~ LVEF + preHHF + diabetes,
  formula_tilting = ~ LVEF + I(LVEF^2) + preHHF + diabetes,
  target_subgroups = list(
    diabetes = function(d) d$diabetes == 1,
    joint = function(d) d$LVEF < 40 & d$diabetes == 1
  ),
  seed = 20261016
)
fit
coef(fit)
vcov(fit)
fit$target$estimates
vcov(fit, type = "target")
```

This example uses a fitted copula as Q and reuses its draws for the target. `sumstatsPULSE$n = 50000` is the original example's illustrative sample size, not the clinical cohort size. `dataHF$Time = 1` is a default coding convention; it does not establish clinical comparability of the original outcome horizons. All 42 original rows are retained in `dataHF`. The `mu` vectors give each trial's LVEF mean, LVEF raw second moment, and the proportions with a history of heart-failure hospitalization and with diabetes; the two proportions equal the subgroup proportions implied by `dataHF` to three decimals.

## Common CATE or pooled trial-specific CATEs

`method = "common"` (default) assumes one CATE function transportable across all trials and the target population (A1) and estimates it from the stacked moments of all trials. The trials then jointly identify the CATE: a trial without subgroup reports for a covariate still contributes through its other moments.

`method = "MVLP"` implements Appendix A.7. It estimates a separate CATE for each trial from that trial's reports only, marginalizes every trial-specific CATE over the same target sample, and combines the trial-specific transported effects by minimum-variance linear pooling. It replaces A1-A2 by A1p (every trial-specific CATE, marginalized over the target, gives the same ATE) and A2p (every trial covers the target support), and needs identification within every trial. The input format is the same for both methods.

```r
# DAPA-HF and EMPEROR-Reduced report no LVEF subgroups, so a trial-specific
# CATE in LVEF is not identified for them:
# CIMA(..., formula_cate = ~ LVEF + preHHF + diabetes, method = "MVLP")
# stops with "J_m_theta for trial DAPA-HF is singular: rank 3 of 4".
pooled <- CIMA(
  reports, mu,
  summaryX = sumstatsPULSE,
  formula_cate = ~ preHHF + diabetes,
  formula_tilting = ~ LVEF + I(LVEF^2) + preHHF + diabetes,
  method = "MVLP",
  target_subgroups = list(diabetes = function(d) d$diabetes == 1),
  seed = 20261016
)
pooled
pooled$target$trial_estimates      # trial-specific transported effects, with pooling weights
pooled$target$estimates            # pooled effects
pooled$pooling$weights             # one weight vector per target estimand
pooled$diagnostics$identification  # moment-Jacobian rank of every trial
coef(pooled)                       # one column of CATE coefficients per trial
vcov(pooled, type = "trial")       # joint covariance of the trial-specific effects
```

The trial-specific transported effects share the target sample (and the Q sample when it differs), so they are correlated. Their joint covariance V follows A.7, including these shared-sample terms, and the pooled effect uses `w = V^{-1} 1 / (1' V^{-1} 1)`, with variance `1 / (1' V^{-1} 1)`. Every target estimand, i.e., the overall effect and each target subgroup at each horizon, is pooled with its own weights, and the joint covariance of the pooled estimates is returned. Weights can be negative when trial-specific effects are strongly correlated. Pooled subgroup effects additionally require the subgroup version of A1p. A singular V stops the fit unless a singularity method is chosen, in which case Moore-Penrose weights are used and Wald inference is withheld. Since the weights require V, MVLP computes the covariance even with `variance = FALSE` but then does not report it. In the HF example, the HFrEF and HFpEF trials enroll different LVEF ranges, so A2p is doubtful; the example only illustrates the interface.

For MVLP, `coef()` returns a term-by-trial matrix, `predict()` returns one CATE column per trial, `fit$CATE$theta` and `vcov(fit)` use stacked `trial:term` coefficients, and `fit$ATE$by_trial` lists the trial-specific overall effects. A supplied `W` contributes only its trial-diagonal blocks, and two-step GMM inverts each trial's own moment covariance block.

## Aggregate reports

Each row is an overall or subgroup additive contrast for one trial, comparison and horizon. Arm 1 is treatment; arm 2 is control.

| Columns | Meaning / default |
|---|---|
| `trial`, `subgroup`, `stratum` | Required identifiers. `subgroup="overall"`, `stratum="TRUE"` denotes an overall report. |
| `event1`, `n1`, `event2`, `n2` | Required binary arm reports; event counts are integers. |
| `mean1`, `se1`, `n1`, `mean2`, `se2`, `n2` | Alternative continuous reports; SEs are standard errors of arm means, **not** outcome SDs. |
| `comparison` | Defaults to `"A vs B"`. Analyze one consistently oriented comparison per fit. |
| `TE`, `seTE` | Preserved when supplied; missing entries are calculated as arm 1 minus arm 2 and the independent-arm SE. Binary TE is a risk difference. |
| `n` | Defaults to `n1+n2`; a supplied value must agree. |
| `Time` | Defaults to 1. Used in the model only when `time_indexed=TRUE`. |

`prepare_agd()` validates and completes reports. Binary data need no mean/SE columns. Zero event counts require no continuity correction; zero reported variances remain zero. For the working control regression only, a zero arm SE uses the smallest positive arm SE as a finite weighting reference, or equal weights if all are zero; this is recorded in the fit.

A stratum can be a level label or a comparison such as `"x >= 40 & x <= 49"`, where `x` denotes the named subgroup variable. Boolean comparisons and `%in% c(...)` are supported; arbitrary executable expressions are not. Empty subgroups under Q fail explicitly.

Total trial size and allocation are inferred from an overall row or a disjoint exhaustive partition verified on Q. Otherwise supply `trial_info = data.frame(trial=..., n=..., allocation=...)`, where allocation is the treatment probability. Never infer a trial's total size from its largest subgroup alone.

## Tilting, input summaries and Q

`mu` is a trial-named list of point estimates corresponding to `formula_tilting`. Named vectors are matched to model-matrix column names. The normalization intercept is fixed at one and can be omitted from `mu`. A trial-named list of formulas supports different reported summaries across trials. Means, squared terms, interactions and subgroup indicators can all be basis functions.

Reported sample proportions do not replace the fitted normalizers `mean_Q(w_s * I_sj)`. Set `match_subgroup_proportions=TRUE` to add independent subgroup-proportion constraints explicitly. Redundant inconsistent constraints raise an error. In particular, the HF covariate summaries are rounded to three decimals, so they agree with the event-analysis denominators only up to rounding and the exact consistency check rejects them; the example keeps `match_subgroup_proportions=FALSE` and does not silently replace either source.

For empirical Q, pass a complete covariate data frame. `target=NULL` reuses Q and combines their influences. A separate target is treated as independent unless it is exactly identical to Q. `q_target_shared` can declare exact row-aligned reuse or independent sampling. Partial overlap is not automatically handled.

For a copula, supply one source using the constructor or an equivalent list:

```r
summaryX <- copula_summary(
  n = 2000,
  marginals = list(
    age = list(distribution = "normal", mean = 60, sd = 10),
    marker = list(distribution = "gamma", mean = 2, sd = 1),
    score = list(distribution = "beta", mean = .4, sd = .1),
    diabetes = list(distribution = "bernoulli", prob = .25)
  )
)
Q <- XCopula(summaryX, seed = 1)
```

Omitting `correlation` specifies independence. Otherwise supply a latent Gaussian correlation matrix, optionally named by covariate. Observed Pearson correlations for mixed or nonnormal margins generally differ from latent correlations. Existing `marginals` data frames with a `variable` column remain supported.

| Distribution | Permissible point summaries |
|---|---|
| `normal` | `mean`, `sd` |
| `bernoulli` | `prob` or `mean`, strictly between 0 and 1 |
| `gamma` | `mean`, `sd`; or `shape` and `rate`/`scale` |
| `beta` on (0,1) | `mean`, `sd`; or `shape1`, `shape2` |
| `lognormal` | Arithmetic `mean`, `sd`; or `meanlog`, `sdlog` |
| `exponential` | `mean` or `rate` |

Here `sd` describes the marginal distribution; it is not a sampling SE. Only point estimates are needed, with no sampling variance/covariance matrix.

* `q_method="plugin"` (default) holds these parameters fixed and draws exactly `summaryX$n` iid observations. Analytic Q variance treats these as an empirical sample. It excludes uncertainty in the supplied model parameters.
* `q_method="bootstrap"` treats the parameters as estimates. It simulates n observations, re-estimates marginal parameters and latent dependence, and refits tilting, CATE and target effects, assuming the Gaussian copula model is correct. Its covariance replaces the analytic Q component; trial and independent-target components remain analytic. `bootstrap_reps=200` by default. Failed replicates are recorded; more than 20% failures stop the fit.

Bootstrap marginal refits use sample moments for normal, Gamma and Beta; a sample proportion for Bernoulli; arithmetic moments for arithmetic lognormal inputs or log moments for log-scale inputs; and inverse sample mean for exponential rate. The uncertainty corresponds to these documented refitting rules. Marginal CDFs transform continuous covariates into Gaussian scores before dependence refitting. Mixed continuous/Bernoulli pairs use conditional probit likelihood, and Bernoulli pairs use bivariate threshold likelihood. A positive-definite projection of a refitted pairwise correlation matrix, when needed, is recorded in the returned summary.

Integration uses common Gaussian draws across bootstrap refits, with `bootstrap_mc_size=max(10000, 2*n)` by default. Increase it and check stability when support is sparse or tilts are concentrated.

## GMM estimation

`gmm_method="one-step"` is the default. Its default weight is diagonal inverse `n_s*seTE^2`, or identity if any contrast SE is zero (decided separately for each trial with MVLP); a supplied positive-definite `W` overrides this initial weight.

`gmm_method="two-step"` first fits the model, estimates the full moment sampling covariance, and uses its inverse as the second-stage weight. This covariance includes both the combined Q contribution and joint trial effect/summary contributions. Purely redundant covariance directions are removed by a range inverse. Informative zero-variance directions require one-step estimation or explicit `singularity_method="ridge"`. Inspect `fit$gmm` for the pilot estimate, covariance source, rank and inversion method. Covariance reconstruction is required for two-step weighting even with `variance=FALSE`.

With copula bootstrap, refitted copulas estimate the pilot Q moment covariance. The resulting second-stage weight is held fixed across bootstrap coefficient fits, which is first-order equivalent under valid population GMM moments and avoids a nested bootstrap. The final analytic covariance is recomputed at the second-stage coefficients before its Q component is replaced by bootstrap covariance.

## Variance and diagnostics

`fit$variance` exposes all Jacobians, the combined Q influence `xi_q`, trial influence covariances, and separate Q/trial components. A.6.3 reconstruction uses the tilted trial law, a working control regression and the CATE, then anchors reported contrast and available covariate-summary variances. Treatment-effect/summary cross-covariances are retained.

For binary outcomes, the control fit obeys `max(0,-g) <= gamma0 <= min(1,1-g)` on Q, where `g` is the CATE used by that trial (the common CATE, or the trial's own CATE with MVLP). The default basis is the CATE basis; `formula_control` can select a smaller basis. Unidentified nuisance terms are dropped in formula order while retaining an intercept, and the selection is reported. This regression serves variance estimation only; it does not implement the A.8 arm-mean estimator or transport control risks to the target. Continuous reconstruction assumes constant residual variance within each arm, estimated from reported mean SEs and sample sizes; inspect the returned residual diagnostics.

The working control mean is chosen automatically, trial by trial. Each trial first fits its own trial-specific working mean, which needs no assumption about control risks in other trials. That fit is kept unless it fails or is unstable: unidentified working terms, an ill-conditioned weighted design (scaled condition number above `control$condition_limit`), active risk bounds, or, for continuous outcomes, a zero residual variance. Only for such a trial, CIMA fits a common working mean to the reported control means of all trials, each aggregated under its own tilted law and weighted by inverse squared SEs, while keeping that trial's basis, risk bounds and residual-variance reports. The common fit replaces the trial-specific one when it is no less stable; otherwise the trial-specific fit is retained. This applies to both methods: with MVLP the CATEs stay strictly trial-specific, but the working mean may borrow information from other trials, because it only approximates the missing covariances and does not affect the consistency of either CATE estimator. The variance remains a working approximation in either case. `fit$variance$working_control` records the model, the detected issues and the reason for every trial; `print(fit)` lists trials that use the common mean. In the HF example, DAPA-HF and EMPEROR-Reduced report no LVEF subgroups, so their own working means cannot identify all terms and the common mean is used for them.

Supply `mu_se` for known marginal covariate-summary SEs. Missing marginal variances use the working-law estimates. For complete information, supply `Sigma`, a trial-named list in the order **that trial's TE rows, then its augmented tilting basis h+ including the intercept**. It is the joint **influence covariance**, equal to trial size times the sampling covariance. The intercept row/column must be zero. Full Sigma bypasses working-control reconstruction for that trial. Matrix names, when supplied, must match the returned coordinate names exactly.

`jacobian_diagnostics=TRUE` reports scaled rank, singular values, conditioning and null directions for `J_m_theta` and `H_s`. With FALSE, detailed reports and conditioning warnings are skipped; necessary solver rank checks remain active. `fit$diagnostics$identification` always summarizes the moment-Jacobian rank: one row for the common CATE, or one row per trial with MVLP, where every trial must reach full rank.

`singularity_method` controls singular systems:

* `"none"`: stop with a diagnostic-bearing `cima_singularity` error. This is the default with Jacobian diagnostics enabled.
* `"minimum-norm"`: use a Moore-Penrose solution in the original coefficient coordinates. This is the default when Jacobian diagnostics are disabled; it adds no penalty or tuning parameter. Minimum coefficient norm is parameterization dependent.
* `"ridge"`: regularize a singular system in scaled coordinates, using `control$ridge_strength=1e-8` by default. GMM adds the dimensionless penalty to the column-normalized weighted Gram matrix; H_s adds it after diagonal normalization. Changing a covariate's units therefore does not change fitted effects through an arbitrary absolute ridge value. The applied strength and scales are returned in `fit$regularization`.

An explicit method overrides the conditional default. Regularization is not evidence of identification. When a singular system is regularized, the package returns the selected point solution but withholds ordinary covariance, SEs and Wald intervals; `inference_status` explains this. Separate fixed-operator sensitivity calculations are numerical diagnostics, not sampling covariance for an unidentified estimand. Contradictory calibration targets and infeasible Q support remain errors under every method.

`tilting_diagnostics=FALSE` by default. When enabled, `diagnostics$sufficiency` implements Appendix A.4.1 under unweighted Q: it projects each subgroup indicator and indicator-times-CATE-basis term onto h+. It reports `R_s`, residual norms and unresolved columns. This is a sufficient span check; failure does not disprove A6, and the diagnostic cannot verify population support dominance. Skipped diagnostics return NULL.

## Time and target subgroups

Time modeling is always opt-in. With `time_indexed=TRUE`, use, for example, `~ factor(Time) * (x + z)` for separate horizon coefficients, or `~ factor(Time) + x + z` for shared effect modification. `target_times` selects observed horizons; the package does not extrapolate to unreported times. A formula without Time explicitly pools across horizons.

Repeated horizons within the same trial require a complete joint `Sigma` including cross-time covariance. Marginal arm reports cannot identify that covariance. `variance=FALSE` permits point estimation without it for a one-step common-CATE fit. This requirement also applies to pooled models with repeated reports. Target subgroup functions can define joint subgroups; their effects and all cross-subgroup/cross-time covariances are returned together, with pointwise Wald intervals. With MVLP, every trial must identify all CATE terms, including time terms, from its own reports; each horizon and target subgroup is then pooled separately.

## Migration, data and development

Canonical arguments now use `data_agd`, `formula_cate` and `formula_tilting`; the original camel-case aliases remain accepted. TE/SE-only input is insufficient for the requested variance reconstruction: add arm reports. The old effect-only `V` argument is replaced by joint `Sigma`; use `W` for a custom GMM weight. The old exported implementation helpers (`gmmFun`, `tiltingModel`, `vcovFun`) have been replaced by the documented diagnostics, preparation function, and standard `coef`, `vcov`, `predict` methods.

RDS copies are in `data/` for development and installed under `system.file("extdata", ..., package="CIMAgD")`. Standard `data(dataHF)` and `data(sumstatsPULSE)` also work. Examples and tests depend only on these packaged data; no external raw file or data-creation script is needed.

```r
devtools::test()
# Or build the source package and run R CMD check --no-manual.
```

The [version history](inst/doc/version.md) records implemented changes; the [validation report](inst/doc/validation.md) records the checks against the current version. The package implements the main estimator and Appendices A.1, A.3, A.4.1, A.6.1, A.6.3 and A.7 of the paper. Indirect treatment comparisons (Section 2.4), relative effects (A.2), the misspecification-robust variance (A.6.2) and the conditional arm-mean estimator (A.8) are not implemented; A.8-type working control means are used only for covariance reconstruction.
