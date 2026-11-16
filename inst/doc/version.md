# CIMAgD version history

## 0.1.1 — 2026-11-16

This version follows the November 2026 version of the paper. The statistical scope is the main linear-basis estimator and Appendices A.1, A.3, A.4.1, A.6.1, A.6.3 and A.7. Indirect treatment comparisons (Section 2.4), relative effects (A.2), the misspecification-robust variance (A.6.2) and the conditional arm-mean estimator (A.8) are intentionally not implemented; A.8-type working control means are used only for the A.6.3 covariance reconstruction.

### 1. Check against the updated theory

* Rechecked the normalized tilting representers, linear-basis GMM, moment, tilting and summary Jacobians, Q and trial influence functions, A.6.1 variance, target/Q reuse, subgroup and time-indexed target effects, and the A.6.3 working covariance and variance anchoring against the paper. They agree with it and are numerically unchanged for the common method.
* The paper's new material is the misspecification-robust variance (A.6.2, not implemented), pooling of trial-specific transported effects (A.7, section 2 below), and the common working control mean in A.6.3 (section 3 below).
* Updated appendix references in help pages, comments, method labels and the README to the paper's current numbering: covariance reconstruction is A.6.3 and conditional arm means are A.8.

### 2. Minimum-variance linear pooling (A.7)

* Added `method = c("common", "MVLP")` to `CIMA()`. `"common"` remains the default and is unchanged. The argument is placed after `gmm_method`, so positional calls up to `gmm_method` are unaffected; later arguments are normally named. Both methods use the same input format.
* `"MVLP"` fits a separate CATE for every trial from that trial's moments only, marginalizes each over the same target sample, and pools the trial-specific transported effects with weights `V^{-1} 1 / (1' V^{-1} 1)` and variance `1 / (1' V^{-1} 1)`.
* The trial-specific parameters are stacked with a block-diagonal moment Jacobian and weight. This reproduces the separate trial-specific GMM fits exactly and lets the existing influence-function code return the joint covariance V of the trial-specific effects, including the shared target and Q terms of A.7 (or their combined contribution when Q reuses the target sample). Copula bootstrap replaces the analytic Q component for the trial-specific effects before pooling.
* Every target estimand, i.e., the overall effect and each target subgroup at each horizon, is pooled with weights from its own trial-by-trial covariance block. The joint covariance of the pooled estimates is `L' V L` for the fixed pooling map `L`, because estimating the weights has no first-order effect.
* Identification is checked per trial: each trial's moment Jacobian must have full column rank, and a singular trial stops the fit with a `cima_singularity` error naming that trial unless a singularity method is chosen. `diagnostics$identification` reports the rank and conditioning of every trial.
* Default one-step weights are inverse `n_s * seTE^2`, with the identity fallback decided within each trial. A supplied `W` contributes its trial-diagonal blocks. Two-step GMM inverts each trial's own moment covariance block.
* The pooling weights require V, so MVLP computes the covariance even with `variance = FALSE` and then withholds it. A singular V stops the fit unless a singularity method is chosen; then Moore-Penrose weights are used and Wald inference is withheld.
* Output: `target$estimates` and `target$vcov` hold the pooled effects; `target$trial_estimates` (with pooling weights) and `target$trial_vcov` hold the trial-specific transported effects; `pooling` holds the weights for each estimand; `ATE$by_trial` lists the trial-specific overall effects. `coef()` returns a term-by-trial matrix, `predict()` a column per trial, and `vcov(fit, type = "trial")` the trial-specific covariance. `CATE$theta` and `vcov(fit)` use stacked `trial:term` names.

### 3. Automatic working control mean for covariance reconstruction (A.6.3)

* Each trial first fits its trial-specific working control mean, as before. The fit is flagged when it fails, drops unidentified working terms, has a weighted design with scaled condition number above `control$condition_limit`, needs active risk bounds, or (continuous outcomes) estimates a zero residual variance.
* For a flagged trial only, a common working mean is fitted to the reported control means of all trials, each aggregated under its own tilted law and weighted by inverse squared SEs. It keeps the receiving trial's working basis, its CATE-compatibility bounds and its residual-variance reports. The common fit is used when its issues are no more severe, compared in the order failure, unidentified terms, conditioning, bounds, residual variance; otherwise the trial-specific fit is retained. With a single trial, no common fit is attempted.
* This is automatic rather than a user option and applies to both methods. With MVLP the CATEs remain strictly trial-specific, but the working mean may borrow information from other trials because it only approximates missing covariances. The variance remains a working approximation in either case.
* `variance$working_control` records the model, detected issues and reason for every trial, and `print()` lists trials that use the common mean. In the HF example, DAPA-HF and EMPEROR-Reduced, which report no LVEF subgroups, now use the common working mean; point estimates are unchanged.

### 4. Other changes

* Corrected the README's HF `mu` vectors: the proportions with a history of heart-failure hospitalization are now .261 for DELIVER and .274 for DAPA-HF, as implied by `dataHF` (previously .405 and .474). The packaged data are unchanged. The README no longer states that the HF summaries disagree with the event-analysis denominators beyond rounding.
* `print()` shows the moment-Jacobian rank from the solver even when detailed diagnostics are disabled.
* Added tests for MVLP and for the automatic working control mean; the HF test uses the corrected summaries.

## 0.1.0 — update dated 2026-10-18

This entry records changes from the 2026-10-16 development snapshot of 0.1.0. The package version remains 0.1.0. The statistical scope remains the main linear-basis estimator and Appendices A.1, A.3, A.4.1, A.6.1 and A.6.3 of the paper.

### 1. Optional diagnostics and singular-system handling

* Added `jacobian_diagnostics=TRUE` to `CIMA()`. Setting it to FALSE omits detailed Jacobian reports and conditioning warnings; the solver still performs necessary numerical rank checks.
* Added `singularity_method`, accepting `"none"`, `"minimum-norm"` and `"ridge"`. Its default is `"none"` when Jacobian diagnostics are enabled and `"minimum-norm"` otherwise. An explicit choice always takes precedence. `"none"` stops on singularity.
* Minimum-norm uses Moore-Penrose solutions in the original parameter coordinates. It is the default when diagnostics are disabled because it introduces neither a penalty parameter nor shrinkage in estimable directions. Minimum coefficient norm depends on the chosen parameter units.
* Ridge is applied only when a system is numerically singular. For GMM, write `A = chol(W) J_m_theta`, `D = diag(sqrt(colSums(A^2)))`, and `A_star = A D^{-1}`. The solve uses `(A_star' A_star + lambda I)^{-1} A_star'`, mapped back through `D^{-1}`. Equivalently, the unscaled Gram penalty is `lambda D^2`, so rescaling a nonzero basis column does not change fitted effects. `control$ridge_strength=1e-8` is dimensionless relative to unit column information. Zero columns receive unit numerical scale and coefficient zero.
* For H_s, ridge is applied after diagonal congruence scaling, equivalent to adding `lambda * diag(diag(H_s))` for nonzero columns. Newton stabilization still solves the original calibration equations: it does not silently turn contradictory reported moments into approximate matches. Incompatible redundant targets and infeasible support remain errors.
* Regularization selects a numerical solution; it does not establish the paper's full-rank identification assumptions. If a singular CATE/tilting system or informative singular optimal-weight system is regularized, ordinary covariance, SEs and Wald intervals are withheld. Fixed-operator numerical sensitivity covariance is retained separately, explicitly labeled as non-inferential. Use an identifiable formulation to obtain the paper's standard inference.

### 2. Optional sufficient-span diagnostic

Added `tilting_diagnostics=FALSE`. Enabling it computes the Q-based projection diagnostic in Appendix A.4.1. Skipped diagnostics return NULL and are not printed as a successful check. Residual directions identify what the tilting summaries do not preserve automatically; they do not disprove A6.

### 3. Gaussian copula inputs and marginal families

* Added `copula_summary(n, marginals, correlation=NULL, label=NULL)`, with named marginal lists as a concise input format. Existing marginal tables and the earlier single-source summary format remain supported.
* Added Gamma, Beta, lognormal and exponential margins to normal and Bernoulli margins. Means and SDs can specify normal, Gamma, Beta and lognormal margins. Native distribution parameters are also supported where appropriate. The SD is the distribution's SD, not the SE of its estimated mean.
* Correlations remain latent Gaussian correlations. Marginal parameter validation, probability-integral transformations and numerically stable tail calculations now cover all supported families.
* Plugin uncertainty conditions on the supplied parameters and uses n iid draws. Parametric bootstrap simulates n source observations, refits marginal parameters and latent dependence, then refits downstream estimation under the assumed Gaussian copula model. Neither method requires a sampling covariance matrix for the supplied point estimates.
* Refit rules are explicit: sample mean/SD for normal; sample proportion for Bernoulli; moment matching for Gamma and Beta; arithmetic moment matching for arithmetic lognormal inputs or log-scale moments for meanlog/sdlog inputs; inverse sample mean for exponential rate. Continuous latent correlations use transformed Gaussian scores, continuous/Bernoulli pairs use conditional probit likelihood, and Bernoulli pairs use threshold likelihood. These rules define the propagated uncertainty; unavailable original estimator-specific uncertainty is not inferred.

### 4. Optional two-step GMM

* Added `gmm_method="one-step"` by default and `"two-step"` as an option. The second-stage weight uses the estimated covariance of the complete moment vector, including trial effect/summary cross-covariances and the combined Q influence from Appendix A.6.1.
* Covariance reconstruction is required for two-step weighting even when final variance reporting is disabled. At the second-stage coefficients, the final influence covariance is recomputed using the fixed second-stage weight.
* Singular covariance directions that contain no parameter information are omitted using a range inverse. Informative zero-variance directions are not silently discarded: they require a one-step fit or an explicit ridge choice. The returned GMM record includes both stages, covariance rank, inversion method and any ridge strength.
* With copula bootstrap, the pilot's Q moment covariance comes from refitted copulas, replacing the analytic integration-sample Q term. The second-stage weight is held fixed in subsequent bootstrap coefficient fits. Under valid population GMM moments, estimating this weight has no first-order effect; this avoids nested bootstrap computation. Trial and independent-target variance components remain analytic.

### 5. Computation

Cached common tilting bases and repeated subgroup indicators, reused the CATE basis for exact Q/target reuse, and removed repeated tall-matrix rank decompositions where an existing decomposition suffices. The optional span diagnostic avoids allocating its large residual matrices by default. No convergence tolerance, covariance term or bootstrap refit is removed to obtain a speed improvement.

## 0.1.0 — initial revision dated 2026-10-16

Relative to 0.0.1, this revision introduced normalized subgroup representers, modular input validation, linear-basis GMM, joint Q/trial/target variance propagation, A.6.3 working-control covariance reconstruction, subgroup and optional time-indexed effects, single-source copula propagation, packaged HF/PULSE data, and tests against the paper's formulas. Legacy main-function argument aliases were retained; the effect-only covariance argument and exported implementation helpers were retired. The 2026-10-18 entry above supersedes the initial diagnostic defaults and extends numerical solving and copula support.
