# CIMAgD version history

## 0.1.0 — update dated 2026-10-18

This entry records changes from the 2026-10-16 development snapshot of 0.1.0. The package version remains 0.1.0. The statistical scope remains the main linear-basis estimator and Appendices A.1, A.3, A.4.1, A.6.1 and A.6.2 of the paper.

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

### 6. Documentation

Replaced the planning document with this version history. The validation document maps checks to these changes. Package help, examples, comments and README now refer to “the paper,” retain relevant Appendix references, and describe assumptions for general research use. Added executable examples for the copula specification constructor and the main estimator.

### 7. Packaged data

Verified equivalence of the RDA objects and development/installed RDS copies. Runtime examples and tests depend only on those packaged objects. Removed the redundant raw CSV and its creation script. All 42 HF reports and the PULSE numerical summary values are retained. Dataset descriptions and the PULSE label now use general research terminology.

## 0.1.0 — initial revision dated 2026-10-16

Relative to 0.0.1, this revision introduced normalized subgroup representers, modular input validation, linear-basis GMM, joint Q/trial/target variance propagation, A.6.2 working-control covariance reconstruction, subgroup and optional time-indexed effects, single-source copula propagation, packaged HF/PULSE data, and tests against the paper's formulas. Legacy main-function argument aliases were retained; the effect-only covariance argument and exported implementation helpers were retired. The 2026-10-18 entry above supersedes the initial diagnostic defaults and extends numerical solving and copula support.
