#' Causally interpretable meta-analysis of aggregate trial reports
#'
#' Fits a linear-in-parameters additive CATE with normalized exponentially tilted
#' subgroup moments and averages it over a target population. Implements the
#' main estimator and Appendices A.1, A.3, A.4.1, A.6.1, A.6.3 and A.7 of the
#' paper. `method = "common"` (default) estimates one CATE common to all trials
#' from their stacked moments. `method = "MVLP"` estimates a separate CATE for
#' each trial, marginalizes each over the target population, and combines the
#' trial-specific transported effects by minimum-variance linear pooling (A.7).
#' Relative effects (A.2), the misspecification-robust variance (A.6.2), the
#' conditional arm-mean estimator (A.8) and indirect treatment comparisons are
#' not implemented; A.8-type working control means serve only the A.6.3
#' covariance reconstruction.
#'
#' @param data_agd Contrast reports; see [prepare_agd()].
#' @param mu Named list of trial covariate moment vectors. Named vectors must
#'   match formula_tilting model-matrix columns; omit the fixed intercept or give
#'   it as 1. Unnamed vectors follow model-matrix order.
#' @param Q Empirical base-distribution data frame. Supply either Q or summaryX.
#' @param target Target covariate data frame. NULL reuses Q, with joint inference.
#' @param formula_cate One-sided linear basis formula, e.g. ~ x + z + x:z.
#' @param formula_tilting One-sided formula with intercept, or a trial-named list.
#' @param summaryX Single-source Gaussian copula specification; see [XCopula()].
#' @param q_method Copula uncertainty: plugin (default; fixed parameters and n
#'   iid draws), or bootstrap (refit margins and latent correlations each draw).
#' @param outcome auto, binary, or continuous.
#' @param trial_info Optional data frame with trial, n (total trial size), and
#'   allocation (treatment randomization probability). Otherwise inferred from
#'   overall reports or a disjoint exhaustive partition verified on Q.
#' @param Sigma Optional trial-named list of full joint INFLUENCE covariance
#'   matrices, ordered as that trial's TE rows in input order followed by h+ in
#'   model-matrix order (including the fixed intercept). This is n_s times the
#'   sampling covariance, and the intercept row/column must be zero. Missing
#'   trials use A.6.3 reconstruction. Repeated horizons require joint Sigma.
#' @param mu_se Optional trial-named vectors of covariate-summary SEs, excluding
#'   the intercept, in original tilting-basis order. NA uses working variances.
#' @param W Optional positive-definite GMM weight in report-row order. Default
#'   diagonal entries are inverse n_s*seTE^2; if any SE is zero, the default is
#'   identity (decided within each trial for MVLP). MVLP uses only the
#'   trial-diagonal blocks of a supplied W.
#' @param formula_control Optional working control regression formula or named
#'   list, used only for reconstructed variance (A.6.3). Default is the CATE
#'   basis. Each trial first fits its own working control mean. A trial whose
#'   own fit fails or is unstable (unidentified working terms, ill-conditioning,
#'   active risk bounds, or zero continuous residual variance) automatically
#'   uses a common working mean fitted to every trial's control-arm reports, if
#'   that fit is no less stable; other trials keep their own fits. The choice is
#'   recorded in `variance$working_control`.
#' @param match_subgroup_proportions Add independent subgroup proportions to
#'   trial tilting moments. Default FALSE; conflicting counts/mu cause errors.
#' @param time_indexed Enable Time in the CATE basis; FALSE by default, even if
#'   Time differs. TRUE with a formula omitting Time explicitly pools effects.
#' @param target_times Observed horizons at which to estimate target effects;
#'   defaults to all observed times if time_indexed, otherwise 1.
#' @param target_subgroups Named list of functions of target returning logical
#'   membership, including joint subgroups. Overall is always included. With
#'   MVLP, each subgroup effect is pooled across trials with its own weights.
#' @param q_target_shared NULL detects exact Q/target reuse. TRUE declares exact
#'   row-aligned reuse; FALSE declares independent samples. Partial overlap is
#'   not inferred or handled by the independent-sample variance formula.
#' @param jacobian_diagnostics Report Jacobian rank/conditioning diagnostics;
#'   TRUE by default. Necessary numerical rank checks remain active when FALSE.
#' @param singularity_method Singular-system handling: none, minimum-norm, or
#'   ridge. NULL selects none when jacobian_diagnostics=TRUE, minimum-norm
#'   otherwise. Explicit none always stops on singularity. Minimum-norm uses
#'   Moore-Penrose solutions in the original parameter coordinates. Ridge uses
#'   dimensionless penalties after scaling, only for singular systems.
#' @param tilting_diagnostics Compute the A.4.1 sufficient-span diagnostic;
#'   FALSE by default. It cannot establish that A6 fails.
#' @param gmm_method One-step (default) or two-step. Two-step uses the estimated
#'   full moment sampling covariance to construct its second-stage weight; with
#'   MVLP, each trial inverts its own moment covariance block. Redundant
#'   covariance directions are removed; informative zero-variance directions
#'   require ridge or a one-step analysis. Covariance estimation is needed for
#'   two-step even if variance=FALSE.
#' @param method `"common"` (default): one CATE shared by all trials (A1-A6).
#'   `"MVLP"`: trial-specific CATEs, each identified from its own trial's
#'   reports, transported to the target and combined by minimum-variance linear
#'   pooling (A.7: A1p, A2p, and A4-A6 within each trial). Every trial must
#'   identify the CATE basis; the moment-Jacobian rank is checked per trial.
#'   Pooling weights need the joint covariance of the trial-specific estimates,
#'   which is therefore computed even when variance=FALSE.
#' @param variance Calculate variance (default TRUE). FALSE permits point-only
#'   one-step common-CATE fits when cross-time covariance or feasible working
#'   risk models are absent. Two-step and MVLP fits still compute the covariance
#'   for their weights but do not report it.
#' @param bootstrap_reps Number of copula refit replicates, default 200.
#' @param bootstrap_mc_size Numerical integration size for bootstrap, default
#'   max(10000, 2*summaryX$n). Common Gaussian draws are reused for integration
#'   across refits to reduce numerical noise. Only bootstrap uses this setting.
#' @param seed Optional local RNG seed.
#' @param control Named numerical controls: rank_tol=1e-9, condition_limit=1e8,
#'   tilt_tol=1e-10, span_tol=1e-8, max_iter=200, ridge_strength=1e-8. The latter
#'   is dimensionless; GMM adds it to the column-normalized Gram matrix and H
#'   adds it to its diagonal-normalized Gram matrix. condition_limit also flags
#'   ill-conditioned working control regressions.
#' @param dataAgD,formulaCATE,formulaTilting Legacy argument aliases.
#' @param V Retired treatment-effect-only covariance argument. Use joint Sigma.
#' @return A CIMA object with CATE, ATE, joint target estimates, normalized
#'   reports, tilting fits, Jacobians, influence/variance components, diagnostics
#'   and copula metadata. `diagnostics$identification` gives the moment-Jacobian
#'   rank (one row per trial for MVLP) and `variance$working_control` the working
#'   control model used for each trial. For MVLP, `CATE$theta` stacks the
#'   trial-specific coefficients (named `trial:term`) and [coef()] returns them
#'   as a term-by-trial matrix; `target$estimates` and `target$vcov` hold the
#'   pooled target effects; `target$trial_estimates` and `target$trial_vcov` the
#'   trial-specific transported effects with their joint covariance; `pooling`
#'   the minimum-variance weights for each estimand; and `ATE$by_trial` the
#'   trial-specific overall effects. Reconstruction is a working approximation
#'   under raw complete-outcome randomized-trial influence functions, not exact
#'   inference for arbitrary adjusted or censored estimates. Continuous
#'   reconstruction assumes arm-specific constant conditional residual variances.
#'   If a singular system is regularized, ordinary identification-based
#'   covariance and Wald intervals are withheld; numerical sensitivity covariance
#'   is stored separately.
#' @examples
#' reports <- data.frame(
#'   trial = "s1", subgroup = c("overall", "x", "x"),
#'   stratum = c("TRUE", "x == 0", "x == 1"),
#'   event1 = c(30, 10, 20), n1 = c(100, 50, 50),
#'   event2 = c(20, 5, 15), n2 = c(100, 50, 50)
#' )
#' fit <- CIMA(reports, mu = list(s1 = c(x = .5)),
#'   Q = data.frame(x = rep(0:1, each = 100)),
#'   formula_cate = ~ x, formula_tilting = ~ x)
#' coef(fit)
#' fit$target$estimates
#'
#' # Trial-specific CATEs pooled by minimum-variance linear pooling
#' two <- rbind(reports, data.frame(
#'   trial = "s2", subgroup = c("overall", "x", "x"),
#'   stratum = c("TRUE", "x == 0", "x == 1"),
#'   event1 = c(45, 20, 25), n1 = c(150, 90, 60),
#'   event2 = c(30, 12, 18), n2 = c(150, 90, 60)))
#' pooled <- CIMA(two, mu = list(s1 = c(x = .5), s2 = c(x = .4)),
#'   Q = data.frame(x = rep(0:1, each = 100)),
#'   formula_cate = ~ x, formula_tilting = ~ x, method = "MVLP")
#' coef(pooled)
#' pooled$target$trial_estimates
#' pooled$target$estimates
#' @export
CIMA <- function(data_agd = NULL, mu, Q = NULL, target = NULL,
                 formula_cate = NULL, formula_tilting = NULL,
                 summaryX = NULL, q_method = c("plugin", "bootstrap"),
                 outcome = c("auto", "binary", "continuous"), trial_info = NULL,
                 Sigma = NULL, mu_se = NULL, W = NULL, formula_control = NULL,
                 match_subgroup_proportions = FALSE, time_indexed = FALSE,
                 target_times = NULL, target_subgroups = list(), q_target_shared = NULL,
                 jacobian_diagnostics = TRUE, singularity_method = NULL,
                 tilting_diagnostics = FALSE, gmm_method = c("one-step", "two-step"),
                 method = c("common", "MVLP"),
                 variance = TRUE, bootstrap_reps = 200L, bootstrap_mc_size = NULL,
                 seed = NULL, control = list(), dataAgD = NULL,
                 formulaCATE = NULL, formulaTilting = NULL, V = NULL) {
  cl <- match.call()
  if (!is.null(V)) stop("V contains only effect covariance. Supply the complete joint influence covariance via Sigma, or let CIMA reconstruct it from arm reports.", call. = FALSE)
  alias <- function(primary, old, name) {
    if (!is.null(primary) && !is.null(old)) stop("Supply only one spelling of ", name, ".", call. = FALSE)
    primary %||% old
  }
  data_agd <- alias(data_agd, dataAgD, "data_agd")
  formula_cate <- alias(formula_cate, formulaCATE, "formula_cate")
  formula_tilting <- alias(formula_tilting, formulaTilting, "formula_tilting")
  if (is.null(data_agd) || is.null(formula_cate) || is.null(formula_tilting))
    stop("data_agd, formula_cate and formula_tilting are required.", call. = FALSE)
  for (nm in c("variance", "time_indexed", "match_subgroup_proportions", "jacobian_diagnostics", "tilting_diagnostics")) {
    value <- get(nm)
    if (!is.logical(value) || length(value) != 1L || is.na(value)) stop(nm, " must be TRUE or FALSE.", call. = FALSE)
  }
  if (!is.null(q_target_shared) && (!is.logical(q_target_shared) || length(q_target_shared) != 1L || is.na(q_target_shared))) stop("q_target_shared must be NULL, TRUE or FALSE.", call. = FALSE)
  defaults <- list(rank_tol = 1e-9, condition_limit = 1e8, tilt_tol = 1e-10, span_tol = 1e-8, max_iter = 200L, ridge_strength = 1e-8)
  if (length(setdiff(names(control), names(defaults)))) stop("Unknown numerical control.", call. = FALSE)
  ctl <- utils::modifyList(defaults, control)
  for (nm in setdiff(names(ctl), "max_iter")) check_scalar(ctl[[nm]], nm, .Machine$double.eps)
  check_scalar(ctl$max_iter, "max_iter", 1, TRUE)
  if (ctl$rank_tol >= 1) stop("rank_tol must be less than one.", call. = FALSE)
  singularity_method <- singularity_method %||% if (jacobian_diagnostics) "none" else "minimum-norm"
  singularity_method <- match.arg(singularity_method, c("none", "minimum-norm", "ridge"))
  ctl$singularity_method <- singularity_method
  ctl$jacobian_diagnostics <- jacobian_diagnostics
  gmm_method <- match.arg(gmm_method)
  method <- match.arg(method)
  trial_specific <- method == "MVLP"
  q_method <- match.arg(q_method); outcome <- match.arg(outcome)
  d <- prepare_agd(data_agd, outcome)
  outcome <- attr(d, "outcome")
  if (is.null(summaryX) == is.null(Q)) stop("Supply exactly one of empirical Q or copula summaryX.", call. = FALSE)
  if (is.null(summaryX) && q_method == "bootstrap") stop("q_method='bootstrap' requires summaryX.", call. = FALSE)
  for (nm in c("Sigma", "mu_se")) {
    value <- get(nm)
    if (!is.null(value) && (!is.list(value) || is.null(names(value)) || anyDuplicated(names(value)) || any(!names(value) %in% d$trial)))
      stop(nm, " must be a uniquely trial-named list.", call. = FALSE)
  }
  if (!is.list(target_subgroups)) stop("target_subgroups must be a named list of functions.", call. = FALSE)
  with_seed(seed, {
    spec_q <- integration_noise <- NULL
    if (!is.null(summaryX)) {
      spec_q <- validate_summary(summaryX)
      size <- spec_q$n
      if (q_method == "bootstrap") {
        check_scalar(bootstrap_reps, "bootstrap_reps", 2, TRUE)
        size <- bootstrap_mc_size %||% max(10000L, 2L * spec_q$n)
        check_scalar(size, "bootstrap_mc_size", 2, TRUE)
      }
      integration_noise <- matrix(stats::rnorm(size * nrow(spec_q$marginals)), size)
      Q <- copula_from_noise(spec_q, integration_noise)
    }
    Q <- as.data.frame(Q)
    if (nrow(Q) < 2L || anyNA(Q)) stop("Q needs at least two complete observations.", call. = FALSE)
    target_default <- is.null(target)
    if (target_default) target <- Q
    target <- as.data.frame(target)
    if (nrow(target) < 2L || anyNA(target)) stop("target needs at least two complete observations.", call. = FALSE)
    identical_samples <- identical(Q, target)
    shared <- q_target_shared %||% identical_samples
    if (target_default && !shared) stop("target=NULL reuses Q; supply a separate target to declare independence.", call. = FALSE)
    if (shared && !identical_samples) stop("Shared target/Q requires identical, row-aligned data; partial overlap is not supported.", call. = FALSE)
    times <- if (time_indexed) sort(unique(d$Time)) else 1
    target_times <- target_times %||% times
    if (!is.numeric(target_times) || any(!is.finite(target_times)) || anyDuplicated(target_times) || !length(target_times) || any(!target_times %in% times))
      stop("target_times must be distinct observed horizons (1 when time_indexed=FALSE).", call. = FALSE)
    core <- fit_cima_core(d, Q, mu, formula_cate, formula_tilting, trial_info,
      match_subgroup_proportions, time_indexed, times, W, ctl, diagnostics = tilting_diagnostics,
      method = method)
    fit <- core$fit; context <- core$context
    two_step <- gmm_method == "two-step"
    # MVLP needs the covariance of the trial-specific estimates for its weights.
    keep_covariance <- variance || trial_specific
    vc_pilot <- if (keep_covariance || two_step) variance_components(context, fit, Sigma, mu_se, formula_control, ctl) else NULL
    gmm <- list(method = gmm_method, steps = if (two_step) 2L else 1L,
                initial_theta = fit$theta, initial_W = fit$W, optimal = NULL)
    bootstrap_q <- !is.null(spec_q) && q_method == "bootstrap"
    prepared_bootstrap <- NULL
    if (two_step) {
      covariance <- vc_pilot$moment_covariance
      if (bootstrap_q) {
        prepared_bootstrap <- prepare_copula_bootstrap(spec_q, integration_noise, bootstrap_reps,
          d, mu, formula_cate, formula_tilting, trial_info, match_subgroup_proportions,
          time_indexed, times, W, ctl, fit$theta, method)
        covariance <- prepared_bootstrap$moment_covariance + vc_pilot$trial_moment_covariance
      }
      optimal <- optimal_weight_blocks(covariance, context$moments$J_m_theta, fit$blocks, ctl)
      optimal$source <- if (is.null(prepared_bootstrap)) "analytic Q and trial covariance" else "refitted copula bootstrap and trial covariance"
      fit <- fit_gmm(d, context$moments, context$trial_info, optimal$W, ctl, weight_psd = TRUE, blocks = fit$blocks)
      fit$weight_note <- paste("Two-step GMM:", optimal$method)
      gmm$optimal <- optimal
    }
    vc <- if (two_step && keep_covariance) variance_components(context, fit, Sigma, mu_se, formula_control, ctl) else if (keep_covariance) vc_pilot else NULL
    target_fit <- fit_target(target, context$basis, fit, target_times, target_subgroups, vc, shared)
    bootstrap <- NULL
    if (bootstrap_q && keep_covariance) {
      bootstrap <- bootstrap_copula(spec_q, integration_noise, bootstrap_reps,
        d, mu, formula_cate, formula_tilting, trial_info, match_subgroup_proportions,
        time_indexed, times, if (two_step) fit$W else W, ctl, target, shared, target_times, target_subgroups,
        prepared = prepared_bootstrap, weight_psd = two_step, theta_reference = fit$theta, method = method)
      vc$q_covariance_plugin <- vc$q_covariance
      vc$q_covariance <- bootstrap$theta_covariance
      vc$vcov_theta <- symmetrize(vc$trial_covariance + vc$q_covariance)
      theta_target_trial <- target_fit$J_psi_theta %*% vc$trial_covariance %*% t(target_fit$J_psi_theta)
      independent_target <- if (shared) matrix(0, nrow(target_fit$vcov), ncol(target_fit$vcov)) else
        crossprod(center_columns(target_fit$phi_target)) / nrow(target)^2
      target_fit$vcov <- symmetrize(bootstrap$target_covariance + theta_target_trial + independent_target)
      dimnames(target_fit$vcov) <- dimnames(bootstrap$target_covariance)
      target_fit$components <- list(copula_bootstrap = bootstrap$target_covariance,
        trials = theta_target_trial, independent_target = independent_target)
      target_fit$estimates <- add_wald(target_fit$estimates, target_fit$vcov)
      vc$q_method <- "bootstrap replaces analytic Q covariance"
      vc$moment_covariance_plugin <- vc$moment_covariance
      vc$moment_covariance <- bootstrap$moment_covariance + vc$trial_moment_covariance
    }
    regularized <- fit$regularization$used || any(vapply(context$tilting,
      function(x) x$regularization$used, logical(1))) || (two_step && gmm$optimal$ridge > 0)
    pooling <- NULL
    if (trial_specific) {
      # Weights use the full trial-specific covariance, before any inference is withheld.
      pooling <- pool_target(target_fit, names(fit$blocks), ctl)
      regularized <- regularized || pooling$regularized
    }
    inference <- if (regularized) "Singular system: selected numerical solution; standard Wald inference unavailable" else
      if (variance) paste0("Paper's full-rank influence-function inference",
        if (trial_specific) " with A.7 minimum-variance linear pooling" else "") else "Variance not requested"
    if (regularized && variance) {
      vc$sensitivity_vcov_theta <- vc$vcov_theta
      vc$sensitivity_note <- "Fixed-operator numerical sensitivity only; not a sampling covariance for an unidentified or regularized estimand."
      vc$vcov_theta[,] <- NA_real_
      target_fit <- withhold_inference(target_fit)
      if (trial_specific) pooling <- withhold_inference(pooling)
      warning("A singular system was regularized. Standard covariance and Wald intervals are unavailable; see numerical sensitivity covariance and inference_status.", call. = FALSE)
    }
    if (!variance && trial_specific) {
      # The covariance served only the pooling weights.
      target_fit <- withhold_inference(target_fit, keep = FALSE)
      pooling <- withhold_inference(pooling, keep = FALSE)
      vc <- NULL
    }
    P <- length(fit$theta)
    theta_cov <- if (variance) vc$vcov_theta else matrix(NA_real_, P, P)
    dimnames(theta_cov) <- list(names(fit$theta), names(fit$theta))
    theta_se <- sqrt(pmax(diag(theta_cov), 0))
    CATE <- list(theta = fit$theta, vcov = theta_cov, se = theta_se,
      ci95 = cbind(lower = fit$theta - stats::qnorm(.975) * theta_se,
                   upper = fit$theta + stats::qnorm(.975) * theta_se),
      J_m_theta = context$moments$J_m_theta, Dtheta = context$moments$J_m_theta)
    target_out <- target_fit
    if (trial_specific) {
      CATE$coefficients <- matrix(fit$theta, ncol = length(fit$blocks),
        dimnames = list(context$basis$columns, names(fit$blocks)))
      CATE$design <- fit$design
      trial_table <- target_fit$estimates
      key <- paste(trial_table$subgroup, trial_table$Time, sep = "@")
      trial_table$weight <- pooling$map[cbind(seq_len(nrow(trial_table)), match(key, colnames(pooling$map)))]
      target_out <- list(estimates = pooling$estimates, vcov = pooling$vcov,
        J_psi_theta = pooling$J_psi_theta, phi_target = pooling$phi_target,
        components = pooling$components, shared_q = shared, sensitivity_vcov = pooling$sensitivity_vcov,
        trial_estimates = trial_table, trial_vcov = target_fit$vcov,
        trial_components = target_fit$components, trial_J_psi_theta = target_fit$J_psi_theta,
        trial_phi_target = target_fit$phi_target, trial_sensitivity_vcov = target_fit$sensitivity_vcov)
    }
    overall <- which(target_out$estimates$subgroup == "overall")
    ate <- target_out$estimates[overall, , drop = FALSE]
    ATE <- list(ATE = ate$estimate, var = ate$se^2, se = ate$se,
      ci95 = if (nrow(ate) == 1L) c(ate$lower, ate$upper) else cbind(ate$lower, ate$upper), Time = ate$Time)
    if (trial_specific) {
      ATE$by_trial <- trial_table[trial_table$subgroup == "overall",
        c("trial", "Time", "estimate", "se", "lower", "upper", "weight")]
      rownames(ATE$by_trial) <- NULL
    }
    result <- list(method = method, CATE = CATE, ATE = ATE,
      target = target_out, pooling = if (trial_specific) list(
        method = "Minimum-variance linear pooling of trial-specific transported effects (A.7)",
        weights = pooling$weights, status = pooling$status, map = pooling$map,
        covariance_blocks = if (variance) pooling$covariance_blocks else NULL,
        note = paste("Weights are computed separately for every target estimand from its trial-by-trial",
          "covariance; estimating them has no first-order effect on the pooled estimator.")) else NULL,
      data = d, trial_info = context$trial_info,
      moments = list(tau_hat = d$TE, m_hat = fit$m_hat, residual = fit$residual,
        foc = fit$foc, objective = fit$objective, alpha = context$moments$alpha,
        c_sj = context$moments$c_sj, W = fit$W, weight_note = fit$weight_note),
      tilting = context$tilting, variance = vc, gmm = gmm, inference_status = inference,
      regularization = list(method = singularity_method, used = regularized,
        CATE = fit$regularization, tilting = lapply(context$tilting, `[[`, "regularization"),
        pooling = if (trial_specific) pooling$status else NULL),
      diagnostics = list(J_m_theta = fit$diagnostics, identification = fit$identification,
        H_s = if (jacobian_diagnostics) lapply(context$tilting, `[[`, "diagnostics") else NULL,
        sufficiency = core$sufficiency),
      copula = if (is.null(spec_q)) NULL else list(summaryX = spec_q, method = q_method,
        original_n = spec_q$n, integration_n = nrow(Q), bootstrap = bootstrap),
      basis = context$basis, time_indexed = time_indexed, call = cl)
    class(result) <- "CIMA"
    result
  })
}

# Withhold ordinary covariance and Wald intervals, keeping numerical sensitivity
# covariance when requested.
withhold_inference <- function(x, keep = TRUE) {
  if (keep) x$sensitivity_vcov <- x$vcov
  x$vcov[,] <- NA_real_
  x$estimates[c("se", "lower", "upper")] <- NA_real_
  if (!keep) x$components <- NULL
  x
}

fit_cima_core <- function(d, Q, mu, formula_cate, formula_tilting, trial_info,
                         match_proportions, time_indexed, times, W, control, diagnostics,
                         weight_psd = FALSE, method = "common") {
  indicators <- indicator_matrix(Q, d)
  if (any(colSums(indicators) == 0)) stop("No Q support for subgroup: ", paste(d$row_id[colSums(indicators) == 0], collapse = ", "), call. = FALSE)
  info <- resolve_trial_info(d, indicators, trial_info)
  tilt <- fit_tilting(Q, d, indicators, mu, formula_tilting, info, match_proportions, control)
  basis <- make_basis_spec(formula_cate, Q, times, time_indexed)
  moments <- build_moments(d, indicators, tilt, basis)
  fit <- fit_gmm(d, moments, info, W, control, weight_psd, cate_blocks(d, basis$columns, method))
  sufficiency <- NULL
  if (diagnostics) sufficiency <- list()
  if (diagnostics) for (s in names(tilt)) {
    idx <- which(d$trial == s)
    # Separate horizon checks retain labels and avoid irrelevant time-basis rows.
    sufficiency[[s]] <- lapply(unique(moments$row_times[idx]), function(tt) {
      jj <- idx[moments$row_times[idx] == tt]
      diagnose_tilting(tilt[[s]]$h_plus, basis$basis_q[[tt]], indicators[, jj, drop = FALSE], control$span_tol)
    })
    names(sufficiency[[s]]) <- unique(moments$row_times[idx])
  }
  list(fit = fit, sufficiency = sufficiency, context = list(data = d, Q = Q,
    indicators = indicators, trial_info = info, tilting = tilt, basis = basis,
    moments = moments, outcome = attr(d, "outcome")))
}

#' Methods for a CIMA fit
#' @param x,object A fitted CIMA object.
#' @param digits Number of printed decimal places.
#' @param type Covariance of theta, of the (pooled) target estimates, or, for
#'   MVLP, of the trial-specific transported estimates (`"trial"`).
#' @param newdata Covariate data frame for CATE prediction.
#' @param time One observed horizon for a time-indexed fit.
#' @param ... Unused.
#' @return print invisibly returns x; coef returns coefficients (a term-by-trial
#'   matrix for MVLP); vcov returns a covariance matrix; predict returns fitted
#'   additive CATE values (one column per trial for MVLP).
#' @export
print.CIMA <- function(x, digits = 4, ...) {
  trial_specific <- identical(x$method, "MVLP")
  cat("CIMAgD fit:", if (trial_specific) "trial-specific CATEs with minimum-variance linear pooling (MVLP)" else "common CATE", "\n")
  if (trial_specific) cat("\nPooled target estimates\n")
  print(x$target$estimates, digits = digits, row.names = FALSE)
  if (trial_specific) {
    cat("\nTrial-specific transported estimates\n")
    print(x$target$trial_estimates, digits = digits, row.names = FALSE)
    cat("\nCATE coefficients by trial\n")
    print(round(x$CATE$coefficients, digits))
  } else {
    cat("\nCATE coefficients\n")
    print(round(cbind(estimate = x$CATE$theta, se = x$CATE$se), digits))
  }
  id <- x$diagnostics$identification
  if (!is.null(id)) {
    if (trial_specific) cat("\nMoment Jacobian rank by trial:",
      paste0(id$trial, " ", id$rank, "/", id$parameters, collapse = ", "), "\n")
    else cat("\nMoment Jacobian rank:", id$rank, "/", id$parameters, "\n")
  }
  if (!is.null(x$diagnostics$sufficiency)) {
    n_unresolved <- sum(vapply(x$diagnostics$sufficiency, function(s)
      sum(vapply(s, function(t) length(t$unresolved), integer(1))), integer(1)))
    cat("Required aggregate basis columns outside the tilting span:", n_unresolved, "\n")
  }
  wc <- x$variance$working_control
  if (!is.null(wc) && any(wc$model == "common"))
    cat("Common working control mean (A.6.3) used for:", paste(wc$trial[wc$model == "common"], collapse = ", "), "\n")
  cat("GMM:", x$gmm$method, "\n")
  if (x$regularization$used) cat(x$inference_status, "\n")
  invisible(x)
}

#' @rdname print.CIMA
#' @export
coef.CIMA <- function(object, ...) {
  if (identical(object$method, "MVLP")) object$CATE$coefficients else object$CATE$theta
}

#' @rdname print.CIMA
#' @export
vcov.CIMA <- function(object, type = c("theta", "target", "trial"), ...) {
  type <- match.arg(type)
  if (type == "theta") return(object$CATE$vcov)
  if (type == "target") return(object$target$vcov)
  if (!identical(object$method, "MVLP")) stop("type='trial' is available only for method='MVLP'.", call. = FALSE)
  object$target$trial_vcov
}

#' @rdname print.CIMA
#' @export
predict.CIMA <- function(object, newdata, time = NULL, ...) {
  if (object$time_indexed) {
    if (is.null(time) || length(time) != 1L || !time %in% object$basis$times)
      stop("Supply one observed time for prediction.", call. = FALSE)
  } else time <- 1
  b <- new_basis(object$basis, as.data.frame(newdata), time)
  if (!identical(object$method, "MVLP")) return(as.vector(b %*% object$CATE$theta))
  out <- b %*% object$CATE$coefficients
  dimnames(out) <- list(NULL, colnames(object$CATE$coefficients))
  out
}
