# A.6.3 working nuisance regressions only; this is not the A.8 arm-mean estimator.
# Rows of A are tilted subgroup means of the working basis. By default they are the
# receiving trial's own reports. A pooled (common) fit uses every trial's reported
# control means, each aggregated under its own tilted law, but keeps the receiving
# trial's basis, CATE-compatibility bounds and residual-variance reports.
fit_working_control <- function(b, alpha, w, g, reports, outcome, control, pooled = NULL) {
  A_own <- crossprod(alpha, b) / nrow(b)
  if (is.null(pooled)) {
    A <- A_own; nu0 <- reports$nu0; se <- reports$se_nu0
  } else {
    A <- as.matrix(pooled$A); nu0 <- pooled$nu0; se <- pooled$se0
    if (ncol(A) != ncol(b) || nrow(A) != length(nu0) || length(se) != length(nu0))
      stop("Pooled working control rows do not match the working basis.", call. = FALSE)
  }
  reference <- if (any(se > 0)) min(se[se > 0]) else 1
  precision <- 1 / ifelse(se > 0, se^2, reference^2)
  a <- A * sqrt(precision)
  y <- nu0 * sqrt(precision)
  # Choose a smaller identifiable working basis in formula order, retaining an
  # intercept when present. Projecting onto the aggregate row space can exclude
  # even a feasible constant risk; column selection avoids that restriction.
  selected <- integer()
  for (k in seq_len(ncol(a))) {
    candidate <- a[, c(selected, k), drop = FALSE]
    if (diagnose_jacobian(candidate, control$rank_tol)$full_rank)
      selected <- c(selected, k)
  }
  if (!length(selected)) stop("Working control regression has no estimable basis.", call. = FALSE)
  column_scale <- sqrt(colSums(a[, selected, drop = FALSE]^2))
  z <- svd(sweep(a[, selected, drop = FALSE], 2, column_scale, "/"))
  keep <- z$d > max(z$d, 0) * control$rank_tol
  if (!any(keep)) stop("Working control regression has no estimable basis.", call. = FALSE)
  condition <- max(z$d[keep]) / min(z$d[keep])
  transform <- matrix(0, ncol(a), sum(keep))
  transform[selected, ] <- sweep(sweep(z$v[, keep, drop = FALSE], 2, z$d[keep], "/"),
                                 1, column_scale, "/")
  u <- z$u[, keep, drop = FALSE]
  unconstrained <- as.vector(crossprod(u, y))
  b_reduced <- b %*% transform
  constrained <- FALSE
  if (outcome == "binary") {
    if (any(abs(g) > 1 + 1e-8)) stop(binary_cate_message(), call. = FALSE)
    lower <- pmax(0, -g); upper <- pmin(1, 1 - g)
    predicted <- as.vector(b_reduced %*% unconstrained)
    if (any(predicted < lower - 1e-10 | predicted > upper + 1e-10)) {
      constrained <- TRUE
      fit <- tryCatch(quadprog::solve.QP(Dmat = diag(2, sum(keep)),
        dvec = 2 * unconstrained, Amat = t(rbind(b_reduced, -b_reduced)),
        bvec = c(lower, -upper)), error = function(e) e)
      if (inherits(fit, "error")) stop("No feasible working control-risk regression in the chosen basis: ",
        conditionMessage(fit), ". Change formula_control or provide Sigma.", call. = FALSE)
      coefficient <- fit$solution
    } else coefficient <- unconstrained
  } else coefficient <- unconstrained
  beta <- as.vector(transform %*% coefficient); names(beta) <- colnames(b)
  gamma0 <- as.vector(b %*% beta); gamma1 <- gamma0 + g
  if (outcome == "binary") {
    if (min(c(gamma0, gamma1)) < -1e-7 || max(c(gamma0, gamma1)) > 1 + 1e-7)
      stop("Working risk constraints failed numerical verification.", call. = FALSE)
    gamma0 <- pmin(1, pmax(0, gamma0)); gamma1 <- pmin(1, pmax(0, gamma1))
    v0 <- gamma0 * (1 - gamma0); v1 <- gamma1 * (1 - gamma1)
    residual_model <- "Bernoulli conditional variance"
    residual_diagnostics <- NULL
  } else {
    estimate_residual <- function(gamma, se_arm, n_arm) {
      means <- as.vector(crossprod(alpha, gamma) / nrow(alpha))
      centered <- sweep(matrix(gamma, nrow(alpha), ncol(alpha)), 2, means, "-")
      mean_variance <- colMeans(alpha * centered^2)
      residual <- n_arm * se_arm^2 - mean_variance
      variance <- max(0, stats::weighted.mean(residual, n_arm))
      list(variance = variance, subgroup_residual = residual,
           at_boundary = variance == 0, subgroup_mean_variance = mean_variance)
    }
    r0 <- estimate_residual(gamma0, reports$se_nu0, reports$n2)
    r1 <- estimate_residual(gamma1, reports$se_nu1, reports$n1)
    v0 <- rep(r0$variance, length(g)); v1 <- rep(r1$variance, length(g))
    residual_model <- "Arm-specific constant nonnegative conditional residual variance"
    residual_diagnostics <- list(control = r0, treatment = r1)
  }
  aggregate_means <- as.vector(A_own %*% beta)
  list(beta = beta, gamma0 = gamma0, gamma1 = gamma1, v0 = v0, v1 = v1,
       model = if (is.null(pooled)) "trial-specific" else "common",
       aggregate_means = aggregate_means, residual = aggregate_means - reports$nu0,
       fitted_rows = if (is.null(pooled)) reports$row_id else pooled$row_id,
       fitted_residual = if (is.null(pooled)) NULL else as.vector(A %*% beta) - nu0,
       rank = sum(keep), basis_dimension = ncol(b), condition = condition,
       retained_columns = colnames(b)[selected], dropped_columns = colnames(b)[-selected],
       identification = "Independent columns in formula order (intercept retained)",
       constrained = constrained, zero_se_weight_fallback = any(se == 0),
       residual_model = residual_model, residual_diagnostics = residual_diagnostics)
}

binary_cate_message <- function() paste("Fitted binary CATE leaves [-1, 1] on Q; no compatible working arm risks.",
  "Simplify the CATE, supply Sigma, or use variance=FALSE.")

# Stability flags of a working control fit, ordered by severity for comparison.
working_control_issues <- function(x, outcome, control) {
  if (inherits(x, "error"))
    return(list(flags = "estimation failure", severity = c(1, 0, 0, 0, 0), message = conditionMessage(x)))
  dropped <- length(x$dropped_columns)
  ill <- !is.finite(x$condition) || x$condition > (control$condition_limit %||% 1e8)
  boundary <- isTRUE(x$constrained)
  zero_variance <- outcome == "continuous" && isTRUE(x$residual_diagnostics$control$at_boundary ||
                                                     x$residual_diagnostics$treatment$at_boundary)
  flags <- c(if (dropped) "unidentified working terms", if (ill) "ill-conditioned",
             if (boundary) "active risk bounds", if (zero_variance) "zero residual variance")
  list(flags = flags, severity = c(0, dropped, ill, boundary, zero_variance), message = NULL)
}

lexicographic_le <- function(a, b) {
  differ <- which(a != b)
  !length(differ) || a[differ[1L]] < b[differ[1L]]
}

# Automatic A.6.3 working control model. The trial-specific fit is kept whenever
# it is stable. A trial whose own fit fails, drops unidentified terms, is
# ill-conditioned, needs active risk bounds, or (continuous) has zero residual
# variance receives the common working mean instead, provided that fit is no less
# stable. Other trials are unaffected. The working mean only approximates missing
# covariances; it does not enter the CATE estimator.
select_working_control <- function(s, b, alpha, w, g, reports, outcome, control, pooled_rows) {
  if (outcome == "binary" && any(abs(g) > 1 + 1e-8)) stop(binary_cate_message(), call. = FALSE)
  specific <- tryCatch(fit_working_control(b, alpha, w, g, reports, outcome, control), error = identity)
  specific_issues <- working_control_issues(specific, outcome, control)
  chosen <- specific; common_issues <- NULL
  reason <- "trial-specific fit is stable"
  rows <- if (length(specific_issues$flags)) pooled_rows() else NULL
  if (!is.null(rows) && !is.null(rows$trial) && all(rows$trial == s)) {
    rows <- NULL
    reason <- "no other trial reports control arms; trial-specific fit retained"
  }
  if (!is.null(rows)) {
    common <- tryCatch(fit_working_control(b, alpha, w, g, reports, outcome, control, pooled = rows),
                       error = identity)
    common_issues <- working_control_issues(common, outcome, control)
    if (!inherits(common, "error") && lexicographic_le(common_issues$severity, specific_issues$severity)) {
      chosen <- common
      reason <- if (length(common_issues$flags)) "common fit is no less stable than the trial-specific fit" else
        "common fit resolves the trial-specific issues"
    } else reason <- "common fit is not more stable; trial-specific fit retained"
  }
  if (inherits(chosen, "error"))
    stop("No working control model for trial ", s, ". Trial-specific fit: ", specific_issues$message,
         if (is.null(common_issues)) " No common fit is available from other trials." else
           paste(" Common fit:", common_issues$message), call. = FALSE)
  describe <- function(x) if (is.null(x)) NA_character_ else if (length(x$flags)) paste(x$flags, collapse = "; ") else "none"
  chosen$selection <- list(model = chosen$model, reason = reason,
    trial_specific_issues = specific_issues$flags, trial_specific_error = specific_issues$message,
    common_issues = common_issues$flags, common_error = common_issues$message,
    summary = data.frame(trial = s, model = chosen$model, trial_specific_issues = describe(specific_issues),
                         common_issues = describe(common_issues), reason = reason))
  chosen
}

working_joint_covariance <- function(indicators, w, h, fit, allocation) {
  p <- colMeans(indicators * w)
  k <- sweep(indicators, 2, p, "/")
  n_q <- nrow(k)
  hc <- sweep(h, 2, colMeans(h * w), "-")
  arm <- function(gamma, variance, e) {
    nu <- as.vector(crossprod(k, gamma * w) / n_q)
    r <- k * sweep(matrix(gamma, n_q, ncol(k)), 2, nu, "-")
    covariance <- (crossprod(k, k * (w * variance)) + crossprod(r, r * w)) / (n_q * e)
    cross_mu <- crossprod(r, hc * w) / n_q
    list(covariance = covariance, cross_mu = cross_mu, nu = nu)
  }
  a1 <- arm(fit$gamma1, fit$v1, allocation)
  a0 <- arm(fit$gamma0, fit$v0, 1 - allocation)
  tau_cov <- a1$covariance + a0$covariance
  cross_mu <- a1$cross_mu - a0$cross_mu
  mu_cov <- crossprod(hc, hc * w) / n_q
  sigma <- rbind(cbind(tau_cov, cross_mu), cbind(t(cross_mu), mu_cov))
  list(Sigma = symmetrize(sigma), arm1 = a1, arm0 = a0)
}

anchor_covariance <- function(working, variances) {
  missing <- is.na(variances)
  variances[missing] <- pmax(diag(working)[missing], 0)
  original <- pmax(diag(working), 0)
  scale <- ifelse(original > 0, sqrt(variances / pmax(original, .Machine$double.xmin)), 0)
  sigma <- working * outer(scale, scale)
  unsupported <- which(original == 0 & variances > 0)
  # No correlation is identifiable for a zero working coordinate. An independent
  # diagonal block preserves the reported variance and positive semidefiniteness.
  if (length(unsupported)) diag(sigma)[unsupported] <- variances[unsupported]
  list(Sigma = symmetrize(sigma), scale = scale, unmodeled_variance = unsupported,
       variances = variances, working_marginals_used = which(missing))
}
