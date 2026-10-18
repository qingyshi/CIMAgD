# A.6.2 working nuisance regressions only; no A.7 arm-mean transport estimator.
fit_working_control <- function(b, alpha, w, g, reports, outcome, control) {
  A <- crossprod(alpha, b) / nrow(b)
  se <- reports$se_nu0
  reference <- if (any(se > 0)) min(se[se > 0]) else 1
  precision <- 1 / ifelse(se > 0, se^2, reference^2)
  a <- A * sqrt(precision)
  y <- reports$nu0 * sqrt(precision)
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
  transform <- matrix(0, ncol(a), sum(keep))
  transform[selected, ] <- sweep(sweep(z$v[, keep, drop = FALSE], 2, z$d[keep], "/"),
                                 1, column_scale, "/")
  u <- z$u[, keep, drop = FALSE]
  unconstrained <- as.vector(crossprod(u, y))
  b_reduced <- b %*% transform
  constrained <- FALSE
  if (outcome == "binary") {
    if (any(abs(g) > 1 + 1e-8)) stop("Fitted binary CATE leaves [-1, 1] on Q; no compatible working arm risks. Simplify the CATE, supply Sigma, or use variance=FALSE.", call. = FALSE)
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
  list(beta = beta, gamma0 = gamma0, gamma1 = gamma1, v0 = v0, v1 = v1,
       aggregate_means = as.vector(A %*% beta), residual = as.vector(A %*% beta) - reports$nu0,
       rank = sum(keep), basis_dimension = ncol(b),
       retained_columns = colnames(b)[selected], dropped_columns = colnames(b)[-selected],
       identification = "Independent columns in formula order (intercept retained)",
       constrained = constrained, zero_se_weight_fallback = any(se == 0),
       residual_model = residual_model, residual_diagnostics = residual_diagnostics)
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
