# Rank decisions use dimensionless coordinates; the public MP solution minimizes
# Euclidean coefficient norm in the original parameterization.
scaled_svd <- function(a, tol) {
  scale <- sqrt(colSums(a^2)); scale[scale == 0] <- 1
  normalized <- sweep(a, 2, scale, "/")
  z <- svd(normalized, nv = ncol(a))
  keep <- z$d > max(z$d, 0) * tol
  list(scale = scale, normalized = normalized, svd = z, keep = keep, rank = sum(keep))
}

mp_operator <- function(a, decomposition) {
  z <- decomposition$svd; keep <- decomposition$keep
  out <- matrix(0, ncol(a), nrow(a))
  if (any(keep)) out <- sweep(z$v[, which(keep), drop = FALSE] %*%
    (t(z$u[, keep, drop = FALSE]) / z$d[keep]), 1, decomposition$scale, "/")
  if (decomposition$rank < ncol(a)) {
    null <- sweep(z$v[, seq.int(decomposition$rank + 1L, ncol(a)), drop = FALSE],
                  1, decomposition$scale, "/")
    null <- qr.Q(qr(null, LAPACK = TRUE))
    out <- out - null %*% crossprod(null, out)
  }
  out
}

jacobian_report <- function(a, label, control, rank) {
  method <- control$singularity_method %||% "none"
  enabled <- control$jacobian_diagnostics %||% TRUE
  if (rank < ncol(a) && method == "none") {
    report <- diagnose_jacobian(a, control$rank_tol, control$condition_limit)
    report$rank <- rank; report$full_rank <- FALSE; report$status <- "singular"
    stop(structure(list(message = paste0(label, " is singular: rank ", rank, " of ", ncol(a),
      ". Remove redundant terms or choose an explicit singularity method."), call = NULL,
      matrix_name = label, diagnostics = report), class = c("cima_singularity", "error", "condition")))
  }
  if (!enabled) return(NULL)
  report <- diagnose_jacobian(a, control$rank_tol, control$condition_limit)
  if (report$status == "ill_conditioned") warning(label, " is ill-conditioned (scaled condition ",
    signif(report$condition, 4), ").", call. = FALSE)
  report
}

least_squares_operator <- function(a, label, control) {
  z <- scaled_svd(a, control$rank_tol)
  method <- control$singularity_method %||% "none"
  report <- jacobian_report(a, label, control, z$rank)
  used <- z$rank < ncol(a)
  lambda <- if (used && method == "ridge") control$ridge_strength else 0
  if (lambda > 0) {
    inverse <- solve(crossprod(z$normalized) + diag(lambda, ncol(a)), t(z$normalized))
    inverse <- sweep(inverse, 1, z$scale, "/")
  } else inverse <- mp_operator(a, z)
  list(inverse = inverse, diagnostics = report, rank = z$rank, method = method,
       used = used, ridge = lambda, scales = z$scale)
}

# A congruence-scaled ridge H + lambda*diag(diag(H)) is unchanged by
# changing units of h. Zero columns get unit scale and have no fitted effect.
symmetric_inverse <- function(a, label, control, report = TRUE) {
  scale <- sqrt(pmax(diag(a), 0)); scale[scale == 0] <- 1
  correlation <- symmetrize(a / outer(scale, scale))
  eig <- eigen(correlation, symmetric = TRUE)
  rank <- sum(eig$values > max(eig$values, 0) * control$rank_tol)
  method <- control$singularity_method %||% "none"
  if (rank < ncol(a) && method == "none")
    jacobian_report(correlation, label, control, rank)
  dg <- if (report) jacobian_report(correlation, label, control, rank) else NULL
  lambda <- if (rank < ncol(a) && method == "ridge") control$ridge_strength else 0
  if (lambda > 0) inverse <- solve(correlation + diag(lambda, ncol(a))) / outer(scale, scale)
  else if (rank == ncol(a)) inverse <- solve(correlation) / outer(scale, scale)
  else inverse <- mp_operator(a, scaled_svd(a, control$rank_tol))
  list(inverse = symmetrize(inverse), diagnostics = dg, rank = rank,
       method = method, used = rank < ncol(a), ridge = lambda, scales = scale)
}

optimal_weight <- function(covariance, jacobian, control) {
  covariance <- validate_covariance(covariance, nrow(jacobian), "Moment sampling covariance")
  sd <- sqrt(pmax(diag(covariance), 0)); sd[sd == 0] <- 1
  correlation <- symmetrize(covariance / outer(sd, sd))
  e <- eigen(correlation, symmetric = TRUE)
  keep <- e$values > max(e$values, 0) * control$rank_tol
  if (!any(keep)) stop("Two-step GMM needs a nonzero moment covariance; use one-step GMM for deterministic moments.", call. = FALSE)
  lambda <- 0
  if (!all(keep)) {
    D <- sweep(jacobian, 1, sd, "/")
    columns <- sqrt(colSums(D^2)); columns[columns == 0] <- 1
    informative_null <- max(abs(crossprod(e$vectors[, !keep, drop = FALSE],
      sweep(D, 2, columns, "/")))) > sqrt(control$rank_tol)
    if (informative_null && control$singularity_method != "ridge")
      stop("Moment covariance has informative zero-variance directions. Use one-step GMM or singularity_method='ridge'.", call. = FALSE)
    if (informative_null) lambda <- control$ridge_strength
  }
  inverse <- if (lambda > 0) solve(correlation + diag(lambda, nrow(correlation))) else
    tcrossprod(sweep(e$vectors[, keep, drop = FALSE], 2, sqrt(e$values[keep]), "/"))
  W <- symmetrize(inverse / outer(sd, sd))
  list(W = W, covariance = covariance, rank = sum(keep), ridge = lambda,
       method = if (lambda > 0) "scaled ridge" else if (all(keep)) "inverse" else "inverse on covariance range")
}
