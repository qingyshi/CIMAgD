make_basis_spec <- function(formula, Q, times, time_indexed) {
  f <- as_rhs_formula(formula, "formula_cate")
  if (!time_indexed && "Time" %in% all.vars(f)) stop("Set time_indexed=TRUE to use Time in formula_cate.", call. = FALSE)
  stacked <- do.call(rbind, lapply(times, function(tt) {
    d <- Q
    if (time_indexed) d$Time <- tt
    d
  }))
  mf <- stats::model.frame(f, stacked, na.action = stats::na.pass)
  b <- stats::model.matrix(f, mf)
  if (nrow(b) != nrow(stacked) || any(!is.finite(b)) || !ncol(b))
    stop("formula_cate must produce a finite nonempty basis.", call. = FALSE)
  xlev <- lapply(mf[vapply(mf, is.factor, logical(1))], levels)
  blocks <- lapply(seq_along(times), function(k) b[(k - 1L) * nrow(Q) + seq_len(nrow(Q)), , drop = FALSE])
  names(blocks) <- as.character(times)
  list(formula = f, xlev = xlev, contrasts = attr(b, "contrasts"), columns = colnames(b),
       basis_q = blocks, times = times, time_indexed = time_indexed)
}

new_basis <- function(spec, data, time) {
  if (spec$time_indexed) data$Time <- time
  mf <- stats::model.frame(spec$formula, data, xlev = spec$xlev, na.action = stats::na.pass)
  b <- stats::model.matrix(spec$formula, mf, contrasts.arg = spec$contrasts)
  if (!identical(colnames(b), spec$columns) || nrow(b) != nrow(data) || any(!is.finite(b)))
    stop("New covariates do not match the fitted CATE basis.", call. = FALSE)
  b
}

build_moments <- function(d, indicators, tilting, spec) {
  n_q <- nrow(indicators); J <- nrow(d)
  alpha <- matrix(0, n_q, J, dimnames = dimnames(indicators))
  c_sj <- numeric(J)
  for (s in names(tilting)) {
    idx <- which(d$trial == s)
    wi <- indicators[, idx, drop = FALSE] * tilting[[s]]$weights
    c_sj[idx] <- colMeans(wi)
    if (any(c_sj[idx] <= 1e-12)) stop("Empty or negligible tilted Q support for ",
      paste(d$row_id[idx[c_sj[idx] <= 1e-12]], collapse = ", "), call. = FALSE)
    alpha[, idx] <- sweep(wi, 2, c_sj[idx], "/")
  }
  row_times <- if (spec$time_indexed) as.character(d$Time) else rep("1", J)
  J_m_theta <- matrix(0, J, length(spec$columns), dimnames = list(d$row_id, spec$columns))
  for (tt in names(spec$basis_q)) {
    idx <- which(row_times == tt)
    if (length(idx)) J_m_theta[idx, ] <- crossprod(alpha[, idx, drop = FALSE], spec$basis_q[[tt]]) / n_q
  }
  list(alpha = alpha, c_sj = c_sj, J_m_theta = J_m_theta, row_times = row_times)
}

fit_gmm <- function(d, moments, trial_info, W, control, weight_psd = FALSE) {
  J <- nrow(d); D <- moments$J_m_theta
  diagnostics <- jacobian_report(D, "J_m_theta", control, scaled_svd(D, control$rank_tol)$rank)
  weight_note <- "User-supplied positive-definite W."
  if (is.null(W)) {
    n_s <- trial_info$n[match(d$trial, trial_info$trial)]
    variances <- n_s * d$seTE^2
    if (any(variances == 0)) {
      W <- diag(J); weight_note <- "Identity W: at least one reported contrast SE is zero."
    } else {
      W <- diag(1 / variances, J); weight_note <- "Diagonal inverse n_s * seTE^2."
    }
  }
  W <- validate_covariance(W, J, "W", positive = !weight_psd)
  if (weight_psd) {
    e <- eigen(W, symmetric = TRUE)
    L <- sqrt(pmax(e$values, 0)) * t(e$vectors)
  } else L <- chol(W)
  a <- L %*% D
  solution <- least_squares_operator(a, "Weighted J_m_theta", control)
  theta <- as.vector(solution$inverse %*% L %*% d$TE)
  names(theta) <- colnames(D)
  J_theta_m <- -solution$inverse %*% L
  dimnames(J_theta_m) <- list(colnames(D), d$row_id)
  m_hat <- as.vector(D %*% theta)
  residual <- m_hat - d$TE
  list(theta = theta, W = W, weight_note = weight_note, m_hat = m_hat,
       residual = residual, foc = as.vector(crossprod(D, W %*% residual)),
       objective = as.numeric(crossprod(residual, W %*% residual)),
       J_theta_m = J_theta_m, diagnostics = diagnostics,
       regularization = solution[c("method", "used", "rank", "ridge", "scales")])
}
