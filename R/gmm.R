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

# CATE parameter blocks. The common method uses one CATE for every moment. MVLP
# (A.7) uses one CATE per trial: stacking the trial parameters with a
# block-diagonal moment Jacobian and weight reproduces the separate trial-specific
# GMM fits exactly, while the shared Q/target influences remain jointly available.
cate_blocks <- function(d, columns, method = "common") {
  p <- length(columns)
  if (method == "common")
    return(list(common = list(name = "common", rows = seq_len(nrow(d)),
                              columns = seq_len(p), labels = columns)))
  trials <- unique(d$trial)
  stats::setNames(lapply(seq_along(trials), function(k) list(name = trials[k],
    rows = which(d$trial == trials[k]), columns = (k - 1L) * p + seq_len(p),
    labels = paste0(trials[k], ":", columns))), trials)
}

is_trial_blocks <- function(blocks) !identical(names(blocks), "common")

# Moment Jacobian with respect to the stacked parameter (J x P).
cate_design <- function(B, blocks) {
  labels <- unlist(lapply(blocks, `[[`, "labels"), use.names = FALSE)
  D <- matrix(0, nrow(B), length(labels), dimnames = list(rownames(B), labels))
  for (blk in blocks) D[blk$rows, blk$columns] <- B[blk$rows, , drop = FALSE]
  D
}

default_weight <- function(d, trial_info, blocks) {
  J <- nrow(d)
  variances <- trial_info$n[match(d$trial, trial_info$trial)] * d$seTE^2
  W <- matrix(0, J, J)
  notes <- character()
  # The identity fallback is decided within each separately fitted CATE block.
  for (blk in blocks) {
    r <- blk$rows
    if (any(variances[r] == 0)) {
      W[r, r] <- diag(length(r)); notes[blk$name] <- "Identity W: at least one reported contrast SE is zero."
    } else {
      W[r, r] <- diag(1 / variances[r], length(r)); notes[blk$name] <- "Diagonal inverse n_s * seTE^2."
    }
  }
  note <- if (is_trial_blocks(blocks)) paste0(names(notes), ": ", notes, collapse = " ") else unname(notes)
  list(W = W, note = note)
}

fit_gmm <- function(d, moments, trial_info, W, control, weight_psd = FALSE, blocks = NULL) {
  J <- nrow(d); B <- moments$J_m_theta
  blocks <- blocks %||% cate_blocks(d, colnames(B))
  by_trial <- is_trial_blocks(blocks)
  D <- cate_design(B, blocks)
  if (is.null(W)) {
    default <- default_weight(d, trial_info, blocks)
    W <- default$W; weight_note <- default$note
  } else {
    W <- validate_covariance(W, J, "W", positive = !weight_psd)
    weight_note <- if (by_trial) "User-supplied W; MVLP uses its trial-diagonal blocks." else "User-supplied positive-definite W."
    if (by_trial) {
      diagonal <- matrix(0, J, J)
      for (blk in blocks) diagonal[blk$rows, blk$rows] <- W[blk$rows, blk$rows]
      W <- diagonal
    }
  }
  theta <- stats::setNames(numeric(ncol(D)), colnames(D))
  J_theta_m <- matrix(0, ncol(D), J, dimnames = list(colnames(D), d$row_id))
  diagnostics <- regularization <- list()
  identification <- vector("list", length(blocks))
  for (k in seq_along(blocks)) {
    blk <- blocks[[k]]; r <- blk$rows
    Db <- B[r, , drop = FALSE]
    label <- if (by_trial) paste0("J_m_theta for trial ", blk$name) else "J_m_theta"
    rank <- scaled_svd(Db, control$rank_tol)$rank
    report <- jacobian_report(Db, label, control, rank)
    diagnostics[blk$name] <- list(report)
    Wb <- W[r, r, drop = FALSE]
    if (weight_psd) {
      e <- eigen(Wb, symmetric = TRUE)
      L <- sqrt(pmax(e$values, 0)) * t(e$vectors)
    } else L <- chol(Wb)
    solution <- least_squares_operator(L %*% Db, paste("Weighted", label), control)
    theta[blk$columns] <- as.vector(solution$inverse %*% L %*% d$TE[r])
    J_theta_m[blk$columns, r] <- -solution$inverse %*% L
    regularization[[blk$name]] <- solution[c("method", "used", "rank", "ridge", "scales")]
    identification[[k]] <- data.frame(block = blk$name, moments = length(r), parameters = ncol(Db),
      rank = rank, condition = if (is.null(report)) NA_real_ else report$condition,
      status = if (is.null(report)) if (rank < ncol(Db)) "singular" else "not assessed" else report$status)
  }
  m_hat <- as.vector(D %*% theta)
  residual <- m_hat - d$TE
  identification <- do.call(rbind, identification); rownames(identification) <- NULL
  if (by_trial) names(identification)[1L] <- "trial"
  list(theta = theta, W = W, weight_note = weight_note, m_hat = m_hat,
       residual = residual, foc = as.vector(crossprod(D, W %*% residual)),
       objective = as.numeric(crossprod(residual, W %*% residual)),
       J_theta_m = J_theta_m, blocks = blocks, design = D, identification = identification,
       diagnostics = if (by_trial) { if (all(vapply(diagnostics, is.null, logical(1)))) NULL else diagnostics } else diagnostics$common,
       regularization = if (by_trial) list(method = control$singularity_method %||% "none",
         used = any(vapply(regularization, `[[`, logical(1), "used")), by_trial = regularization) else regularization$common)
}
