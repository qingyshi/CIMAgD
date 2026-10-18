normalize_mu <- function(mu, h, trial) {
  if (!is.numeric(mu) || any(!is.finite(mu))) stop("mu must be finite numeric for trial ", trial, ".", call. = FALSE)
  nm <- colnames(h)
  if (!is.null(names(mu))) {
    expected <- if ("(Intercept)" %in% names(mu)) nm else nm[-1L]
    if (anyDuplicated(names(mu)) || !setequal(names(mu), expected))
      stop("Named mu does not match tilting basis for trial ", trial, ": ", paste(expected, collapse = ", "), call. = FALSE)
    mu <- mu[expected]
  }
  if (length(mu) == ncol(h) - 1L) mu <- c(1, mu)
  if (length(mu) != ncol(h) || abs(mu[1L] - 1) > 1e-12)
    stop("mu must match h (or h+ with normalization 1) for trial ", trial, ".", call. = FALSE)
  stats::setNames(as.numeric(mu), nm)
}

fit_tilt <- function(h, mu_plus, trial, control) {
  if (colnames(h)[1L] != "(Intercept)" || any(h[, 1L] != 1))
    stop("formula_tilting must include an intercept.", call. = FALSE)
  p <- ncol(h); n_q <- nrow(h)
  center <- colMeans(h); center[1L] <- 0
  scale <- sqrt(colMeans(sweep(h, 2, center, "-")^2))
  scale[scale == 0] <- 1
  transform <- diag(1 / scale, p)
  if (p > 1L) transform[1L, -1L] <- -center[-1L] / scale[-1L]
  z <- h %*% transform
  target <- as.vector(mu_plus %*% transform)
  colnames(z) <- colnames(h)
  z_decomposition <- scaled_svd(z, control$rank_tol)
  if (z_decomposition$rank < p && (control$singularity_method %||% "none") == "none")
    require_identified(z, paste0("Tilting basis for trial ", trial, " (H_s cannot be invertible)"),
                      control$rank_tol, control$condition_limit)
  if (z_decomposition$rank < p) {
    null <- z_decomposition$svd$v[, seq.int(z_decomposition$rank + 1L, p), drop = FALSE]
    discrepancy <- as.vector((target / z_decomposition$scale) %*% null)
    if (max(abs(discrepancy)) > 10 * control$tilt_tol)
      stop("Incompatible redundant tilting targets for trial ", trial, ".", call. = FALSE)
  }
  if (p > 1L) {
    for (k in 2:p) if (diff(range(z[, k])) > control$tilt_tol &&
                      (target[k] <= min(z[, k]) + control$tilt_tol ||
                      target[k] >= max(z[, k]) - control$tilt_tol))
      stop("No finite interior tilting solution for trial ", trial, ": target for ",
           colnames(h)[k], " is on/outside Q support.", call. = FALSE)
  }
  evaluate <- function(lambda) {
    lp <- if (p > 1L) as.vector(z[, -1L, drop = FALSE] %*% lambda) else rep(0, n_q)
    a <- max(lp); exp_lp <- exp(lp - a); prob <- exp_lp / sum(exp_lp)
    log_normalizer <- a + log(mean(exp_lp))
    mean_z <- colSums(z * prob)
    list(objective = log_normalizer - sum(target[-1L] * lambda), prob = prob,
         gradient = mean_z[-1L] - target[-1L], mean_z = mean_z,
         log_normalizer = log_normalizer)
  }
  lambda <- rep(0, p - 1L); current <- evaluate(lambda); iterations <- 0L
  if (p > 1L) for (iter in seq_len(control$max_iter)) {
    iterations <- iter
    if (max(abs(current$gradient)) <= control$tilt_tol) break
    zc <- sweep(z[, -1L, drop = FALSE], 2, current$mean_z[-1L], "-")
    curvature <- crossprod(zc, zc * current$prob)
    step <- tryCatch(as.vector(symmetric_inverse(curvature, "Tilting Newton curvature",
      control, report = FALSE)$inverse %*% current$gradient), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) stop("Tilting curvature collapsed for trial ", trial,
      "; check joint feasibility and Q support.", call. = FALSE)
    rate <- 1; accepted <- FALSE
    for (ls in seq_len(50L)) {
      proposed <- evaluate(lambda - rate * step)
      if (is.finite(proposed$objective) &&
          proposed$objective <= current$objective - 1e-4 * rate * sum(current$gradient * step) + 1e-14) {
        accepted <- TRUE; break
      }
      rate <- rate / 2
    }
    if (!accepted) break
    lambda <- lambda - rate * step; current <- proposed
  }
  error <- if (p > 1L) max(abs(current$gradient)) else 0
  if (!is.finite(error) || error > control$tilt_tol * 10)
    stop("Tilting did not converge for trial ", trial, "; standardized moment error = ", signif(error, 4),
         ". Check joint feasibility or increase max_iter.", call. = FALSE)
  w <- n_q * current$prob
  eta <- as.vector(transform %*% c(-current$log_normalizer, lambda))
  names(eta) <- colnames(h)
  H_s <- crossprod(h, h * w) / n_q
  H_solve <- symmetric_inverse(H_s, paste0("H_s for trial ", trial), control)
  diag_H <- H_solve$diagnostics
  if (!is.null(diag_H)) {
    diag_H$coordinate_system <- "diagonal-normalized h+ Gram matrix"
    diag_H$original <- diagnose_jacobian(H_s, control$rank_tol, control$condition_limit)
  }
  H_inv <- H_solve$inverse
  dimnames(H_inv) <- dimnames(H_s)
  if (z_decomposition$rank < p && control$singularity_method == "minimum-norm")
    eta <- stats::setNames(as.vector(mp_operator(h, scaled_svd(h, control$rank_tol)) %*% log(w)), colnames(h))
  list(eta = eta, weights = w, h_plus = h, mu_plus = mu_plus, H_s = H_s, H_inv = H_inv,
       diagnostics = diag_H, regularization = H_solve[c("method", "used", "rank", "ridge", "scales")], iterations = iterations, moment_error = colMeans(h * w) - mu_plus,
       standardized_error = error, effective_n = sum(w)^2 / sum(w^2), max_weight = max(w))
}

fit_tilting <- function(Q, d, indicators, mu, formula_tilting, trial_info, match_proportions, control) {
  trials <- unique(d$trial)
  out <- stats::setNames(vector("list", length(trials)), trials)
  common <- !(is.list(formula_tilting) && !inherits(formula_tilting, "formula"))
  common_formula <- if (common) as_rhs_formula(formula_tilting, "formula_tilting") else NULL
  common_basis <- if (common) basis_matrix(common_formula, Q, "Tilting formula") else NULL
  for (s in trials) {
    f <- if (is.list(formula_tilting) && !inherits(formula_tilting, "formula"))
      named_trial_item(formula_tilting, s, "formula_tilting") else formula_tilting
    f <- as_rhs_formula(f, "formula_tilting")
    if ("Time" %in% all.vars(f)) stop("Tilting uses baseline covariates; Time is not allowed.", call. = FALSE)
    h <- if (common) common_basis else basis_matrix(f, Q, "Tilting formula")
    if (colnames(h)[1L] != "(Intercept)") stop("formula_tilting must include an intercept.", call. = FALSE)
    mu_plus <- normalize_mu(named_trial_item(mu, s, "mu"), h, s)
    original_names <- colnames(h)
    if (match_proportions) {
      idx <- which(d$trial == s)
      n_s <- trial_info$n[match(s, trial_info$trial)]
      for (j in idx) {
        v <- indicators[, j]; p_j <- d$n[j] / n_s
        co <- minimum_norm(h, v, control$rank_tol)
        if (max(abs(v - h %*% co)) < control$rank_tol * 10) {
          if (abs(sum(co * mu_plus) - p_j) > 1e-6)
            stop("Reported subgroup proportion conflicts with mu or another proportion: ", d$row_id[j], call. = FALSE)
        } else {
          nm <- paste0("prop:", d$subgroup[j], ":", d$stratum[j])
          h <- cbind(h, v); colnames(h)[ncol(h)] <- nm
          mu_plus <- c(mu_plus, stats::setNames(p_j, nm))
        }
      }
    }
    out[[s]] <- fit_tilt(h, mu_plus, s, control)
    out[[s]]$formula <- f
    out[[s]]$original_names <- original_names
  }
  out
}
