# Joint target subgroup / time influence functions (A.1, A.3, A.6.1). With
# trial-specific CATE blocks (MVLP), every trial is marginalized over the same
# target sample, so the returned covariance contains the cross-trial terms of A.7.
fit_target <- function(target, spec, fit, times, subgroups, variance, shared) {
  groups <- c(list(overall = function(d) rep(TRUE, nrow(d))), subgroups)
  if (is.null(names(groups)) || any(!nzchar(names(groups))) || anyDuplicated(names(groups)))
    stop("target_subgroups must have unique names other than overall.", call. = FALSE)
  membership <- lapply(groups, function(f) {
    if (!is.function(f)) stop("Each target subgroup must be a function of the target data frame.", call. = FALSE)
    v <- f(target)
    if (!is.logical(v) || length(v) != nrow(target) || anyNA(v) || !any(v))
      stop("Each target subgroup must return nonmissing logical membership with at least one member.", call. = FALSE)
    v
  })
  blocks <- fit$blocks %||% cate_blocks(data.frame(trial = character()), names(fit$theta))
  by_trial <- is_trial_blocks(blocks)
  P <- length(fit$theta)
  bases <- lapply(times, function(tt) if (shared) spec$basis_q[[as.character(tt)]] else new_basis(spec, target, tt))
  names(bases) <- as.character(times)
  rows <- list(); derivatives <- influences <- list(); k <- 0L
  for (blk in blocks) for (tt in times) {
    b <- bases[[as.character(tt)]]
    g <- as.vector(b %*% fit$theta[blk$columns])
    for (name in names(groups)) {
      k <- k + 1L; ii <- membership[[name]]; psi <- mean(g[ii])
      rows[[k]] <- data.frame(trial = blk$name, subgroup = name, Time = tt, estimate = psi, n = sum(ii))
      derivative <- numeric(P); derivative[blk$columns] <- colMeans(b[ii, , drop = FALSE])
      derivatives[[k]] <- derivative
      influences[[k]] <- as.numeric(ii) / mean(ii) * (g - psi)
    }
  }
  table <- do.call(rbind, rows); rownames(table) <- NULL
  if (!by_trial) table$trial <- NULL
  J_psi_theta <- do.call(rbind, derivatives)
  colnames(J_psi_theta) <- names(fit$theta)
  phi_target <- do.call(cbind, influences)
  names_estimates <- paste(table$subgroup, table$Time, sep = "@")
  if (by_trial) names_estimates <- paste0(table$trial, ":", names_estimates)
  rownames(J_psi_theta) <- colnames(phi_target) <- names_estimates
  cov <- matrix(NA_real_, nrow(table), nrow(table), dimnames = list(names_estimates, names_estimates))
  source <- NULL
  if (!is.null(variance)) {
    trial <- J_psi_theta %*% variance$trial_covariance %*% t(J_psi_theta)
    if (shared) {
      zeta <- phi_target + variance$theta_if_q %*% t(J_psi_theta)
      source_cov <- crossprod(center_columns(zeta)) / nrow(target)^2
      cov <- symmetrize(source_cov + trial)
      source <- list(shared_target_q = source_cov, trials = trial)
    } else {
      target_cov <- crossprod(center_columns(phi_target)) / nrow(target)^2
      q_cov <- J_psi_theta %*% variance$q_covariance %*% t(J_psi_theta)
      cov <- symmetrize(target_cov + q_cov + trial)
      source <- list(target = target_cov, Q = q_cov, trials = trial)
    }
    dimnames(cov) <- list(names_estimates, names_estimates)
    source <- lapply(source, function(x) { dimnames(x) <- dimnames(cov); x })
  }
  table <- add_wald(table, cov)
  list(estimates = table, vcov = cov, J_psi_theta = J_psi_theta,
       phi_target = phi_target, components = source, shared_q = shared)
}

add_wald <- function(table, cov) {
  table$se <- sqrt(pmax(diag(cov), 0))
  table$lower <- table$estimate - stats::qnorm(.975) * table$se
  table$upper <- table$estimate + stats::qnorm(.975) * table$se
  table
}

# A.7: pool the trial-specific transported estimates of each target estimand
# (overall or subgroup, at each horizon) with minimum-variance weights computed
# from that estimand's trial-by-trial covariance block. Weight estimation has no
# first-order effect, so pooled covariances are L' C L for the fixed pooling map L.
pool_target <- function(trial_fit, trials, control) {
  est <- trial_fit$estimates
  C <- trial_fit$vcov
  if (anyNA(C)) stop("MVLP pooling weights require the covariance of the trial-specific transported estimates.", call. = FALSE)
  key <- paste(est$subgroup, est$Time, sep = "@")
  keys <- unique(key)
  L <- matrix(0, nrow(est), length(keys), dimnames = list(rownames(C), keys))
  weights <- matrix(NA_real_, length(trials), length(keys), dimnames = list(trials, keys))
  status <- stats::setNames(character(length(keys)), keys)
  regularized <- FALSE
  blocks <- list()
  for (k in seq_along(keys)) {
    idx <- which(key == keys[k])
    label <- if (length(keys) == 1L) "the overall target effect" else keys[k]
    mv <- minimum_variance_weights(C[idx, idx, drop = FALSE], control, label)
    L[idx, k] <- mv$weights
    weights[match(est$trial[idx], trials), k] <- mv$weights
    status[k] <- mv$status
    regularized <- regularized || mv$regularized
    blocks[[keys[k]]] <- C[idx, idx, drop = FALSE]
  }
  first <- match(keys, key)
  table <- data.frame(subgroup = est$subgroup[first], Time = est$Time[first],
                      estimate = as.vector(crossprod(L, est$estimate)), n = est$n[first])
  vcov <- symmetrize(crossprod(L, C %*% L))
  dimnames(vcov) <- list(keys, keys)
  table <- add_wald(table, vcov)
  components <- if (is.null(trial_fit$components)) NULL else
    lapply(trial_fit$components, function(x) { out <- symmetrize(crossprod(L, x %*% L)); dimnames(out) <- dimnames(vcov); out })
  list(estimates = table, vcov = vcov, map = L, weights = weights, status = status,
       regularized = regularized, covariance_blocks = blocks, components = components,
       J_psi_theta = crossprod(L, trial_fit$J_psi_theta), phi_target = trial_fit$phi_target %*% L)
}
