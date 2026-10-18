# Joint target subgroup / time influence functions (A.1, A.3, A.6.1).
fit_target <- function(target, spec, theta, times, subgroups, variance, shared) {
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
  rows <- list(); derivatives <- influences <- list(); k <- 0L
  for (tt in times) {
    b <- if (shared) spec$basis_q[[as.character(tt)]] else new_basis(spec, target, tt)
    g <- as.vector(b %*% theta)
    for (name in names(groups)) {
      k <- k + 1L; ii <- membership[[name]]; psi <- mean(g[ii])
      rows[[k]] <- data.frame(subgroup = name, Time = tt, estimate = psi, n = sum(ii))
      derivatives[[k]] <- colMeans(b[ii, , drop = FALSE])
      influences[[k]] <- as.numeric(ii) / mean(ii) * (g - psi)
    }
  }
  table <- do.call(rbind, rows); rownames(table) <- NULL
  J_psi_theta <- do.call(rbind, derivatives)
  phi_target <- do.call(cbind, influences)
  names_estimates <- paste(table$subgroup, table$Time, sep = "@")
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
  }
  table$se <- sqrt(pmax(diag(cov), 0))
  table$lower <- table$estimate - stats::qnorm(.975) * table$se
  table$upper <- table$estimate + stats::qnorm(.975) * table$se
  list(estimates = table, vcov = cov, J_psi_theta = J_psi_theta,
       phi_target = phi_target, components = source, shared_q = shared)
}
