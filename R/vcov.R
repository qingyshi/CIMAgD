# Appendix A.6.1: all quantities below are on the influence-covariance scale.
# The parameter may be one common CATE or the stacked trial-specific CATEs used by
# MVLP (A.7); each moment row uses the CATE of its own parameter block.
variance_components <- function(context, fit, Sigma, mu_se, formula_control, control) {
  d <- context$data; Q <- context$Q; moments <- context$moments
  alpha <- moments$alpha; tilting <- context$tilting; spec <- context$basis
  n_q <- nrow(Q); J <- nrow(d); p <- length(fit$theta)
  blocks <- fit$blocks %||% cate_blocks(d, colnames(moments$J_m_theta))
  block_of_row <- integer(J)
  for (k in seq_along(blocks)) block_of_row[blocks[[k]]$rows] <- k
  cate_cache <- list()
  cate_on_q <- function(k, tt) {
    key <- paste(k, tt, sep = "@")
    if (is.null(cate_cache[[key]]))
      cate_cache[[key]] <<- as.vector(spec$basis_q[[tt]] %*% fit$theta[blocks[[k]]$columns])
    cate_cache[[key]]
  }
  phi_q_m <- matrix(0, n_q, J, dimnames = dimnames(alpha))
  for (k in seq_along(blocks)) for (tt in names(spec$basis_q)) {
    idx <- which(block_of_row == k & moments$row_times == tt)
    if (length(idx)) phi_q_m[, idx] <- alpha[, idx, drop = FALSE] *
      sweep(matrix(cate_on_q(k, tt), n_q, length(idx)), 2, fit$m_hat[idx], "-")
  }
  xi_q <- phi_q_m
  trial_covariance <- matrix(0, p, p)
  trial_moment_covariance <- matrix(0, J, J)
  Gamma_s <- J_m_eta <- J_m_mu <- trial_details <- list()
  selection <- list()
  pooled_cache <- list()
  for (s in names(tilting)) {
    idx <- which(d$trial == s); ds <- d[idx, , drop = FALSE]; tilt <- tilting[[s]]
    h <- tilt$h_plus; w <- tilt$weights; mu_plus <- tilt$mu_plus
    n_s <- context$trial_info$n[match(s, context$trial_info$trial)]
    allocation <- context$trial_info$allocation[match(s, context$trial_info$trial)]
    eta_block <- crossprod(phi_q_m[, idx, drop = FALSE], h) / n_q
    J_eta <- matrix(0, J, ncol(h)); J_eta[idx, ] <- eta_block
    J_mu <- J_eta %*% tilt$H_inv
    dimnames(J_eta) <- dimnames(J_mu) <- list(d$row_id, colnames(h))
    xi_q <- xi_q - sweep(h * w, 2, mu_plus, "-") %*% t(J_mu)
    J_m_eta[[s]] <- J_eta; J_m_mu[[s]] <- J_mu
    E_s <- diag(J)[, idx, drop = FALSE]
    D_s <- cbind(-E_s, J_mu)
    coordinate_names <- c(paste0("tau:", ds$row_id), colnames(h))
    nuisance <- working <- anchor <- NULL
    if (!is.null(Sigma) && s %in% names(Sigma)) {
      sigma <- validate_covariance(Sigma[[s]], length(idx) + ncol(h), paste0("Sigma[['", s, "']]"))
      if (!is.null(rownames(sigma))) {
        if (!identical(rownames(sigma), coordinate_names) || !identical(colnames(sigma), coordinate_names))
          stop("Sigma dimnames must match tau row order followed by h+ for trial ", s, ".", call. = FALSE)
      }
      k0 <- length(idx) + 1L
      if (max(abs(sigma[k0, ])) > 1e-10)
        stop("Sigma fixed normalization coordinate must have zero variance/covariance.", call. = FALSE)
      method <- "Supplied complete joint influence covariance"
      selection[[s]] <- data.frame(trial = s, model = "supplied Sigma", trial_specific_issues = NA_character_,
                                   common_issues = NA_character_, reason = "complete joint Sigma supplied")
    } else {
      if (length(unique(ds$Time)) > 1L) stop("Repeated horizons within trial ", s,
        " require a joint Sigma with cross-time covariance (variance=FALSE permits a one-step common-CATE point fit).", call. = FALSE)
      time_s <- moments$row_times[idx[1L]]
      f0 <- NULL
      if (is.null(formula_control)) b0 <- spec$basis_q[[time_s]] else {
        f0 <- if (is.list(formula_control) && !inherits(formula_control, "formula"))
          named_trial_item(formula_control, s, "formula_control") else formula_control
        b0 <- basis_matrix(as_rhs_formula(f0, "formula_control"), Q, "Control formula")
      }
      g <- cate_on_q(block_of_row[idx[1L]], time_s)
      # The common working mean expresses every trial's control reports in the
      # receiving trial's working basis; it is fitted only when needed.
      pooled_rows <- function() {
        key <- if (is.null(f0)) "cate-basis" else paste(deparse(f0), collapse = " ")
        if (is.null(pooled_cache[[key]])) {
          A <- if (is.null(f0)) moments$J_m_theta else crossprod(alpha, b0) / n_q
          pooled_cache[[key]] <<- list(A = A, nu0 = d$nu0, se0 = d$se_nu0, row_id = d$row_id, trial = d$trial)
        }
        pooled_cache[[key]]
      }
      nuisance <- select_working_control(s, b0, alpha[, idx, drop = FALSE], w, g, ds,
                                         context$outcome, control, pooled_rows)
      selection[[s]] <- nuisance$selection$summary
      working <- working_joint_covariance(context$indicators[, idx, drop = FALSE], w, h, nuisance, allocation)
      se_mu <- rep(NA_real_, ncol(h)); se_mu[1L] <- 0
      if (!is.null(mu_se) && s %in% names(mu_se)) {
        provided <- mu_se[[s]]
        original_names <- tilt$original_names[-1L]
        if (!is.numeric(provided) || any(provided < 0, na.rm = TRUE) || any(is.infinite(provided)))
          stop("mu_se must contain nonnegative SEs (or NA).", call. = FALSE)
        if (!is.null(names(provided))) {
          if (!setequal(names(provided), original_names)) stop("mu_se names must match original non-intercept moments.", call. = FALSE)
          provided <- provided[original_names]
        }
        if (length(provided) != length(original_names)) stop("mu_se length mismatch for ", s, ".", call. = FALSE)
        se_mu[match(original_names, colnames(h))] <- provided
      }
      anchor <- anchor_covariance(working$Sigma, n_s * c(ds$seTE, se_mu)^2)
      sigma <- anchor$Sigma
      method <- paste0("A.6.3 working-law reconstruction (", nuisance$model,
                       " working control mean) with marginal variance anchoring")
    }
    dimnames(sigma) <- list(coordinate_names, coordinate_names)
    gamma <- symmetrize(D_s %*% sigma %*% t(D_s))
    Gamma_s[[s]] <- gamma
    trial_moment_covariance <- trial_moment_covariance + gamma / n_s
    vc_s <- fit$J_theta_m %*% gamma %*% t(fit$J_theta_m) / n_s
    trial_covariance <- trial_covariance + vc_s
    trial_details[[s]] <- list(Sigma = sigma, D_s = D_s, n = n_s,
      method = method, control = nuisance, working = working, anchoring = anchor, vcov_theta = vc_s)
  }
  xi_q <- center_columns(xi_q)
  Gamma_q <- crossprod(xi_q) / n_q
  theta_if_q <- xi_q %*% t(fit$J_theta_m)
  q_covariance <- crossprod(theta_if_q) / n_q^2
  vcov_theta <- symmetrize(trial_covariance + q_covariance)
  dimnames(vcov_theta) <- list(names(fit$theta), names(fit$theta))
  working_control <- do.call(rbind, selection); rownames(working_control) <- NULL
  list(vcov_theta = vcov_theta, q_covariance = q_covariance,
       trial_covariance = symmetrize(trial_covariance), xi_q = xi_q,
       theta_if_q = theta_if_q, phi_q_m = phi_q_m, Gamma_q = Gamma_q, Gamma_s = Gamma_s,
       J_m_eta = J_m_eta, J_m_mu = J_m_mu, trials = trial_details,
       working_control = working_control,
       trial_moment_covariance = trial_moment_covariance,
       moment_covariance = symmetrize(Gamma_q / n_q + trial_moment_covariance))
}
