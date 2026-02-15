#' Create correlation matrix for the treatment effects
#'
#' @param dataAgD A data.frame with treatment effects and standard errors.
#' @param IMat A list of indicator matrices.
#' @param weightsMat A list of tilting weights matrices.
#'
#' @returns A list of correlation matrices.
#'
CMatFun <- function(dataAgD, IMat, weightsMat) {
  trials <- unique(as.character(dataAgD$trial))
  out <- lapply(trials, function(s) {
    idx <- which(as.character(dataAgD$trial) == s)
    ds <- dataAgD[idx, , drop = FALSE]
    Js <- nrow(ds)
    C <- diag(1, Js)

    if (Js == 1L) {
      return(C)
    }

    I_s <- IMat[[s]]
    w_s <- weightsMat[[s]][, 1]

    prop <- ds$prop
    prop[is.na(prop) | prop <= 0] <- 1
    if ("subgroup" %in% names(ds)) {
      prop[ds$subgroup == "overall"] <- 1
    }
    se <- as.numeric(ds$seTE)

    for (j in seq_len(Js - 1L)) {
      for (k in (j + 1L):Js) {
        gj <- as.character(ds$subgroup[j])
        gk <- as.character(ds$subgroup[k])
        sj <- as.character(ds$stratum[j])
        sk <- as.character(ds$stratum[k])

        corr_jk <- 0

        if (gj == "overall" && gk != "overall") {
          corr_jk <- prop[k] * se[k] / se[j]
        } else if (gk == "overall" && gj != "overall") {
          corr_jk <- prop[j] * se[j] / se[k]
        } else if (gj == gk && sj != sk) {
          corr_jk <- 0
        } else {
          overlap <- mean(w_s * I_s[, j] * I_s[, k])
          denom <- prop[j] * prop[k]
          corr_jk <- if (denom > 0) overlap / denom else 0
        }

        corr_jk <- max(min(corr_jk, 0.999), -0.999)
        C[j, k] <- corr_jk
        C[k, j] <- corr_jk
      }
    }

    as.matrix(Matrix::nearPD(C, corr = TRUE)$mat)
  })
  names(out) <- trials
  out
}


# build reported sampling covariance blocks for treatment effects
tau_vcov_blocks <- function(dataAgD, IMat, weightsMat, V = NULL) {
  trials <- unique(as.character(dataAgD$trial))
  if (!is.null(V)) {
    out <- lapply(seq_along(trials), function(i) {
      s <- trials[i]
      if (!is.null(names(V)) && s %in% names(V)) {
        as.matrix(V[[s]])
      } else {
        as.matrix(V[[i]])
      }
    })
    names(out) <- trials
    return(out)
  }

  C_list <- CMatFun(dataAgD = dataAgD, IMat = IMat, weightsMat = weightsMat)
  out <- lapply(trials, function(s) {
    idx <- which(as.character(dataAgD$trial) == s)
    se <- as.numeric(dataAgD$seTE[idx])
    D <- diag(se, nrow = length(se))
    V_s <- D %*% C_list[[s]] %*% D
    as.matrix(Matrix::nearPD(V_s)$mat)
  })
  names(out) <- trials
  out
}


#' Variance-covariance matrix function
#'
#' @param theta_hat A fitted parameter vector.
#' @param W A weighting matrix.
#' @param Dtheta A derivative matrix for theta.
#' @param Q A matrix of covariates.
#' @param designMat A design matrix.
#' @param tiltingMat A tilting matrix.
#' @param IMat A list of indicator matrices.
#' @param weightsMat A list of tilting weights matrices.
#' @param dataAgD A data.frame with treatment effects and standard errors.
#' @param mu A list of tilting targets for each trial.
#' @param V A list of variance-covariance matrices of treatment effects.
#'
#' @returns The variance-covariance matrix.
#' @export vcovFun
#'
vcovFun <- function(theta_hat, W, Dtheta, Q, designMat, tiltingMat,
                    IMat, weightsMat, dataAgD, mu, V = NULL) {
  Q <- as.data.frame(Q)
  n_q <- nrow(Q)
  J <- nrow(dataAgD)
  eps <- 1e-8
  trials <- unique(as.character(dataAgD$trial))
  n_map <- trial_n_map(dataAgD)

  alpha_q <- alpha_matrix(dataAgD, IMat, weightsMat)
  g_q <- as.numeric(as.matrix(designMat) %*% theta_hat)
  tau_hat <- as.numeric(dataAgD$TE)
  M <- alpha_q * g_q - matrix(tau_hat, nrow = n_q, ncol = J, byrow = TRUE)

  # n^{-1} * Omega first term from Q
  Omega_ninv <- crossprod(M) / (n_q^2)

  Vtau_blocks_report <- tau_vcov_blocks(dataAgD = dataAgD, IMat = IMat, weightsMat = weightsMat, V = V)

  for (s in trials) {
    idx <- which(as.character(dataAgD$trial) == s)
    n_s <- n_map[s]

    h_plus <- as.matrix(tiltingMat)
    p_h <- ncol(h_plus)
    mu_s <- mu[[s]]
    mu_plus <- if (length(mu_s) == p_h - 1L) c(1, mu_s) else mu_s

    w_s <- weightsMat[[s]][, 1]
    alpha_s <- alpha_q[, idx, drop = FALSE]

    # J^m_{eta_s}
    Jm_eta_s_block <- crossprod(alpha_s * g_q, h_plus) / n_q
    Jm_eta_s <- matrix(0, nrow = J, ncol = p_h)
    Jm_eta_s[idx, ] <- Jm_eta_s_block

    # J^{eta_s}_{mu_s}
    Hwh <- crossprod(h_plus, h_plus * w_s) / n_q
    Jeta_mu_s <- -solve(Hwh + diag(eps, p_h))
    Jmu_to_m_s <- Jm_eta_s %*% Jeta_mu_s

    centered_wh <- h_plus * w_s - matrix(mu_plus, nrow = n_q, ncol = p_h, byrow = TRUE)
    Phi_wmu_to_m <- centered_wh %*% t(Jmu_to_m_s)

    # Q-sampling pieces in n^{-1}Omega
    Omega_ninv <- Omega_ninv +
      crossprod(Phi_wmu_to_m) / (n_q^2) +
      2 * crossprod(M, Phi_wmu_to_m) / (n_q^2)

    # Trial-reporting pieces in n^{-1}Omega
    centered_h <- h_plus - matrix(mu_plus, nrow = n_q, ncol = p_h, byrow = TRUE)
    Sigma_mu_s <- crossprod(centered_h, centered_h * w_s) / n_q

    Omega_ninv[idx, idx] <- Omega_ninv[idx, idx] + Vtau_blocks_report[[s]]
    Omega_ninv <- Omega_ninv + (Jmu_to_m_s %*% Sigma_mu_s %*% t(Jmu_to_m_s)) / n_s
  }

  Jtheta_m <- -solve(crossprod(Dtheta, W %*% Dtheta) + diag(eps, ncol(Dtheta))) %*% crossprod(Dtheta, W)
  var_theta <- Jtheta_m %*% Omega_ninv %*% t(Jtheta_m)
  var_theta <- (var_theta + t(var_theta)) / 2
  e <- eigen(var_theta, symmetric = TRUE)
  e$values[e$values < eps] <- eps
  var_theta <- e$vectors %*% diag(e$values, nrow = length(e$values)) %*% t(e$vectors)

  list(
    vcov_theta = var_theta,
    Omega_ninv = Omega_ninv,
    Jtheta_m = Jtheta_m
  )
}
