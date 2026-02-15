# internal helper: trial sample size per trial
trial_n_map <- function(dataAgD) {
  trials <- unique(as.character(dataAgD$trial))
  out <- stats::setNames(rep(1, length(trials)), trials)

  if (!("n" %in% names(dataAgD))) {
    return(out)
  }

  for (s in trials) {
    ds <- dataAgD[dataAgD$trial == s, , drop = FALSE]
    n_s <- NA_real_
    if ("subgroup" %in% names(ds)) {
      n_s <- ds$n[which(ds$subgroup == "overall")[1]]
    }
    if (is.na(n_s) || !is.finite(n_s) || n_s <= 0) {
      n_s <- suppressWarnings(max(ds$n, na.rm = TRUE))
    }
    if (!is.finite(n_s) || n_s <= 0) {
      n_s <- 1
    }
    out[s] <- n_s
  }
  out
}


# internal helper: build alpha(x, eta) matrix
alpha_matrix <- function(dataAgD, IMat, weightsMat) {
  trials <- unique(as.character(dataAgD$trial))
  J <- nrow(dataAgD)
  n_q <- nrow(IMat[[1]])
  A <- matrix(0, nrow = n_q, ncol = J)

  for (s in trials) {
    idx <- which(as.character(dataAgD$trial) == s)
    I_s <- IMat[[s]]
    w_s <- weightsMat[[s]]
    prop_s <- dataAgD$prop[idx]
    prop_s[is.na(prop_s) | prop_s <= 0] <- 1
    A[, idx] <- sweep(I_s * w_s, 2, prop_s, "/")
  }
  A
}


# internal helper: sample moment vector
moment_vector <- function(theta, designMat, IMat, weightsMat, dataAgD) {
  g <- as.numeric(designMat %*% theta)
  alpha <- alpha_matrix(dataAgD = dataAgD, IMat = IMat, weightsMat = weightsMat)
  colMeans(alpha * g)
}


# internal helper: GMM residual moment vector
moment_residual <- function(theta, designMat, IMat, weightsMat, dataAgD) {
  moment_vector(theta, designMat, IMat, weightsMat, dataAgD) - as.numeric(dataAgD$TE)
}


# internal helper: Jacobian matrix D for theta
DthetaMat <- function(designMat, IMat, weightsMat, dataAgD) {
  alpha <- alpha_matrix(dataAgD = dataAgD, IMat = IMat, weightsMat = weightsMat)
  crossprod(alpha, designMat) / nrow(designMat)
}


#' Weighting matrix W
#'
#' @param dataAgD A data.frame with treatment effects and standard errors.
#' @param V A list of variance-covariance matrices.
#'
#' @returns The weighting matrix.
#'
WMat <- function(dataAgD, V = NULL) {
  eps <- 1e-8
  trials <- unique(as.character(dataAgD$trial))
  n_map <- trial_n_map(dataAgD)
  J <- nrow(dataAgD)
  W <- matrix(0, nrow = J, ncol = J)

  for (i in seq_along(trials)) {
    s <- trials[i]
    idx <- which(as.character(dataAgD$trial) == s)
    n_s <- n_map[s]

    if (!is.null(V)) {
      Vs <- if (!is.null(names(V)) && s %in% names(V)) V[[s]] else V[[i]]
      V_report <- as.matrix(Vs)
    } else {
      se <- as.numeric(dataAgD$seTE[idx])
      V_report <- diag(se^2, nrow = length(se))
    }

    # Keep user row order from dataAgD by inserting trial block at its native indices.
    W[idx, idx] <- solve(n_s * V_report + diag(eps, nrow(V_report)))
  }

  W
}


#' Generalized method of moments
#'
#' @param theta A parameter vector.
#' @param W A weighting matrix.
#' @param designMat A design matrix.
#' @param IMat A list of indicator matrices.
#' @param weightsMat A list of tilting weights matrices.
#' @param dataAgD A data.frame with treatment effects and standard errors.
#'
#' @returns The GMM loss.
#' @export gmmFun
#'
gmmFun <- function(theta, W, designMat, IMat, weightsMat, dataAgD) {
  m <- moment_residual(theta, designMat, IMat, weightsMat, dataAgD)
  as.numeric(t(m) %*% W %*% m)
}
