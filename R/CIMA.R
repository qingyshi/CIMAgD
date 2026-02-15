
#' Causally-interpretable meta-analysis using aggregate data
#'
#' @param dataAgD A data.frame with treatment effects and standard errors.
#' The data.frame must contain columns: `trial`, `subgroup`, `stratum`, `TE`,
#' `seTE`, `n`, and `prop`. Columns `trial`, `subgroup`, `stratum`, `TE`, and
#' `seTE` must not contain missing values. Columns `n` and `prop` may contain
#' missing values.
#' @param mu A list of tilting targets for each trial.
#' @param Q A data.frame of covariates that represent a auxiliary distribution.
#' @param target A data.frame of covariates for the target population.
#' @param formulaCATE A character string of the formula for CATE.
#' @param formulaTilting A character string of the formula for tilting.
#' @param V A list of variance-covariance matrices of treatment effects. Maybe missing.
#'
#' @returns A object of CIMA.
#' @export CIMA
#'
CIMA <- function(dataAgD, mu, Q, target, formulaCATE, formulaTilting, V = NULL) {
  req_cols <- c("trial", "subgroup", "stratum", "TE", "seTE", "n", "prop")
  miss_cols <- setdiff(req_cols, names(dataAgD))
  if (length(miss_cols) > 0) {
    stop("Missing columns in `dataAgD`: ", paste(miss_cols, collapse = ", "))
  }
  key_cols <- c("trial", "subgroup", "stratum", "TE", "seTE")
  has_na <- vapply(key_cols, function(x) any(is.na(dataAgD[[x]])), logical(1))
  if (any(has_na)) {
    stop("Missing values are not allowed in: ", paste(key_cols[has_na], collapse = ", "))
  }
  for (s in unique(as.character(dataAgD$trial))) {
    idx <- as.character(dataAgD$trial) == s
    n0 <- suppressWarnings(max(dataAgD$n[idx & dataAgD$subgroup == "overall"], na.rm = TRUE))
    if (is.finite(n0) && n0 > 0) {
      fill_idx <- idx & (is.na(dataAgD$prop) | dataAgD$prop <= 0) & !is.na(dataAgD$n) & dataAgD$n > 0
      dataAgD$prop[fill_idx] <- dataAgD$n[fill_idx] / n0
    }
  }
  dataAgD$prop[is.na(dataAgD$prop) | dataAgD$prop <= 0] <- 1

  f_cate <- if (inherits(formulaCATE, "formula")) formulaCATE else stats::as.formula(formulaCATE)
  f_tilt <- if (inherits(formulaTilting, "formula")) formulaTilting else stats::as.formula(formulaTilting)

  Q <- as.data.frame(Q)
  target <- as.data.frame(target)

  designMat <- stats::model.matrix(f_cate, data = Q)
  designMat_target <- stats::model.matrix(f_cate, data = target)
  tiltingMat <- stats::model.matrix(f_tilt, data = Q)

  IMat <- IMatFun(Q = Q, dataAgD = dataAgD)
  tiltingPars <- tiltingModel(mu = mu, dataAgD = dataAgD, tiltingMat = tiltingMat)
  weightsMat <- weightsMatFun(tiltingMat = tiltingMat, dataAgD = dataAgD, tiltingPars = tiltingPars)

  Dtheta <- DthetaMat(designMat = designMat, IMat = IMat, weightsMat = weightsMat, dataAgD = dataAgD)
  W <- WMat(dataAgD = dataAgD, V = V)

  tau_hat <- as.numeric(dataAgD$TE)

  # GMM estimation
  M1 <- crossprod(Dtheta, W %*% Dtheta)
  M2 <- crossprod(Dtheta, W %*% tau_hat)
  theta_hat <- try(solve(M1, M2), silent = TRUE)
  if (inherits(theta_hat, "try-error")) {
    theta_hat <- stats::optim(
      par = rep(0, ncol(designMat)),
      fn = gmmFun,
      W = W,
      designMat = designMat,
      IMat = IMat,
      weightsMat = weightsMat,
      dataAgD = dataAgD,
      method = "BFGS"
    )$par
  } else {
    theta_hat <- as.numeric(theta_hat)
  }
  names(theta_hat) <- colnames(designMat)

  psi_hat <- mean(as.numeric(designMat_target %*% theta_hat))

  vc <- vcovFun(
    theta_hat = theta_hat, W = W, Dtheta = Dtheta,
    Q = Q, designMat = designMat, tiltingMat = tiltingMat, IMat = IMat,
    weightsMat = weightsMat, dataAgD = dataAgD, mu = mu, V = V
  )

  vcov_theta <- vc$vcov_theta
  Jpsi <- colMeans(designMat_target)
  g_target <- as.numeric(designMat_target %*% theta_hat)
  phi_target <- g_target - psi_hat
  var_psi_target <- sum(phi_target^2) / (nrow(target)^2)
  var_psi <- as.numeric(var_psi_target + Jpsi %*% vcov_theta %*% Jpsi)

  results <- list(
    CATE = list(
      theta = theta_hat,
      vcov = vcov_theta,
      se = sqrt(diag(vcov_theta)),
      ci95 = cbind(
        lower = theta_hat - 1.96 * sqrt(diag(vcov_theta)),
        upper = theta_hat + 1.96 * sqrt(diag(vcov_theta))
      ),
      Dtheta = Dtheta
    ),
    ATE = list(
      ATE = psi_hat,
      var = var_psi,
      se = sqrt(var_psi),
      ci95 = c(psi_hat - 1.96 * sqrt(var_psi), psi_hat + 1.96 * sqrt(var_psi))
    ),
    moments = list(
      tau_hat = tau_hat,
      m_hat = moment_vector(theta_hat, designMat, IMat, weightsMat, dataAgD),
      residual = moment_residual(theta_hat, designMat, IMat, weightsMat, dataAgD),
      foc = as.numeric(crossprod(Dtheta, W %*% moment_residual(theta_hat, designMat, IMat, weightsMat, dataAgD)))
    ),
    tilting = list(
      formula = f_tilt,
      mu = mu,
      pars = tiltingPars
    ),
    call = match.call()
  )

  class(results) <- "CIMA"
  results
}


#' Print method for CIMA fit
#'
#' @param x A fitted object returned by [CIMA()].
#' @param digits Number of digits to print.
#' @param ... Unused.
#'
#' @returns Invisibly returns `x`.
#' @export
#'
print.CIMA <- function(x, digits = 4, ...) {
  fmt <- function(z) formatC(z, format = "f", digits = digits)

  cat("CIMAgD Fit\n")
  cat("==========\n")

  if (!is.null(x$ATE) && !is.null(x$ATE$ATE)) {
    ate <- as.numeric(x$ATE$ATE)
    ci <- x$ATE$ci95
    if (is.null(ci) && !is.null(x$ATE$se)) {
      ci <- c(ate - 1.96 * x$ATE$se, ate + 1.96 * x$ATE$se)
    }
    cat("\nATE\n")
    cat("---\n")
    if (!is.null(ci) && length(ci) == 2L) {
      cat("Estimate:", fmt(ate), "   95% CI:", paste0("[", fmt(ci[1]), ", ", fmt(ci[2]), "]"), "\n")
    } else {
      cat("Estimate:", fmt(ate), "\n")
    }
  }

  if (!is.null(x$CATE) && !is.null(x$CATE$theta)) {
    theta <- x$CATE$theta
    ci_mat <- x$CATE$ci95
    if (is.null(ci_mat) && !is.null(x$CATE$se)) {
      ci_mat <- cbind(
        lower = theta - 1.96 * x$CATE$se,
        upper = theta + 1.96 * x$CATE$se
      )
    }

    cat("\nCATE\n")
    cat("----\n")
    if (!is.null(ci_mat) && nrow(ci_mat) == length(theta)) {
      out <- data.frame(
        Estimate = as.numeric(theta),
        CI_Lower = as.numeric(ci_mat[, 1]),
        CI_Upper = as.numeric(ci_mat[, 2]),
        row.names = names(theta),
        check.names = FALSE
      )
      print(round(out, digits))
    } else {
      out <- data.frame(
        Estimate = as.numeric(theta),
        row.names = names(theta),
        check.names = FALSE
      )
      print(round(out, digits))
    }
  }

  invisible(x)
}
