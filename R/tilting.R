# internal helper: numeric-stable
safe_exp <- function(x) {
  exp(pmax(pmin(x, 700), -700))
}


# internal helper: estimating equation for tilting
tiltingloss <- function(lambda, tiltingMat, mu) {
  h_plus <- as.matrix(tiltingMat)
  if (length(mu) == ncol(h_plus) - 1L) {
    mu_plus <- c(1, mu)
  } else if (length(mu) == ncol(h_plus)) {
    mu_plus <- mu
  } else {
    stop("`mu` has incompatible length for `tiltingMat`.")
  }

  lp <- as.numeric(h_plus %*% lambda)
  w <- safe_exp(lp)

  colMeans(h_plus * w) - mu_plus
}


#' Exponential tilting model
#'
#' @param mu A list of tilting targets for each trial.
#' @param tiltingMat The tilting design matrix by auxiliary data.
#' @param dataAgD A data.frame with the treatment effects and standard errors.
#'
#' @returns A list of tilting parameters for each trial.
#' @export tiltingModel
#'
tiltingModel <- function(mu, tiltingMat, dataAgD) {
  trials <- unique(as.character(dataAgD$trial))
  pars <- vector("list", length(trials))
  names(pars) <- trials

  for (s in trials) {
    mu_s <- mu[[s]]
    if (is.null(mu_s)) {
      stop(sprintf("Missing `mu` for trial '%s'.", s))
    }
    fit <- try(
      rootSolve::multiroot(
        f = tiltingloss,
        start = rep(0, ncol(tiltingMat)),
        tiltingMat = tiltingMat,
        mu = mu_s,
        rtol = 1e-10,
        atol = 1e-10,
        ctol = 1e-10,
        maxiter = 200
      ),
      silent = TRUE
    )
    if (inherits(fit, "try-error")) {
      stop(sprintf("Tilting root solve failed for trial '%s'.", s))
    }
    pars[[s]] <- fit$root
  }

  pars
}


#' Tilting weights matrix
#'
#' @param tiltingMat The tilting matrix.
#' @param tiltingPars A list of tilting parameters for each trial.
#' @param dataAgD A data.frame with the treatment effects and standard errors.
#'
#' @returns A list of matrices of tilting weights.
#'
weightsMatFun <- function(tiltingMat, tiltingPars, dataAgD) {
  trials <- unique(as.character(dataAgD$trial))
  tiltingMat <- as.matrix(tiltingMat)

  out <- lapply(trials, function(s) {
    Js <- sum(dataAgD$trial == s)
    lp <- as.numeric(tiltingMat %*% tiltingPars[[s]])
    w <- safe_exp(lp)
    matrix(w, nrow = nrow(tiltingMat), ncol = Js)
  })
  names(out) <- trials
  out
}
