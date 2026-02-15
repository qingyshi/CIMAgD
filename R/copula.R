
#' Gaussian copula model for generating covariates
#'
#' @param summaryX A list describing one trial or multiple trials.
#' Supported formats are:
#' 1) single-trial list with `covariate`, `distribution`, `sumstats`, and optional `n`;
#' 2) multi-trial list in `summaryX$trials` (named list of single-trial specs);
#' 3) named list of single-trial specs.
#' @param N Number of samples to generate. Can be scalar or a named vector (for multiple trials).
#'
#' @returns A data.frame for one trial, or a named list of data.frames for multiple trials.
#' @export XCopula
#'
XCopula <- function(summaryX, N = 1000) {
  resolve_qfun <- function(dist_name) {
    qname <- paste0("q", tolower(dist_name))

    # Preferred in installed-package contexts (e.g., R CMD check).
    ns <- tryCatch(asNamespace("CIMAgD"), error = function(e) NULL)
    if (!is.null(ns) && exists(qname, mode = "function", envir = ns, inherits = FALSE)) {
      return(get(qname, mode = "function", envir = ns, inherits = FALSE))
    }

    # Fallback for sourced-file contexts.
    if (exists(qname, mode = "function", envir = environment(), inherits = TRUE)) {
      return(get(qname, mode = "function", envir = environment(), inherits = TRUE))
    }

    stop("Quantile function not found for distribution: ", dist_name)
  }

  is_single_spec <- function(x) {
    is.list(x) && ("sumstats" %in% names(x)) && ("covariate" %in% names(x))
  }

  parse_trials <- function(summaryX) {
    if (is_single_spec(summaryX)) {
      out <- list(single = summaryX)
      return(out)
    }
    if ("trials" %in% names(summaryX) && is.list(summaryX$trials)) {
      return(summaryX$trials)
    }
    if (all(vapply(summaryX, is_single_spec, logical(1)))) {
      return(summaryX)
    }
    stop("`summaryX` has unsupported structure.")
  }

  get_n <- function(spec, trial_name, N) {
    if (length(N) == 1L) {
      if (!is.null(spec$n) && is.finite(spec$n) && spec$n > 0) {
        return(as.integer(spec$n))
      }
      return(as.integer(N))
    }
    if (!is.null(names(N)) && trial_name %in% names(N)) {
      return(as.integer(N[[trial_name]]))
    }
    if (!is.null(spec$n) && is.finite(spec$n) && spec$n > 0) {
      return(as.integer(spec$n))
    }
    stop("Missing sample size for trial '", trial_name, "'.")
  }

  build_corr_param <- function(corr, p) {
    if (is.null(corr)) {
      return(rep(0, p * (p - 1) / 2))
    }
    if (is.matrix(corr)) {
      if (!all(dim(corr) == c(p, p))) {
        stop("Correlation matrix has wrong dimension.")
      }
      return(copula::P2p(corr))
    }
    corr <- as.numeric(corr)
    if (length(corr) == p * (p - 1) / 2) {
      return(corr)
    }
    stop("Correlation input has wrong length.")
  }

  second_to_param2 <- function(ss, mean, dist) {
    p <- length(mean)
    param2 <- rep(NA_real_, p)

    if (!is.null(ss$param2)) {
      param2 <- as.numeric(ss$param2)
    } else if (!is.null(ss$sd)) {
      param2 <- as.numeric(ss$sd)
    } else if (!is.null(ss$var)) {
      param2 <- sqrt(pmax(as.numeric(ss$var), 1e-8))
    } else if (!is.null(ss$second)) {
      second <- as.numeric(ss$second)
      param2 <- sqrt(pmax(second - mean^2, 1e-8))
    }

    for (k in seq_len(p)) {
      if (is.na(param2[k])) {
        if (dist[k] %in% c("binomial")) {
          param2[k] <- 1
        } else {
          param2[k] <- 1
        }
      }
    }
    param2
  }

  gen_one <- function(spec, trial_name, N) {
    covariate <- as.character(spec$covariate)
    dist <- as.character(spec$distribution)
    if (length(dist) == 1L) {
      dist <- rep(dist, length(covariate))
    }

    ss <- spec$sumstats
    ss_mean <- as.numeric(ss$mean)
    if (length(ss_mean) != length(covariate)) {
      stop("`sumstats$mean` length mismatch in trial '", trial_name, "'.")
    }
    ss_param2 <- second_to_param2(ss, ss_mean, dist)
    ss_min <- if (!is.null(ss$min)) ss$min else rep(NA_real_, length(covariate))
    ss_max <- if (!is.null(ss$max)) ss$max else rep(NA_real_, length(covariate))
    ss_corr <- build_corr_param(ss$corr, length(covariate))

    n_trial <- get_n(spec, trial_name, N)
    u_qmc <- randtoolbox::sobol(n = n_trial, dim = length(covariate))
    cop <- copula::normalCopula(ss_corr, dim = length(covariate), dispstr = "un")
    u_cop <- copula::cCopula(u_qmc, cop, inverse = TRUE)

    X <- sapply(seq_along(covariate), function(k) {
      qfun <- resolve_qfun(dist[k])
      qfun(u_cop[, k], ss_mean[k], ss_param2[k], ss_min[k], ss_max[k])
    })
    X <- as.data.frame(X)
    names(X) <- covariate
    X
  }

  trials <- parse_trials(summaryX)
  trial_names <- names(trials)
  if (is.null(trial_names)) {
    trial_names <- paste0("trial", seq_along(trials))
  }

  out <- lapply(seq_along(trials), function(i) {
    gen_one(trials[[i]], trial_names[i], N = N)
  })
  names(out) <- trial_names

  if (length(out) == 1L) {
    out[[1]]
  } else {
    out
  }

}
