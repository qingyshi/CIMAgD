# Normalize point estimates to marginal parameters and latent Gaussian dependence.
margin_specification <- function(x, variable) {
  x <- as.list(x)
  family <- tolower(x$distribution %||% x$dist %||% "")
  if (length(family) != 1L || is.na(family) || !family %in%
      c("normal", "bernoulli", "gamma", "beta", "lognormal", "exponential"))
    stop("Unsupported marginal distribution for ", variable,
         ". Use normal, bernoulli, gamma, beta, lognormal or exponential.", call. = FALSE)
  has <- function(nm) !is.null(x[[nm]]) && length(x[[nm]]) == 1L && !is.na(x[[nm]])
  get_number <- function(nm, positive = FALSE) {
    value <- x[[nm]]
    if (!has(nm) || !is.numeric(value) || !is.finite(value) || (positive && value <= 0))
      stop("Margin '", variable, "' needs a ", if (positive) "positive " else "finite ", nm, ".", call. = FALSE)
    as.numeric(value)
  }
  out <- as.list(stats::setNames(rep(NA_real_, 10),
    c("mean", "sd", "prob", "shape", "rate", "shape1", "shape2", "meanlog", "sdlog", "scale")))
  parameterization <- x$parameterization %||% ""
  if (family == "normal") {
    out$mean <- get_number("mean"); out$sd <- get_number("sd", TRUE)
    parameterization <- "moments"
  } else if (family == "bernoulli") {
    out$prob <- if (has("prob")) get_number("prob") else get_number("mean")
    if (out$prob <= 0 || out$prob >= 1) stop("Bernoulli probabilities must be strictly between zero and one.", call. = FALSE)
    out$mean <- out$prob; out$sd <- sqrt(out$prob * (1 - out$prob)); parameterization <- "probability"
  } else if (family == "gamma") {
    native <- if (nzchar(parameterization)) parameterization == "native" else has("shape")
    if (native) {
      out$shape <- get_number("shape", TRUE)
      out$rate <- if (has("rate")) get_number("rate", TRUE) else 1 / get_number("scale", TRUE)
      out$mean <- out$shape / out$rate; out$sd <- sqrt(out$shape) / out$rate
    } else {
      out$mean <- get_number("mean", TRUE); out$sd <- get_number("sd", TRUE)
      out$shape <- (out$mean / out$sd)^2; out$rate <- out$mean / out$sd^2
    }
    out$scale <- 1 / out$rate; parameterization <- if (native) "native" else "moments"
  } else if (family == "beta") {
    native <- if (nzchar(parameterization)) parameterization == "native" else has("shape1") || has("shape2")
    if (native) {
      out$shape1 <- get_number("shape1", TRUE); out$shape2 <- get_number("shape2", TRUE)
      total <- out$shape1 + out$shape2
      out$mean <- out$shape1 / total
      out$sd <- sqrt(out$shape1 * out$shape2 / (total^2 * (total + 1)))
    } else {
      out$mean <- get_number("mean"); out$sd <- get_number("sd", TRUE)
      if (out$mean <= 0 || out$mean >= 1 || out$sd^2 >= out$mean * (1 - out$mean))
        stop("Beta moments need 0 < mean < 1 and 0 < sd^2 < mean*(1-mean).", call. = FALSE)
      total <- out$mean * (1 - out$mean) / out$sd^2 - 1
      out$shape1 <- out$mean * total; out$shape2 <- (1 - out$mean) * total
    }
    parameterization <- if (native) "native" else "moments"
  } else if (family == "lognormal") {
    native <- if (nzchar(parameterization)) parameterization == "log-moments" else has("meanlog") || has("sdlog")
    if (native) {
      out$meanlog <- get_number("meanlog"); out$sdlog <- get_number("sdlog", TRUE)
      out$mean <- exp(out$meanlog + out$sdlog^2 / 2)
      out$sd <- out$mean * sqrt(expm1(out$sdlog^2))
    } else {
      out$mean <- get_number("mean", TRUE); out$sd <- get_number("sd", TRUE)
      out$sdlog <- sqrt(log1p((out$sd / out$mean)^2))
      out$meanlog <- log(out$mean) - out$sdlog^2 / 2
    }
    parameterization <- if (native) "log-moments" else "moments"
  } else {
    out$rate <- if (has("rate")) get_number("rate", TRUE) else 1 / get_number("mean", TRUE)
    out$mean <- out$sd <- 1 / out$rate; parameterization <- "rate"
  }
  if (any(!is.finite(unlist(out)[!is.na(unlist(out))]))) stop("Marginal parameters overflow for ", variable, ".", call. = FALSE)
  as.data.frame(c(list(variable = variable, distribution = family, parameterization = parameterization), out),
                stringsAsFactors = FALSE)
}

validate_summary <- function(summaryX) {
  if (!is.list(summaryX) || !is.null(summaryX$trials))
    stop("summaryX must describe one source, not multiple trials.", call. = FALSE)
  if (is.null(summaryX$marginals) && !is.null(summaryX$covariate) && !is.null(summaryX$sumstats)) {
    ss <- summaryX$sumstats; variables <- summaryX$covariate; p <- length(variables)
    distribution <- rep_len(summaryX$distribution, p)
    if (any(distribution == "binomial") && !is.null(ss$param2) && any(ss$param2[distribution == "binomial"] != 1))
      stop("Only Bernoulli (binomial size 1) discrete margins are supported.", call. = FALSE)
    distribution[distribution == "binomial"] <- "bernoulli"
    sd <- ss$sd %||% ss$param2
    if (is.null(sd) && !is.null(ss$var)) sd <- sqrt(ss$var)
    if (is.null(sd) && !is.null(ss$second)) sd <- sqrt(ss$second - ss$mean^2)
    sd <- sd %||% rep(NA_real_, p)
    correlation <- ss$corr %||% diag(p)
    if (!is.matrix(correlation)) {
      if (length(correlation) != p * (p - 1L) / 2L) stop("Legacy latent correlation length mismatch.", call. = FALSE)
      r <- diag(p); r[lower.tri(r)] <- correlation; r[upper.tri(r)] <- t(r)[upper.tri(r)]; correlation <- r
    }
    summaryX <- list(n = summaryX$n, marginals = data.frame(variable = variables,
      distribution = distribution, mean = ss$mean, sd = sd), correlation = correlation)
  }
  check_scalar(summaryX$n, "summaryX$n", 2, TRUE)
  margins <- summaryX$marginals
  if (is.data.frame(margins)) {
    if (!nrow(margins) || !all(c("variable", "distribution") %in% names(margins)))
      stop("The marginals table needs variable and distribution columns.", call. = FALSE)
    variables <- as.character(margins$variable)
    specifications <- lapply(seq_len(nrow(margins)), function(k) as.list(margins[k, , drop = FALSE]))
  } else if (is.list(margins) && length(margins) && !is.null(names(margins))) {
    variables <- names(margins); specifications <- margins
  } else stop("summaryX$marginals must be a named list of marginal specifications or a data frame.", call. = FALSE)
  if (anyNA(variables) || any(!nzchar(variables)) || anyDuplicated(variables))
    stop("Marginal variable names must be unique and nonempty.", call. = FALSE)
  m <- do.call(rbind, Map(margin_specification, specifications, variables)); rownames(m) <- NULL
  p <- nrow(m); R <- summaryX$correlation %||% diag(p)
  if (!is.null(rownames(R)) || !is.null(colnames(R))) {
    if (!setequal(rownames(R), variables) || !setequal(colnames(R), variables))
      stop("Correlation names must match marginal variable names.", call. = FALSE)
    R <- R[variables, variables, drop = FALSE]
  }
  R <- validate_covariance(R, p, "summaryX$correlation", positive = TRUE)
  if (max(abs(diag(R) - 1)) > 1e-10) stop("summaryX$correlation must have unit diagonal.", call. = FALSE)
  if (!is.null(summaryX$correlation_scale) && summaryX$correlation_scale != "latent")
    stop("Supply latent Gaussian correlations, not observed Pearson correlations for mixed margins.", call. = FALSE)
  dimnames(R) <- list(variables, variables)
  list(n = summaryX$n, marginals = m, correlation = R, correlation_scale = "latent", label = summaryX$label)
}

#' Construct a single-source Gaussian copula specification
#'
#' @param n Original source sample size.
#' @param marginals Named list of marginal specifications, or a data frame with
#'   variable and distribution columns. Normal: mean, sd. Bernoulli: prob (or
#'   mean). Gamma: mean, sd or shape and rate (or scale). Beta on (0,1): mean, sd
#'   or shape1, shape2. Lognormal: arithmetic mean, sd or meanlog, sdlog.
#'   Exponential: mean or rate. These are point estimates, not their SEs.
#' @param correlation Latent Gaussian correlation matrix, optionally named;
#'   NULL means independence. Observed correlations are not latent correlations.
#' @param label Optional source description.
#' @return A validated summaryX list, with a normalized marginal-parameter table.
#'   No sampling covariance of the supplied point estimates is required.
#' @examples
#' s <- copula_summary(1000, list(
#'   age = list(distribution = "normal", mean = 60, sd = 10),
#'   marker = list(distribution = "gamma", mean = 2, sd = 1),
#'   score = list(distribution = "beta", mean = .4, sd = .1),
#'   diabetes = list(distribution = "bernoulli", prob = .25)))
#' head(XCopula(s, N = 10, seed = 1))
#' @export
copula_summary <- function(n, marginals, correlation = NULL, label = NULL) {
  validate_summary(list(n = n, marginals = marginals, correlation = correlation, label = label))
}

normal_quantile_transform <- function(z, qfun, args) {
  result <- numeric(length(z)); lower <- z <= 0
  for (tail in c(TRUE, FALSE)) {
    idx <- which(lower == tail)
    result[idx] <- do.call(qfun, c(list(p = stats::pnorm(z[idx], lower.tail = tail, log.p = TRUE),
      lower.tail = tail, log.p = TRUE), args))
  }
  result
}

copula_from_noise <- function(spec, noise) {
  z <- noise %*% chol(spec$correlation); m <- spec$marginals
  out <- lapply(seq_len(nrow(m)), function(k) {
    x <- m[k, ]; latent <- z[, k]
    switch(x$distribution,
      normal = x$mean + x$sd * latent,
      bernoulli = as.numeric(latent > stats::qnorm(1 - x$prob)),
      gamma = normal_quantile_transform(latent, stats::qgamma, list(shape = x$shape, rate = x$rate)),
      beta = normal_quantile_transform(latent, stats::qbeta, list(shape1 = x$shape1, shape2 = x$shape2)),
      lognormal = exp(x$meanlog + x$sdlog * latent),
      exponential = -stats::pnorm(latent, lower.tail = FALSE, log.p = TRUE) / x$rate)
  })
  names(out) <- m$variable
  if (any(!is.finite(unlist(out)))) stop("Nonfinite generated covariates; check marginal parameters.", call. = FALSE)
  as.data.frame(out)
}

#' Draw covariates from a single-source Gaussian copula
#'
#' @param summaryX Point-estimate specification, as constructed by
#'   [copula_summary()]. Named marginal lists and marginal tables are accepted.
#' @param N Number of iid draws; defaults to summaryX$n. CIMA's plugin Q always
#'   uses the original n, irrespective of this standalone simulation option.
#' @param seed Optional reproducible seed, restoring the caller's RNG state.
#' @return A covariate data frame. The single-source covariate/distribution/
#'   sumstats format remains accepted; multi-source lists are rejected.
#' @export
XCopula <- function(summaryX, N = NULL, seed = NULL) {
  spec <- validate_summary(summaryX)
  N <- N %||% spec$n; check_scalar(N, "N", 2, TRUE)
  with_seed(seed, copula_from_noise(spec, matrix(stats::rnorm(N * nrow(spec$marginals)), N)))
}

latent_scores <- function(x, margin) {
  family <- margin$distribution
  if (family == "normal") return((x - margin$mean) / margin$sd)
  if (family == "lognormal") return((log(x) - margin$meanlog) / margin$sdlog)
  args <- switch(family, gamma = list(shape = margin$shape, rate = margin$rate),
    beta = list(shape1 = margin$shape1, shape2 = margin$shape2), exponential = list(rate = margin$rate))
  fun <- switch(family, gamma = stats::pgamma, beta = stats::pbeta, exponential = stats::pexp)
  lower <- do.call(fun, c(list(q = x, log.p = TRUE), args))
  upper <- do.call(fun, c(list(q = x, lower.tail = FALSE, log.p = TRUE), args))
  ifelse(lower < log(.5), stats::qnorm(lower, log.p = TRUE),
         stats::qnorm(upper, lower.tail = FALSE, log.p = TRUE))
}

refit_copula <- function(data, spec) {
  out <- spec; m <- spec$marginals; p <- nrow(m)
  out$marginals <- do.call(rbind, lapply(seq_len(p), function(k) {
    x <- data[[m$variable[k]]]; family <- m$distribution[k]
    updated <- list(distribution = family, mean = mean(x), sd = stats::sd(x))
    if (family == "lognormal" && m$parameterization[k] == "log-moments")
      updated <- list(distribution = family, meanlog = mean(log(x)), sdlog = stats::sd(log(x)))
    margin_specification(updated, m$variable[k])
  }))
  scores <- lapply(seq_len(p), function(k) if (m$distribution[k] == "bernoulli") NULL else
    latent_scores(data[[m$variable[k]]], out$marginals[k, ]))
  if (any(!is.finite(unlist(scores)))) stop("Nonfinite latent scores in copula refit.", call. = FALSE)
  R <- diag(p)
  if (p > 1L) for (j in seq_len(p - 1L)) for (k in seq.int(j + 1L, p)) {
    x <- data[[m$variable[j]]]; y <- data[[m$variable[k]]]
    continuous_j <- m$distribution[j] != "bernoulli"; continuous_k <- m$distribution[k] != "bernoulli"
    if (continuous_j && continuous_k) rho <- stats::cor(scores[[j]], scores[[k]]) else {
      if (continuous_j || continuous_k) {
        z <- if (continuous_j) scores[[j]] else scores[[k]]
        binary <- if (continuous_j) y else x
        threshold <- stats::qnorm(1 - mean(binary))
        objective <- function(rho) {
          argument <- (rho * z - threshold) / sqrt(1 - rho^2)
          -sum(ifelse(binary == 1, stats::pnorm(argument, log.p = TRUE),
                      stats::pnorm(argument, lower.tail = FALSE, log.p = TRUE)))
        }
      } else {
        a <- stats::qnorm(1 - mean(x)); b <- stats::qnorm(1 - mean(y))
        counts <- c(sum(x == 0 & y == 0), sum(x == 0 & y == 1), sum(x == 1 & y == 0), sum(x == 1 & y == 1))
        objective <- function(rho) {
          p00 <- stats::integrate(function(z) stats::dnorm(z) * stats::pnorm((b - rho * z) / sqrt(1 - rho^2)),
            -Inf, a, rel.tol = 1e-7)$value
          probs <- c(p00, stats::pnorm(a) - p00, stats::pnorm(b) - p00, 1 - stats::pnorm(a) - stats::pnorm(b) + p00)
          -sum(counts * log(pmax(probs, 1e-15)))
        }
      }
      rho <- stats::optimize(objective, c(-.995, .995), tol = 1e-6)$minimum
    }
    R[j, k] <- R[k, j] <- rho
  }
  original <- R
  if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8)
    R <- as.matrix(Matrix::nearPD(R, corr = TRUE, posd.tol = 1e-7)$mat)
  dimnames(R) <- dimnames(spec$correlation); out$correlation <- R
  out$correlation_projection <- max(abs(R - original))
  validate_summary(out)
  out
}
