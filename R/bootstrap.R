# Moment-based first-stage bootstrap weights avoid nested bootstrap fitting.
prepare_copula_bootstrap <- function(spec, noise, B, d, mu, formula_cate, formula_tilting,
                                     trial_info, match_proportions, time_indexed, times,
                                     W, control, theta_reference, method = "common") {
  summaries <- vector("list", B); moments <- vector("list", B); failures <- character(B)
  control$jacobian_diagnostics <- FALSE
  for (b in seq_len(B)) {
    value <- tryCatch({
      source <- copula_from_noise(spec, matrix(stats::rnorm(spec$n * nrow(spec$marginals)), spec$n))
      refitted <- refit_copula(source, spec)
      core <- fit_cima_core(d, copula_from_noise(refitted, noise), mu, formula_cate,
        formula_tilting, trial_info, match_proportions, time_indexed, times, W,
        control, diagnostics = FALSE, method = method)
      list(summary = refitted, moment = as.vector(core$fit$design %*% theta_reference))
    }, error = function(e) e)
    if (inherits(value, "error")) failures[b] <- conditionMessage(value) else {
      summaries[[b]] <- value$summary; moments[[b]] <- value$moment
    }
  }
  ok <- !nzchar(failures)
  check_bootstrap_success(ok, failures, B, warn = FALSE)
  list(summaries = summaries, failures = failures,
       moment_covariance = stats::cov(do.call(rbind, moments[ok])))
}

check_bootstrap_success <- function(ok, failures, B, warn = TRUE) {
  if (sum(ok) < max(2L, ceiling(.8 * B))) stop("Too many failed copula bootstrap fits (", sum(!ok), "/", B,
    "). First failure: ", failures[which(!ok)[1L]], call. = FALSE)
  if (warn && any(!ok)) warning(sum(!ok), " copula bootstrap fits failed; see copula$bootstrap$failures.", call. = FALSE)
}

bootstrap_copula <- function(spec, noise, B, d, mu, formula_cate, formula_tilting,
                             trial_info, match_proportions, time_indexed, times,
                             W, control, target, shared, target_times, subgroups,
                             prepared = NULL, weight_psd = FALSE, theta_reference, method = "common") {
  draws <- summaries <- moments <- vector("list", B); failures <- character(B)
  control$jacobian_diagnostics <- FALSE
  for (b in seq_len(B)) {
    value <- tryCatch({
      if (is.null(prepared)) {
        source <- copula_from_noise(spec, matrix(stats::rnorm(spec$n * nrow(spec$marginals)), spec$n))
        refitted <- refit_copula(source, spec)
      } else {
        if (nzchar(prepared$failures[b])) stop(prepared$failures[b], call. = FALSE)
        refitted <- prepared$summaries[[b]]
      }
      Q_b <- copula_from_noise(refitted, noise)
      core <- fit_cima_core(d, Q_b, mu, formula_cate, formula_tilting, trial_info,
        match_proportions, time_indexed, times, W, control, diagnostics = FALSE,
        weight_psd = weight_psd, method = method)
      target_b <- if (shared) Q_b else target
      tf <- fit_target(target_b, core$context$basis, core$fit, target_times, subgroups, NULL, shared)
      list(theta = core$fit$theta, target = stats::setNames(tf$estimates$estimate, rownames(tf$J_psi_theta)),
           summary = refitted, moment = as.vector(core$fit$design %*% theta_reference))
    }, error = function(e) e)
    if (inherits(value, "error")) failures[b] <- conditionMessage(value) else {
      draws[[b]] <- c(value$theta, value$target); summaries[[b]] <- value$summary; moments[[b]] <- value$moment
    }
  }
  ok <- !nzchar(failures)
  check_bootstrap_success(ok, failures, B)
  samples <- do.call(rbind, draws[ok]); covariance <- stats::cov(samples)
  p <- length(theta_reference)
  list(draws = samples, covariance = covariance,
       theta_covariance = covariance[seq_len(p), seq_len(p), drop = FALSE],
       target_covariance = covariance[-seq_len(p), -seq_len(p), drop = FALSE],
       moment_covariance = stats::cov(do.call(rbind, moments[ok])),
       refitted_summaries = summaries[ok], requested = B, successful = sum(ok),
       failures = failures[!ok], integration_n = nrow(noise),
       weight_refitted = FALSE,
       note = paste("Conditional on trial reports; bootstrap Q variation replaces analytic Q variation.",
         "Independent target and trial components remain analytic.",
         if (weight_psd) "The estimated second-stage weight is held fixed across replicates, using its first-order equivalence under valid GMM moments." else ""))
}
