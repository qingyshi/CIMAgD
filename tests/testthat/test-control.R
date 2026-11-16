test_that("stable trial-specific working control fits are retained", {
  fit <- do.call(CIMA, two_trials())
  wc <- fit$variance$working_control
  expect_equal(wc$trial, c("s1", "s2"))
  expect_equal(wc$model, rep("trial-specific", 2))
  expect_equal(wc$trial_specific_issues, rep("none", 2))
  expect_true(all(is.na(wc$common_issues)))
  expect_equal(fit$variance$trials$s1$control$model, "trial-specific")
  expect_match(fit$variance$trials$s1$method, "trial-specific working control mean")
})

test_that("an unidentified trial-specific working mean switches to the common fit only for that trial", {
  a <- two_trials()
  a$data_agd <- a$data_agd[c(1:5, 6:8), ]  # s2 reports only x subgroups
  fit <- do.call(CIMA, a)
  wc <- fit$variance$working_control
  expect_equal(wc$model, c("trial-specific", "common"))
  expect_match(wc$trial_specific_issues[2], "unidentified working terms")
  expect_equal(wc$common_issues[2], "none")
  s1 <- fit$variance$trials$s1$control; s2 <- fit$variance$trials$s2$control
  alone <- do.call(CIMA, within(a, { data_agd <- data_agd[1:5, ]; mu <- mu["s1"] }))
  expect_equal(s1$beta, alone$variance$trials$s1$control$beta, tolerance = 1e-12)
  # The common working mean is the inverse-variance least-squares fit to every
  # trial's control reports, each aggregated under its own tilted law.
  d <- fit$data
  oracle <- stats::lm.wfit(fit$CATE$J_m_theta, d$nu0, 1 / d$se_nu0^2)$coefficients
  expect_false(s2$constrained)
  expect_equal(unname(s2$beta), unname(oracle), tolerance = 1e-10)
  expect_equal(s2$fitted_rows, d$row_id)
  expect_equal(s2$gamma1 - s2$gamma0, as.vector(fit$basis$basis_q[[1]] %*% coef(fit)), tolerance = 1e-10)
  expect_output(print(fit), "Common working control mean \\(A.6.3\\) used for: s2")
  expect_gte(min(eigen(fit$variance$trials$s2$Sigma, symmetric = TRUE)$values), -1e-10)
})

test_that("MVLP borrows the common working control mean but keeps trial-specific CATEs", {
  a <- two_trials()
  a$data_agd <- a$data_agd[c(1:5, 6:8), ]
  a$formula_cate <- ~ x
  a$formula_control <- ~ x + z
  fit <- do.call(CIMA, c(a, list(method = "MVLP")))
  expect_equal(fit$variance$working_control$model, c("trial-specific", "common"))
  s2 <- fit$variance$trials$s2$control
  g2 <- as.vector(cbind(1, a$Q$x) %*% coef(fit)[, "s2"])
  expect_equal(s2$gamma1 - s2$gamma0, g2, tolerance = 1e-10)
  A <- crossprod(fit$moments$alpha, model.matrix(~ x + z, a$Q)) / nrow(a$Q)
  oracle <- stats::lm.wfit(A, fit$data$nu0, 1 / fit$data$se_nu0^2)$coefficients
  expect_equal(unname(s2$beta), unname(oracle), tolerance = 1e-10)
  # Borrowing changes only the working covariance, not the trial-specific CATEs.
  b <- a; b$data_agd <- a$data_agd[6:8, ]; b$mu <- a$mu["s2"]
  alone <- do.call(CIMA, b)
  expect_equal(unname(coef(fit)[, "s2"]), unname(coef(alone)), tolerance = 1e-10)
  expect_equal(alone$variance$working_control$model, "trial-specific")
  expect_match(alone$variance$working_control$trial_specific_issues, "unidentified")
})

test_that("active risk bounds trigger the common fit when it is interior", {
  args <- small_example(); d <- prepare_agd(args$data_agd)
  Q <- args$Q; h <- model.matrix(~ x + z, Q)
  ii <- indicator_matrix(Q, d); alpha <- sweep(ii, 2, colMeans(ii), "/")
  d$nu0 <- c(.5, .05, .95, .5, .5); d$se_nu0 <- rep(.02, 5)
  g <- rep(.2, nrow(Q)); ctl <- list(rank_tol = 1e-9, condition_limit = 1e8)
  A <- crossprod(alpha, h) / nrow(Q)
  pooled <- function() list(A = rbind(A, A), nu0 = c(d$nu0, .4, .35, .45, .38, .42),
    se0 = c(d$se_nu0, rep(.002, 5)), row_id = c(d$row_id, paste0("other:", 1:5)))
  chosen <- select_working_control("s1", h, alpha, rep(1, nrow(Q)), g, d, "binary", ctl, pooled)
  expect_equal(chosen$model, "common")
  expect_equal(chosen$selection$trial_specific_issues, "active risk bounds")
  expect_length(chosen$selection$common_issues, 0)
  expect_false(chosen$constrained)
  expect_equal(chosen$gamma1 - chosen$gamma0, g, tolerance = 1e-10)
  # A common fit that fails leaves the feasible trial-specific fit in place.
  broken <- function() list(A = A[, 1:2], nu0 = d$nu0, se0 = d$se_nu0, row_id = d$row_id)
  kept <- select_working_control("s1", h, alpha, rep(1, nrow(Q)), g, d, "binary", ctl, broken)
  expect_equal(kept$model, "trial-specific")
  expect_true(kept$constrained)
  expect_equal(kept$selection$common_issues, "estimation failure")
  expect_match(kept$selection$reason, "trial-specific fit retained")
})

test_that("infeasible working risks fail only when no working model is available", {
  args <- small_example(); d <- prepare_agd(args$data_agd)
  Q <- args$Q; b <- matrix(1, nrow(Q), 1, dimnames = list(NULL, "(Intercept)"))
  ii <- indicator_matrix(Q, d); alpha <- sweep(ii, 2, colMeans(ii), "/")
  g <- ifelse(Q$x == 1, .6, -.6); ctl <- list(rank_tol = 1e-9)
  pooled <- function() list(A = crossprod(alpha, b) / nrow(Q), nu0 = d$nu0, se0 = d$se_nu0, row_id = d$row_id)
  expect_error(select_working_control("s1", b, alpha, rep(1, nrow(Q)), g, d, "binary", ctl, pooled),
               "No working control model for trial s1")
  expect_error(select_working_control("s1", b, alpha, rep(1, nrow(Q)), rep(1.2, nrow(Q)), d, "binary", ctl, pooled),
               "leaves \\[-1, 1\\]")
  expect_true(lexicographic_le(c(0, 0, 0, 1, 0), c(0, 1, 0, 0, 0)))
  expect_false(lexicographic_le(c(0, 1, 0, 0, 0), c(0, 0, 1, 1, 1)))
  expect_true(lexicographic_le(c(0, 0, 0, 1, 0), c(0, 0, 0, 1, 0)))
})
