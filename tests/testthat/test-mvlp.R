subgroup_x1 <- list(x1 = function(d) d$x == 1)

test_that("MVLP reproduces separate trial-specific fits and their variances", {
  a <- two_trials()
  fit <- do.call(CIMA, c(a, list(method = "MVLP", target_subgroups = subgroup_x1)))
  expect_equal(fit$method, "MVLP")
  expect_equal(dim(coef(fit)), c(3L, 2L))
  expect_equal(colnames(coef(fit)), c("s1", "s2"))
  expect_equal(names(fit$CATE$theta), paste0(rep(c("s1", "s2"), each = 3), ":", rownames(coef(fit))))
  for (s in c("s1", "s2")) {
    b <- a; b$data_agd <- a$data_agd[a$data_agd$trial == s, ]; b$mu <- a$mu[s]
    single <- do.call(CIMA, c(b, list(target_subgroups = subgroup_x1)))
    expect_equal(unname(coef(fit)[, s]), unname(coef(single)), tolerance = 1e-10)
    rows <- fit$target$trial_estimates$trial == s
    expect_equal(fit$target$trial_estimates$estimate[rows], single$target$estimates$estimate, tolerance = 1e-10)
    idx <- grep(paste0("^", s, ":"), rownames(vcov(fit, "trial")))
    expect_equal(unname(vcov(fit, "trial")[idx, idx]), unname(vcov(single, "target")), tolerance = 1e-10)
    par <- grep(paste0("^", s, ":"), rownames(vcov(fit)))
    expect_equal(unname(vcov(fit)[par, par]), unname(vcov(single)), tolerance = 1e-10)
  }
  # Trial-specific CATEs share no trial data, so their trial influence terms are separate.
  expect_equal(unname(fit$variance$trial_covariance[1:3, 4:6]), matrix(0, 3, 3))
  expect_gt(max(abs(fit$variance$q_covariance[1:3, 4:6])), 0)
})

test_that("MVLP covariance follows A.7 and pooling uses minimum-variance weights", {
  fit <- do.call(CIMA, c(two_trials(), list(method = "MVLP", target_subgroups = subgroup_x1)))
  tt <- fit$target; vv <- fit$variance
  # Shared target/Q: V_sj = Cov(u_s + v_s, u_j + v_j) / n_0 + 1{s = j} Var(r_s) / n_s.
  zeta <- tt$trial_phi_target + vv$theta_if_q %*% t(tt$trial_J_psi_theta)
  V <- crossprod(center_columns(zeta)) / nrow(zeta)^2 +
    tt$trial_J_psi_theta %*% vv$trial_covariance %*% t(tt$trial_J_psi_theta)
  expect_equal(vcov(fit, "trial"), V, ignore_attr = TRUE, tolerance = 1e-12)
  key <- paste(tt$trial_estimates$subgroup, tt$trial_estimates$Time, sep = "@")
  for (k in unique(key)) {
    idx <- which(key == k)
    x <- unname(solve(V[idx, idx], rep(1, 2)))
    expect_equal(unname(fit$pooling$weights[, k]), x / sum(x), tolerance = 1e-10)
    expect_equal(sum(fit$pooling$weights[, k]), 1)
    expect_equal(tt$estimates$estimate[match(k, rownames(tt$vcov))],
                 sum(x / sum(x) * tt$trial_estimates$estimate[idx]), tolerance = 1e-12)
    expect_equal(tt$vcov[k, k], 1 / sum(x), tolerance = 1e-12)
    expect_equal(tt$trial_estimates$weight[idx], x / sum(x), tolerance = 1e-10)
  }
  L <- fit$pooling$map
  expect_equal(vcov(fit, "target"), crossprod(L, V %*% L), ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(nrow(tt$estimates), 2L)
  expect_equal(fit$ATE$ATE, tt$estimates$estimate[1])
  expect_equal(fit$ATE$by_trial$estimate, tt$trial_estimates$estimate[c(1, 3)])
  independent <- do.call(CIMA, c(two_trials(), list(method = "MVLP", q_target_shared = FALSE)))
  expect_equal(vcov(independent, "trial"), Reduce(`+`, independent$target$trial_components), ignore_attr = TRUE)
  expect_equal(vcov(independent, "target"), Reduce(`+`, independent$target$components), ignore_attr = TRUE)
  expect_gt(abs(vcov(independent, "trial")[1, 2]), 0)
})

test_that("a single trial makes MVLP identical to the common fit", {
  one <- fit_small(method = "MVLP", target_subgroups = subgroup_x1)
  common <- fit_small(target_subgroups = subgroup_x1)
  expect_equal(one$target$estimates, common$target$estimates, tolerance = 1e-12)
  expect_equal(vcov(one, "target"), vcov(common, "target"), tolerance = 1e-12)
  expect_equal(unname(one$pooling$weights), matrix(1, 1, 2))
  expect_equal(unname(coef(one)[, 1]), unname(coef(common)), tolerance = 1e-12)
})

test_that("identification is checked trial by trial", {
  a <- two_trials()
  a$data_agd <- a$data_agd[c(1:5, 6:8), ]
  common <- do.call(CIMA, a)
  expect_equal(common$diagnostics$identification$rank, 3)
  err <- tryCatch(do.call(CIMA, c(a, list(method = "MVLP"))), error = identity)
  expect_s3_class(err, "cima_singularity")
  expect_match(conditionMessage(err), "trial s2 is singular: rank 2 of 3")
  a$formula_cate <- ~ x
  fit <- do.call(CIMA, c(a, list(method = "MVLP")))
  expect_equal(fit$diagnostics$identification$trial, c("s1", "s2"))
  expect_equal(fit$diagnostics$identification$rank, c(2, 2))
  expect_equal(names(fit$diagnostics$J_m_theta), c("s1", "s2"))
  expect_output(print(fit), "Moment Jacobian rank by trial: s1 2/2, s2 2/2")
  regularized <- suppressWarnings(do.call(CIMA, c(within(a, formula_cate <- ~ x + z),
    list(method = "MVLP", singularity_method = "minimum-norm"))))
  expect_true(regularized$regularization$used)
  expect_true(regularized$regularization$CATE$by_trial$s2$used)
  expect_false(regularized$regularization$CATE$by_trial$s1$used)
  expect_true(all(is.na(vcov(regularized, "target"))))
  expect_true(all(is.finite(regularized$target$sensitivity_vcov)))
})

test_that("two-step MVLP inverts each trial's own moment covariance block", {
  a <- two_trials()
  pilot <- do.call(CIMA, c(a, list(method = "MVLP")))
  fit <- do.call(CIMA, c(a, list(method = "MVLP", gmm_method = "two-step")))
  C <- pilot$variance$moment_covariance
  for (r in list(1:5, 6:10)) {
    expect_equal(unname(fit$moments$W[r, r] %*% C[r, r]), diag(5), tolerance = 1e-6)
    D <- fit$CATE$J_m_theta[r, ]; W <- fit$moments$W[r, r]
    s <- unique(fit$data$trial[r])
    expected <- solve(crossprod(D, W %*% D), crossprod(D, W %*% fit$data$TE[r]))
    expect_equal(unname(coef(fit)[, s]), as.vector(expected), tolerance = 1e-8)
  }
  expect_equal(fit$moments$W[1:5, 6:10], matrix(0, 5, 5))
  expect_equal(names(fit$gmm$optimal$by_trial), c("s1", "s2"))
})

test_that("point-only MVLP keeps the pooling weights but reports no covariance", {
  a <- two_trials()
  full <- do.call(CIMA, c(a, list(method = "MVLP")))
  point <- do.call(CIMA, c(a, list(method = "MVLP", variance = FALSE)))
  expect_equal(point$target$estimates$estimate, full$target$estimates$estimate)
  expect_equal(point$pooling$weights, full$pooling$weights)
  expect_true(all(is.na(vcov(point, "target"))))
  expect_true(all(is.na(vcov(point, "trial"))))
  expect_null(point$variance)
  expect_equal(point$inference_status, "Variance not requested")
})

test_that("MVLP methods return trial-indexed coefficients, predictions and covariances", {
  fit <- do.call(CIMA, c(two_trials(), list(method = "MVLP")))
  newdata <- data.frame(x = 0:1, z = 0:1)
  prediction <- predict(fit, newdata)
  expect_equal(dim(prediction), c(2L, 2L))
  expect_equal(colnames(prediction), c("s1", "s2"))
  expect_equal(prediction[, "s2"], as.vector(cbind(1, 0:1, 0:1) %*% coef(fit)[, "s2"]))
  expect_equal(dim(vcov(fit)), c(6L, 6L))
  expect_equal(dim(vcov(fit, "trial")), c(2L, 2L))
  expect_output(print(fit), "Pooled target estimates")
  expect_output(print(fit), "Trial-specific transported estimates")
  expect_error(vcov(fit_small(), "trial"), "only for method='MVLP'")
  expect_equal(fit_small()$method, "common")
})

test_that("continuous and time-indexed MVLP fits pool each horizon", {
  fit <- do.call(CIMA, c(two_trials("continuous"), list(method = "MVLP")))
  expect_equal(nrow(fit$target$trial_estimates), 2L)
  expect_true(all(fit$target$trial_estimates$se > 0))
  a <- small_example("continuous")
  a$data_agd$Time <- 1
  second <- transform(a$data_agd, Time = 2, mean1 = mean1 + 1)
  s1 <- rbind(a$data_agd, second)
  s2 <- transform(s1, trial = "s2", mean1 = mean1 + .5)
  a$data_agd <- rbind(s1, s2); a$mu$s2 <- a$mu$s1
  a$formula_cate <- ~ factor(Time) * (x + z)
  Sigma <- diag(14); Sigma[11, 11] <- 0
  Sigma[1:5, 6:10] <- diag(.3, 5); Sigma[6:10, 1:5] <- diag(.3, 5)
  fit <- do.call(CIMA, c(a, list(time_indexed = TRUE, method = "MVLP", Sigma = list(s1 = Sigma, s2 = Sigma))))
  expect_equal(fit$target$estimates$Time, c(1, 2))
  expect_equal(nrow(fit$target$trial_estimates), 4L)
  expect_equal(fit$target$trial_estimates$estimate, c(1.15, 2.15, 1.65, 2.65), tolerance = 1e-10)
  expect_equal(dim(fit$pooling$weights), c(2L, 2L))
  expect_gt(abs(vcov(fit, "target")[1, 2]), 0)
})

test_that("copula bootstrap MVLP replaces analytic Q variation before pooling", {
  a <- two_trials(); a$Q <- a$target <- NULL; a$summaryX <- small_copula(700)
  fit <- do.call(CIMA, c(a, list(method = "MVLP", q_method = "bootstrap",
    bootstrap_reps = 10, bootstrap_mc_size = 1500, seed = 5)))
  expect_equal(fit$copula$bootstrap$successful, 10L)
  expect_equal(ncol(fit$copula$bootstrap$draws), 6L + 2L)
  expect_equal(vcov(fit, "trial"), Reduce(`+`, fit$target$trial_components), ignore_attr = TRUE)
  L <- fit$pooling$map
  expect_equal(vcov(fit, "target"), crossprod(L, vcov(fit, "trial") %*% L), ignore_attr = TRUE)
  expect_equal(vcov(fit), fit$variance$trial_covariance + fit$copula$bootstrap$theta_covariance, ignore_attr = TRUE)
})
