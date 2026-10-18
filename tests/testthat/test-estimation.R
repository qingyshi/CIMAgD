test_that("known linear effects and normalized moments are recovered", {
  for (outcome in c("binary", "continuous")) {
    fit <- fit_small(outcome, tilting_diagnostics = TRUE)
    truth <- if (outcome == "binary") c(.1, .05, -.02) else c(1, .5, -.2)
    expect_equal(unname(coef(fit)), truth, tolerance = 1e-10)
    expect_equal(colMeans(fit$moments$alpha), rep(1, 5), ignore_attr = TRUE, tolerance = 1e-12)
    expect_lt(max(abs(fit$moments$foc)), 1e-9)
    expect_lt(max(abs(fit$tilting$s1$moment_error)), 1e-8)
    expect_equal(predict(fit, data.frame(x = 0:1, z = 0:1)), c(truth[1], sum(truth)))
    expect_true(all(diag(vcov(fit)) > 0))
    expect_true(all(fit$diagnostics$sufficiency$s1[[1]]$columns$in_span))
    expect_lt(max(abs(fit$variance$Gamma_q)), 1e-16)
  }
})

test_that("ratio normalizers and results are invariant to order", {
  a <- small_example(saturated = FALSE)
  a$data_agd$n1[2] <- a$data_agd$n1[2] + 20
  fit <- do.call(CIMA, a)
  expect_equal(colMeans(fit$moments$alpha), rep(1, 5), ignore_attr = TRUE)
  expect_gt(abs(fit$moments$c_sj[2] - fit$data$n[2] / fit$trial_info$n), 1e-5)
  a <- small_example(); fit <- do.call(CIMA, a)
  a$data_agd <- a$data_agd[c(5, 1, 3, 2, 4), ]
  reordered <- do.call(CIMA, a)
  expect_equal(coef(reordered), coef(fit), tolerance = 1e-10)
  expect_equal(vcov(reordered), vcov(fit), tolerance = 1e-10)
})

test_that("rank and support failures carry useful diagnostics", {
  a <- small_example(); a$formula_cate <- ~ x + z + I(2*x)
  err <- tryCatch(do.call(CIMA, a), error = identity)
  expect_s3_class(err, "cima_singularity")
  expect_equal(err$matrix_name, "J_m_theta")
  expect_false(err$diagnostics$full_rank)
  a <- small_example(); a$formula_tilting <- ~ x + z + I(2*x); a$mu <- list(s1 = c(.6, .7, 1.2))
  expect_error(do.call(CIMA, a), "H_s cannot be invertible", class = "cima_singularity")
  a <- small_example(); a$mu$s1[1] <- 1.1
  expect_error(do.call(CIMA, a), "outside Q support")
  a <- small_example(); a$data_agd$stratum[2] <- "x > 20"
  expect_error(do.call(CIMA, a), "No Q support")
  dg <- diagnose_jacobian(cbind(1:4, 1:4 + c(0, 0, 0, 1e-5)), condition_limit = 100)
  expect_equal(dg$status, "ill_conditioned")
  expect_equal(diagnose_jacobian(matrix(1:6, 2, 3))$rank, 2)
})

test_that("A.4.1 projection identifies the missing interaction", {
  q <- expand.grid(x = 0:1, z = 0:1)
  i <- cbind(overall = 1, x1 = q$x)
  insufficient <- diagnose_tilting(model.matrix(~ x + z, q), model.matrix(~ x + z, q), i)
  sufficient <- diagnose_tilting(model.matrix(~ x * z, q), model.matrix(~ x + z, q), i)
  expect_false(insufficient$sufficient)
  expect_true(sufficient$sufficient)
  expect_true(any(grepl("I\\*z", insufficient$unresolved)))
  expect_equal(insufficient$R_s, crossprod(insufficient$residuals) / nrow(q))
  expect_gte(min(eigen(insufficient$R_s, symmetric = TRUE)$values), -1e-12)
  redundant <- diagnose_tilting(model.matrix(~ x * z + I(2*x), q), model.matrix(~x+z, q), i)
  expect_true(redundant$sufficient)
})

test_that("trial-specific tilting and explicit proportion matching work", {
  a <- small_example(); a$formula_tilting <- list(s1 = ~ x * z)
  a$mu <- list(s1 = c("z" = .7, "x:z" = .4, "x" = .6))
  fit <- do.call(CIMA, c(a, list(match_subgroup_proportions = TRUE)))
  expect_equal(unname(coef(fit)), c(.1, .05, -.02), tolerance = 1e-9)
  a$mu$s1["x"] <- .61
  expect_error(do.call(CIMA, c(a, list(match_subgroup_proportions = TRUE))), "conflicts")
})

test_that("rescaling a covariate preserves estimates and control fits", {
  a <- small_example(); base <- do.call(CIMA,a)
  a$Q$large <- a$Q$x * 1e10
  a$target <- a$Q
  a$formula_cate <- ~ large + z
  fit <- do.call(CIMA,a)
  expect_equal(unname(coef(fit))*c(1,1e10,1),unname(coef(base)),tolerance=1e-9)
  expect_equal(fit$ATE$ATE,base$ATE$ATE,tolerance=1e-10)
  expect_equal(fit$ATE$var,base$ATE$var,tolerance=1e-10)
  expect_equal(fit$variance$trials$s1$control$gamma0,
               base$variance$trials$s1$control$gamma0,tolerance=1e-9)
})
