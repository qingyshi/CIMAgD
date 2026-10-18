test_that("time indexing is opt-in and separate horizons can be identified", {
  a <- small_example("continuous")
  second <- transform(a$data_agd, trial = "s2", Time = 2, mean1 = mean1 + 1)
  a$data_agd$Time <- 1
  a$data_agd <- rbind(a$data_agd, second); a$mu$s2 <- a$mu$s1
  pooled <- do.call(CIMA, a)
  expect_false(pooled$time_indexed)
  expect_equal(nrow(pooled$target$estimates), 1L)
  a$formula_cate <- ~ factor(Time) * (x + z)
  expect_error(do.call(CIMA, a), "time_indexed=TRUE")
  fit <- do.call(CIMA, c(a, list(time_indexed = TRUE)))
  expect_equal(fit$ATE$ATE, c(1.15, 2.15), tolerance = 1e-10)
  expect_equal(dim(vcov(fit, "target")), c(2L, 2L))
  expect_gt(abs(vcov(fit, "target")[1,2]), 0)
  expect_equal(predict(fit, data.frame(x=0,z=0), time=2), 2, tolerance = 1e-10)
  expect_error(predict(fit, data.frame(x=0,z=0), time=3), "observed time")
})

test_that("repeated horizons require joint information and support complete Sigma", {
  a <- small_example("continuous")
  a$data_agd$Time <- 1
  second <- transform(a$data_agd, Time = 2, mean1 = mean1 + 1)
  a$data_agd <- rbind(a$data_agd, second)
  a$formula_cate <- ~ factor(Time) * (x + z)
  a$time_indexed <- TRUE
  expect_error(do.call(CIMA, a), "cross-time covariance")
  point <- do.call(CIMA, c(a, list(variance = FALSE)))
  expect_equal(point$ATE$ATE, c(1.15, 2.15), tolerance = 1e-10)
  Sigma <- diag(14); Sigma[11,11] <- 0
  Sigma[1:5,6:10] <- diag(.3,5); Sigma[6:10,1:5] <- diag(.3,5)
  fitted <- do.call(CIMA, c(a, list(Sigma = list(s1 = Sigma))))
  expect_true(all(is.finite(vcov(fitted))))
  expect_gt(vcov(fitted, "target")[1,2], 0)
})
