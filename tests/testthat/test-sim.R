simulate_trial <- function(n = 2400L, n_target = 800L) {
  grid <- expand.grid(x = 0:1, z = 0:1)
  trial <- grid[sample.int(4L, n, replace=TRUE, prob=c(.1,.2,.3,.4)), ]
  trial$A <- rbinom(n, 1, .5)
  g <- .08 + .06 * trial$x - .03 * trial$z
  control <- .35 + .04 * trial$x + .02 * trial$z
  trial$Y <- rbinom(n, 1, control + trial$A*g)
  reports <- data.frame(trial="s1", subgroup=c("overall","x","x","z","z"),
    stratum=c("TRUE","x == 0","x == 1","x == 0","x == 1"))
  indicators <- cbind(TRUE, trial$x == 0, trial$x == 1, trial$z == 0, trial$z == 1)
  for (a in 1:2) {
    arm <- trial$A == if (a==1) 1 else 0
    reports[[paste0("n",a)]] <- colSums(indicators & arm)
    reports[[paste0("event",a)]] <- colSums(indicators & arm & trial$Y==1)
  }
  target <- grid[sample.int(4L,n_target,replace=TRUE,prob=c(.4,.3,.2,.1)), ]
  mu <- c(mean(trial$x), mean(trial$z), mean(trial$x*trial$z))
  list(data_agd=reports, mu=list(s1=mu), Q=target, target=target,
    formula_cate=~x+z, formula_tilting=~x*z)
}

test_that("seeded binary simulations recover coefficients and calibrate variance", {
  set.seed(66104)
  B <- 60L
  draws <- vapply(seq_len(B), function(k) {
    args <- simulate_trial(); fit <- do.call(CIMA,args)
    c(unname(coef(fit)), estimate=fit$ATE$ATE, variance=fit$ATE$var)
  }, numeric(5))
  truth <- c(.08,.06,-.03)
  target_truth <- .08 + .06*.4 - .03*.3
  expect_lt(max(abs(rowMeans(draws[1:3,,drop=FALSE]) - truth)), .02)
  bias <- mean(draws[4,])-target_truth
  expect_lt(abs(bias), .012)
  ratio <- var(draws[4,]) / mean(draws[5,])
  expect_gt(ratio,.5); expect_lt(ratio,1.8)
  coverage <- mean(abs(draws[4,]-target_truth) <= qnorm(.975)*sqrt(draws[5,]))
  expect_gte(coverage,.83)
})
