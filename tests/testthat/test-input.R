test_that("binary and continuous reports are normalized without inventing summaries", {
  for (outcome in c("binary", "continuous")) {
    x <- small_example(outcome)$data_agd
    d <- prepare_agd(x)
    expect_equal(d$n, d$n1 + d$n2)
    expect_equal(d$TE, d$nu1 - d$nu0)
    expect_equal(d$seTE, sqrt(d$se_nu1^2 + d$se_nu0^2))
    expect_true(all(d$Time == 1))
    expect_true(all(d$comparison == "A vs B"))
    x$TE <- d$TE + .01; x$TE[2] <- NA
    x$seTE <- rep(.03, nrow(x)); x$seTE[3] <- NA
    e <- prepare_agd(x)
    expect_equal(e$TE[1], d$TE[1] + .01)
    expect_equal(e$TE[2], d$TE[2])
    expect_equal(e$seTE[1], .03)
    expect_equal(e$seTE[3], d$seTE[3])
  }
})

test_that("input validation rejects insufficient or inconsistent reports", {
  x <- small_example()$data_agd
  expect_error(prepare_agd(x[, setdiff(names(x), "event1")]), "Incomplete")
  bad <- x; bad$n1[1] <- -1; expect_error(prepare_agd(bad), "positive integers")
  bad <- x; bad$event1[1] <- x$n1[1] + 1; expect_error(prepare_agd(bad), "Event counts")
  bad <- x; bad$n <- x$n1 + x$n2 + 1; expect_error(prepare_agd(bad), "n must equal")
  bad <- x; bad$comparison <- c("B vs A", rep("A vs B", 4)); expect_error(prepare_agd(bad), "one consistently")
  bad <- x; bad$subgroup[1] <- NA; expect_error(prepare_agd(bad), "Missing values")
  expect_error(prepare_agd(rbind(x, x[1, ])), "Duplicate")
  args <- small_example(); args$data_agd$stratum[2] <- "system('false')"
  expect_error(do.call(CIMA, args), "Unsupported subgroup")
  args <- small_example(); args$data_agd$subgroup[2] <- "missing"
  expect_error(do.call(CIMA, args), "missing from Q")
})

test_that("sample sizes use an overall row, verified partition or explicit information", {
  a <- small_example(); a$data_agd <- a$data_agd[-1, ]
  fit <- do.call(CIMA, a)
  expect_equal(fit$trial_info$n, 20000)
  a$data_agd <- a$data_agd[c(1, 3), ]
  expect_error(do.call(CIMA, a), "Cannot identify total")
  a$trial_info <- data.frame(trial = "s1", n = 20000, allocation = .5)
  a$formula_cate <- ~ x
  expect_s3_class(do.call(CIMA, a), "CIMA")
})

test_that("legacy names migrate but effect-only covariance is rejected", {
  a <- small_example(); names(a)[c(1, 5, 6)] <- c("dataAgD", "formulaCATE", "formulaTilting")
  expect_s3_class(do.call(CIMA, a), "CIMA")
  expect_error(do.call(CIMA, c(a, list(V = list(s1 = diag(5))))), "complete joint")
})
