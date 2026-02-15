if (!exists("CIMA", mode = "function")) {
  candidates <- c("R", "../R", "../../R")
  r_dir <- candidates[file.exists(candidates)][1]
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  for (f in r_files) source(f)
}

prepare_dataAgD_real <- function(path_real) {
  if (!file.exists(path_real)) {
    return(NULL)
  }

  dataAgD <- utils::read.csv(path_real)
  dataAgD$TE <- dataAgD$event_trt / dataAgD$n_trt - dataAgD$event_c / dataAgD$n_c
  dataAgD$seTE <- sqrt(
    dataAgD$event_trt * (dataAgD$n_trt - dataAgD$event_trt) / dataAgD$n_trt^3 +
      dataAgD$event_c * (dataAgD$n_c - dataAgD$event_c) / dataAgD$n_c^3
  )
  dataAgD$n <- dataAgD$n_trt + dataAgD$n_c
  dataAgD <- dataAgD[, c("trial", "subgroup", "stratum", "TE", "seTE", "n")]
  dataAgD <- dataAgD[dataAgD$subgroup != "NYHA" & dataAgD$subgroup != "BMI", ]

  dataAgD$prop <- NA_real_
  for (s in unique(dataAgD$trial)) {
    idx <- dataAgD$trial == s
    n0 <- dataAgD$n[idx & dataAgD$subgroup == "overall"][1]
    dataAgD$prop[idx] <- dataAgD$n[idx] / n0
  }
  dataAgD
}

build_fallback_real_like <- function() {
  trials <- c("EMPEROR-Preserved", "DELIVER", "DAPA-HF", "EMPEROR-Reduced")
  rows <- list()
  set.seed(8)
  for (s in trials) {
    rows[[length(rows) + 1L]] <- data.frame(
      trial = s, subgroup = "overall", stratum = "TRUE",
      TE = runif(1, -0.09, -0.03), seTE = runif(1, 0.008, 0.018), n = sample(2500:6500, 1)
    )
    for (sg in c("LVEF", "preHHF", "diabetes")) {
      for (st in c("x == 0", "x == 1")) {
        rows[[length(rows) + 1L]] <- data.frame(
          trial = s, subgroup = sg, stratum = st,
          TE = runif(1, -0.12, 0.02), seTE = runif(1, 0.01, 0.03), n = sample(1200:5000, 1)
        )
      }
    }
  }
  dataAgD <- do.call(rbind, rows)
  dataAgD$prop <- NA_real_
  for (s in unique(dataAgD$trial)) {
    idx <- dataAgD$trial == s
    n0 <- dataAgD$n[idx & dataAgD$subgroup == "overall"][1]
    dataAgD$prop[idx] <- pmax(pmin(dataAgD$n[idx] / n0, 1), 0.1)
  }
  dataAgD
}

test_that("CIMA runs on real HF setup or fallback real-like setup", {
  skip_if_not_installed("copula")
  skip_if_not_installed("randtoolbox")

  mu <- list(
    "EMPEROR-Preserved" = c(54.3, 54.3^2 + 8.8^2, 0.229, 0.491),
    "DELIVER" = c(54.3, 54.3^2 + 8.9^2, 0.405, 0.448),
    "DAPA-HF" = c(31.1, 31.1^2 + 6.8^2, 0.474, 0.451),
    "EMPEROR-Reduced" = c(27.2, 27.2^2 + 6.1^2, 0.309, 0.498)
  )

  path_real <- "/Users/qingys/Qingys Lab/Projects/Chapter 2 CIMAgD/Project 2.1 Main Method/Part 3 HF example/dataAgD.csv"
  dataAgD <- prepare_dataAgD_real(path_real)
  if (is.null(dataAgD)) {
    dataAgD <- build_fallback_real_like()
  }
  expect_setequal(names(dataAgD), c("trial", "subgroup", "stratum", "TE", "seTE", "n", "prop"))
  expect_false(anyNA(dataAgD[, c("trial", "subgroup", "stratum", "TE", "seTE")]))

  summaryX_target <- list(
    trial = "PULSE",
    n = 50000,
    covariate = c("LVEF", "preHHF", "diabetes"),
    distribution = c("normal", "binomial", "binomial"),
    sumstats = list(
      mean = c(45.4, 0.091, 0.265),
      param2 = c(14.0, 1, 1),
      min = c(NA, NA, NA),
      max = c(NA, NA, NA),
      corr = c(-0.22, -0.046, 0.7)
    )
  )
  set.seed(20260208)
  target <- XCopula(summaryX_target, N = summaryX_target$n)
  Q <- target[sample.int(nrow(target), floor(nrow(target) / 10)), ]

  fit <- CIMA(
    dataAgD = dataAgD, mu = mu, Q = Q, target = target,
    formulaCATE = "~ LVEF + preHHF + diabetes",
    formulaTilting = "~ LVEF + I(LVEF^2) + preHHF + diabetes"
  )

  designMat <- stats::model.matrix(~ LVEF + preHHF + diabetes, data = Q)
  tiltingMat <- stats::model.matrix(~ LVEF + I(LVEF^2) + preHHF + diabetes, data = Q)
  IMat <- IMatFun(Q = Q, dataAgD = dataAgD)
  tiltingPars <- tiltingModel(mu = mu, dataAgD = dataAgD, tiltingMat = tiltingMat)
  weightsMat <- weightsMatFun(tiltingMat = tiltingMat, dataAgD = dataAgD, tiltingPars = tiltingPars)
  W <- WMat(dataAgD = dataAgD)
  foc <- as.numeric(crossprod(fit$CATE$Dtheta, W %*% (fit$moments$m_hat - fit$moments$tau_hat)))

  tilt_target_err <- vapply(names(mu), function(s) {
    mu_s <- mu[[s]]
    mu_plus <- if (length(mu_s) == ncol(tiltingMat) - 1L) c(1, mu_s) else mu_s
    w_s <- weightsMat[[s]][, 1]
    max(abs(colMeans(tiltingMat * w_s) - mu_plus))
  }, FUN.VALUE = numeric(1))

  expect_true(is.finite(fit$ATE$ATE))
  expect_true(is.finite(fit$ATE$se))
  expect_true(fit$ATE$se > 0)
  expect_true(all(is.finite(fit$CATE$theta)))
  expect_true(all(diag(fit$CATE$vcov) > 0))
  expect_lt(max(abs(foc)), 1e-6)
  expect_lt(max(tilt_target_err), 1e-5)

  # Order-invariance check for user-provided dataAgD row order.
  set.seed(20260209)
  perm <- sample.int(nrow(dataAgD))
  dataAgD_perm <- dataAgD[perm, , drop = FALSE]
  fit_perm <- CIMA(
    dataAgD = dataAgD_perm, mu = mu, Q = Q, target = target,
    formulaCATE = "~ LVEF + preHHF + diabetes",
    formulaTilting = "~ LVEF + I(LVEF^2) + preHHF + diabetes"
  )
  expect_equal(unname(fit_perm$CATE$theta), unname(fit$CATE$theta), tolerance = 1e-8)
  expect_equal(as.numeric(fit_perm$ATE$ATE), as.numeric(fit$ATE$ATE), tolerance = 1e-8)
})
