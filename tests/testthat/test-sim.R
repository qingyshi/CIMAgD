if (!exists("CIMA", mode = "function")) {
  candidates <- c("R", "../R", "../../R")
  r_dir <- candidates[file.exists(candidates)][1]
  r_files <- list.files(r_dir, pattern = "[.]R$", full.names = TRUE)
  for (f in r_files) source(f)
}

simulate_data_binary <- function(N, eta, beta, gamma, theta) {
  covariates <- as.data.frame(simstudy::genCorGen(
    N, nvars = 3, params1 = eta[1:3], dist = "binary", rho = eta[4],
    corstr = "cs", wide = TRUE
  ))
  colnames(covariates) <- c("id", "X1", "X2", "X3")

  pr_r <- plogis(
    beta[1] + beta[2] * covariates$X1 + beta[3] * covariates$X2 + beta[4] * covariates$X3
  )
  R <- stats::rbinom(N, 1, pr_r)
  data_all <- cbind(R, covariates)
  data_trial <- data_all[data_all$R == 1, ]
  data_target <- data_all[data_all$R == 0, ]

  pr_s <- t(mapply(
    function(x1, x2, x3) {
      lp <- c(
        0,
        gamma[1, 1] + gamma[1, 2] * x1 + gamma[1, 3] * x2 + gamma[1, 4] * x3,
        gamma[2, 1] + gamma[2, 2] * x1 + gamma[2, 3] * x2 + gamma[2, 4] * x3,
        gamma[3, 1] + gamma[3, 2] * x1 + gamma[3, 3] * x2 + gamma[3, 4] * x3,
        gamma[4, 1] + gamma[4, 2] * x1 + gamma[4, 3] * x2 + gamma[4, 4] * x3
      )
      p <- exp(lp) / sum(exp(lp))
      names(p) <- paste0("S", 1:5)
      p
    },
    data_trial$X1, data_trial$X2, data_trial$X3
  ))
  S <- apply(pr_s, 1, function(p) sample(paste0("S", 1:5), 1, prob = p))
  data_trial <- cbind(data_trial, S)
  data_target$S <- NA
  data_all <- rbind(data_trial, data_target)

  data_all$A <- ave(
    rep(1, nrow(data_all)),
    data_all$R,
    data_all$S,
    FUN = function(x) stats::rbinom(length(x), 1, 0.5)
  )

  data_all$`Y(1)` <- stats::rbinom(
    nrow(data_all), 1,
    plogis(theta[1, 1] + theta[1, 2] * data_all$X1 + theta[1, 3] * data_all$X2 + theta[1, 4] * data_all$X3)
  )
  data_all$`Y(0)` <- stats::rbinom(
    nrow(data_all), 1,
    plogis(theta[2, 1] + theta[2, 2] * data_all$X1 + theta[2, 3] * data_all$X2 + theta[2, 4] * data_all$X3)
  )
  data_all$Y <- data_all$A * data_all$`Y(1)` + (1 - data_all$A) * data_all$`Y(0)`
  data_all
}

make_agd <- function(dt_trial) {
  overall <- data.frame(
    trial = paste0("S", 1:5),
    subgroup = "overall",
    stratum = "TRUE"
  )
  subgroup <- expand.grid(
    trial = paste0("S", 1:5),
    subgroup = c("X1", "X2", "X3"),
    stratum = c("x == 0", "x == 1")
  )
  dataAgD <- rbind(overall, subgroup)

  agd <- mapply(
    function(trial_i, subgroup_i, stratum_i) {
      trial_ind <- dt_trial$S == trial_i
      if (subgroup_i == "overall" || stratum_i == "TRUE") {
        subgroup_ind <- rep(TRUE, sum(trial_ind))
      } else {
        subgroup_ind <- eval(parse(text = paste0("function (x) ", stratum_i)))(
          dt_trial[trial_ind, subgroup_i]
        )
      }
      dt_sub <- dt_trial[trial_ind, ][subgroup_ind, ]
      p_t <- mean(dt_sub$Y[dt_sub$A == 1])
      p_c <- mean(dt_sub$Y[dt_sub$A == 0])
      n_t <- sum(dt_sub$A == 1)
      n_c <- sum(dt_sub$A == 0)
      c(
        TE = p_t - p_c,
        seTE = sqrt(p_t * (1 - p_t) / n_t + p_c * (1 - p_c) / n_c),
        N = n_t + n_c
      )
    },
    dataAgD$trial, dataAgD$subgroup, dataAgD$stratum
  )
  agd <- t(agd)
  dataAgD$TE <- agd[, "TE"]
  dataAgD$seTE <- agd[, "seTE"]
  dataAgD$n <- agd[, "N"]
  dataAgD <- dataAgD[dataAgD$seTE > 0, ]

  dataAgD$prop <- NA_real_
  for (s in unique(dataAgD$trial)) {
    idx <- dataAgD$trial == s
    n0 <- dataAgD$n[idx & dataAgD$subgroup == "overall"][1]
    dataAgD$prop[idx] <- dataAgD$n[idx] / n0
  }
  dataAgD
}

make_mu <- function(dt_trial) {
  out <- vector("list", 5)
  names(out) <- paste0("S", 1:5)
  for (s in names(out)) {
    d <- dt_trial[dt_trial$S == s, ]
    out[[s]] <- c(
      mean(d$X1), mean(d$X2), mean(d$X3),
      mean(d$X1 * d$X2), mean(d$X2 * d$X3),
      mean(d$X1 * d$X3), mean(d$X1 * d$X2 * d$X3)
    )
  }
  out
}

test_that("CIMA recovers ATE and CATE direction in simulation", {
  skip_if_not_installed("simstudy")
  set.seed(20260208)

  N <- 10000
  eta <- c(0.5, 0.5, 0.5, 0.5)
  beta <- c(log(2), log(0.5), log(2), log(0.5))
  gamma <- matrix(
    c(
      log(2), log(0.5), log(2), log(0.5),
      log(2), log(0.8), log(1.25), log(0.8),
      log(2), log(0.5), log(2), log(0.5),
      log(2), log(0.8), log(1.25), log(0.8)
    ),
    nrow = 4, byrow = TRUE
  )
  theta <- matrix(
    c(log(0.5), log(0.5), log(0.2), log(0.5),
      log(2), log(2), log(5), log(2)),
    nrow = 2, byrow = TRUE
  )

  dt <- simulate_data_binary(N, eta, beta, gamma, theta)
  dt_trial <- dt[dt$R == 1, ]
  dt_target <- dt[dt$R == 0, ]

  dataAgD <- make_agd(dt_trial)
  mu <- make_mu(dt_trial)

  Q <- dt_target
  fit <- CIMA(
    dataAgD = dataAgD, mu = mu, Q = Q, target = dt_target,
    formulaCATE = "~ X1 + X2 + X3",
    formulaTilting = "~ X1 + X2 + X3 + X1:X2 + X2:X3 + X1:X3 + X1:X2:X3"
  )

  ate_true <- mean(dt_target$`Y(1)` - dt_target$`Y(0)`)
  designMat <- stats::model.matrix(~ X1 + X2 + X3, data = Q)
  tiltingMat <- stats::model.matrix(~ X1 + X2 + X3 + X1:X2 + X2:X3 + X1:X3 + X1:X2:X3, data = Q)
  expect_setequal(names(dataAgD), c("trial", "subgroup", "stratum", "TE", "seTE", "n", "prop"))
  expect_false(anyNA(dataAgD[, c("trial", "subgroup", "stratum", "TE", "seTE")]))

  IMat <- IMatFun(Q = Q, dataAgD = dataAgD)
  tiltingPars <- tiltingModel(mu = mu, dataAgD = dataAgD, tiltingMat = tiltingMat)
  weightsMat <- weightsMatFun(tiltingMat = tiltingMat, dataAgD = dataAgD, tiltingPars = tiltingPars)
  W <- WMat(dataAgD = dataAgD)
  obj_hat <- gmmFun(fit$CATE$theta, W, designMat, IMat, weightsMat, dataAgD)
  obj_zero <- gmmFun(rep(0, ncol(designMat)), W, designMat, IMat, weightsMat, dataAgD)
  foc <- as.numeric(crossprod(fit$CATE$Dtheta, W %*% (fit$moments$m_hat - fit$moments$tau_hat)))

  # Tilting target check: E_Q[w_s h_s^+(X)] ~= mu_s^+
  tilt_target_err <- vapply(names(mu), function(s) {
    mu_s <- mu[[s]]
    mu_plus <- if (length(mu_s) == ncol(tiltingMat) - 1L) c(1, mu_s) else mu_s
    w_s <- weightsMat[[s]][, 1]
    max(abs(colMeans(tiltingMat * w_s) - mu_plus))
  }, FUN.VALUE = numeric(1))

  expect_true(is.finite(fit$ATE$ATE))
  expect_true(all(is.finite(fit$CATE$theta)))
  expect_true(all(diag(fit$CATE$vcov) > 0))

  expect_lt(abs(fit$ATE$ATE - ate_true), 0.1)
  expect_lt(obj_hat, obj_zero)
  expect_lt(max(abs(foc)), 1e-6)
  expect_lt(max(tilt_target_err), 1e-5)

  # Order-invariance check: user-provided row order of dataAgD should not change fit.
  set.seed(20260209)
  perm <- sample.int(nrow(dataAgD))
  dataAgD_perm <- dataAgD[perm, , drop = FALSE]
  fit_perm <- CIMA(
    dataAgD = dataAgD_perm, mu = mu, Q = Q, target = dt_target,
    formulaCATE = "~ X1 + X2 + X3",
    formulaTilting = "~ X1 + X2 + X3 + X1:X2 + X2:X3 + X1:X3 + X1:X2:X3"
  )
  expect_equal(unname(fit_perm$CATE$theta), unname(fit$CATE$theta), tolerance = 1e-8)
  expect_equal(as.numeric(fit_perm$ATE$ATE), as.numeric(fit$ATE$ATE), tolerance = 1e-8)
})
