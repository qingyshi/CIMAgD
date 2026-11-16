hf_mu <- list(
  "EMPEROR-Preserved" = c(54.3, 54.3^2 + 8.8^2, .229, .491),
  "DELIVER" = c(54.3, 54.3^2 + 8.9^2, .261, .448),
  "DAPA-HF" = c(31.1, 31.1^2 + 6.8^2, .274, .451),
  "EMPEROR-Reduced" = c(27.2, 27.2^2 + 6.1^2, .309, .498))

test_that("packaged HF/PULSE data drive the actual example without fallback", {
  e <- new.env()
  utils::data("dataHF", package = "CIMAgD", envir = e)
  utils::data("sumstatsPULSE", package = "CIMAgD", envir = e)
  expect_equal(e$dataHF, readRDS(system.file("extdata", "dataHF.rds", package = "CIMAgD")))
  expect_equal(e$sumstatsPULSE, readRDS(system.file("extdata", "sumstatsPULSE.rds", package = "CIMAgD")))
  expect_equal(nrow(e$dataHF), 42L)
  expect_setequal(names(e$dataHF), c("trial", "comparison", "subgroup", "stratum",
    "TE", "seTE", "n", "Time", "event1", "n1", "event2", "n2"))
  expect_true(all(e$dataHF$Time == 1))
  expect_equal(e$dataHF$TE, with(e$dataHF, event1/n1-event2/n2))
  # PULSE has no NYHA/BMI summaries; all source rows remain in dataHF.
  reports <- subset(e$dataHF, !subgroup %in% c("NYHA", "BMI"))
  # The binary summaries in mu are the dataHF subgroup proportions, to three decimals.
  total <- with(subset(reports, subgroup == "overall"), stats::setNames(n, trial))
  for (k in 1:2) {
    rows <- subset(reports, subgroup == c("preHHF", "diabetes")[k] & stratum == "x == 1")
    expect_equal(round(rows$n / total[rows$trial], 3),
                 vapply(hf_mu[rows$trial], `[`, numeric(1), k + 2L), ignore_attr = TRUE)
  }
  Q <- XCopula(e$sumstatsPULSE, N = 4000, seed = 20261016)
  fit <- CIMA(reports, hf_mu, Q, target = Q,
    formula_cate = ~ LVEF + preHHF + diabetes, tilting_diagnostics = TRUE,
    formula_tilting = ~ LVEF + I(LVEF^2) + preHHF + diabetes,
    target_subgroups = list(diabetes = function(d) d$diabetes == 1))
  expect_equal(fit$diagnostics$J_m_theta$rank, 4)
  expect_true(all(is.finite(coef(fit))))
  expect_true(all(is.finite(vcov(fit))))
  expect_lt(max(abs(fit$moments$foc)), 1e-6)
  expect_true(all(vapply(fit$tilting, function(s) s$standardized_error < 1e-8, logical(1))))
  expect_true(all(vapply(fit$tilting, function(s) s$diagnostics$full_rank, logical(1))))
  expect_true(any(vapply(fit$diagnostics$sufficiency, function(s) !s[[1]]$sufficient, logical(1))))
  expect_equal(nrow(fit$target$estimates), 2L)
  expect_gte(min(eigen(vcov(fit), symmetric=TRUE)$values), -1e-10)
  # Trials without LVEF subgroups cannot identify their own working control mean.
  wc <- fit$variance$working_control
  expect_equal(wc$model, c("trial-specific", "trial-specific", "common", "common"))
  expect_true(all(grepl("unidentified", wc$trial_specific_issues[3:4])))
  # MVLP needs every trial to identify its own CATE.
  err <- tryCatch(CIMA(reports, hf_mu, Q, target = Q, method = "MVLP",
    formula_cate = ~ LVEF + preHHF + diabetes,
    formula_tilting = ~ LVEF + I(LVEF^2) + preHHF + diabetes), error = identity)
  expect_s3_class(err, "cima_singularity")
  expect_match(conditionMessage(err), "trial DAPA-HF is singular: rank 3 of 4")
  pooled <- CIMA(reports, hf_mu, Q, target = Q, method = "MVLP",
    formula_cate = ~ preHHF + diabetes,
    formula_tilting = ~ LVEF + I(LVEF^2) + preHHF + diabetes)
  expect_equal(pooled$diagnostics$identification$rank, rep(3, 4))
  expect_equal(nrow(pooled$target$trial_estimates), 4L)
  expect_equal(sum(pooled$pooling$weights), 1)
  expect_true(all(is.finite(vcov(pooled, "trial"))))
})
