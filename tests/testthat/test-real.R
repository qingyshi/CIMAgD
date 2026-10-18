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
  mu <- list(
    "EMPEROR-Preserved" = c(54.3, 54.3^2 + 8.8^2, .229, .491),
    "DELIVER" = c(54.3, 54.3^2 + 8.9^2, .405, .448),
    "DAPA-HF" = c(31.1, 31.1^2 + 6.8^2, .474, .451),
    "EMPEROR-Reduced" = c(27.2, 27.2^2 + 6.1^2, .309, .498))
  Q <- XCopula(e$sumstatsPULSE, N = 4000, seed = 20261016)
  fit <- CIMA(reports, mu, Q, target = Q,
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
})
