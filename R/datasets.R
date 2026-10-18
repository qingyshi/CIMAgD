#' Heart failure aggregate arm reports
#'
#' Treatment and control event counts with additive TE, binomial seTE, total n
#' and Time=1. All 42 reports are retained. Time=1 is an illustrative coding
#' convention, not verification of a common clinical outcome horizon. NYHA and
#' BMI are absent from the accompanying PULSE summaries; omit those subgroup
#' reports when fitting a model using only the PULSE covariates.
#' @format A data frame with 42 rows and 12 columns.
#' @source Aggregate trial reports supplied with the paper's heart-failure example.
"dataHF"

#' Single-source PULSE copula summaries
#'
#' Illustrative point estimates: LVEF mean 45.4, SD 14; history of heart-failure
#' hospitalization probability .091; diabetes probability .265; latent Gaussian
#' correlations -.22, -.046, .7. The source size n=50000 is an illustrative
#' simulation size and does not represent the clinical PULSE cohort size.
#' @format A list with label, n, marginals, correlation and correlation_scale.
#' @source Covariate summaries supplied with the paper's PULSE example.
"sumstatsPULSE"
