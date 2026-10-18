test_that("diagnostic switches and conditional solver defaults are independent", {
  f <- fit_small()
  expect_equal(f$regularization$method, "none")
  expect_null(f$diagnostics$sufficiency)
  expect_false(f$regularization$used)
  quiet <- fit_small(jacobian_diagnostics = FALSE)
  expect_equal(quiet$regularization$method, "minimum-norm")
  expect_null(quiet$diagnostics$J_m_theta)
  expect_null(quiet$diagnostics$H_s)
  expect_equal(coef(quiet), coef(f), tolerance = 1e-12)
  expect_equal(vcov(quiet), vcov(f), tolerance = 1e-12)
  span <- fit_small(jacobian_diagnostics = FALSE, tilting_diagnostics = TRUE)
  expect_true(span$diagnostics$sufficiency$s1[[1]]$sufficient)
  expect_output(print(quiet), "GMM: one-step")
  a <- small_example(); a$formula_cate <- ~ x + z + I(2*x)
  expect_error(do.call(CIMA, c(a, list(jacobian_diagnostics=FALSE, singularity_method="none"))),
    "singular", class="cima_singularity")
})

test_that("minimum-norm returns the original-coordinate MP coefficient solution", {
  a <- small_example(); a$formula_cate <- ~ x + z + I(2*x)
  expect_warning(f <- do.call(CIMA, c(a, list(jacobian_diagnostics=FALSE))), "Wald")
  expect_equal(unname(coef(f)), c(.1,.01,-.02,.02), tolerance=1e-9)
  A <- chol(f$moments$W) %*% f$CATE$J_m_theta
  y <- chol(f$moments$W) %*% f$data$TE
  s <- svd(A); keep <- s$d > max(s$d)*1e-9
  oracle <- s$v[,keep,drop=FALSE] %*% (crossprod(s$u[,keep,drop=FALSE],y)/s$d[keep])
  expect_equal(unname(coef(f)), as.vector(oracle), tolerance=1e-9)
  expect_true(all(is.na(vcov(f))))
  expect_true(all(is.na(f$CATE$ci95)))
  expect_true(all(is.finite(f$variance$sensitivity_vcov_theta)))
  expect_true(f$regularization$used)
  expect_equal(f$moments$m_hat, f$data$TE, tolerance=1e-9)
})

test_that("ridge point solutions are equivariant to predictor units", {
  a <- small_example(); a$formula_cate <- ~ x + z + I(2*x)
  f <- do.call(CIMA,c(a,list(singularity_method="ridge",variance=FALSE)))
  a$Q$big <- a$Q$x*1e8; a$target <- a$Q
  a$formula_cate <- ~ big + z + I(2*big)
  rescaled <- do.call(CIMA,c(a,list(singularity_method="ridge",variance=FALSE)))
  expect_equal(predict(rescaled,a$Q),predict(f,a$Q),tolerance=1e-7)
  expect_equal(rescaled$regularization$CATE$ridge,1e-8)
  expect_equal(unname(coef(rescaled))*c(1,1e8,1,1e8),unname(coef(f)),tolerance=1e-7)
  A <- chol(f$moments$W)%*%f$CATE$J_m_theta
  norms <- sqrt(colSums(A^2)); z <- sweep(A,2,norms,"/")
  expected <- solve(crossprod(z)+diag(1e-8,4),crossprod(z,chol(f$moments$W)%*%f$data$TE))/norms
  expect_equal(unname(coef(f)),as.vector(expected),tolerance=1e-7)
})

test_that("redundant tilting functions preserve feasible calibration and reject conflicts", {
  a <- small_example()
  a$formula_tilting <- ~ x*z + I(2*x)
  a$mu <- list(s1=c(x=.6,z=.7,"I(2 * x)"=1.2,"x:z"=.4))
  baseline <- fit_small()
  for (method in c("minimum-norm","ridge")) {
    fit <- do.call(CIMA,c(a,list(singularity_method=method,variance=FALSE)))
    expect_equal(coef(fit),coef(baseline),tolerance=1e-8)
    expect_lt(max(abs(fit$tilting$s1$moment_error)),1e-8)
    expect_equal(as.vector(exp(fit$tilting$s1$h_plus%*%fit$tilting$s1$eta)),fit$tilting$s1$weights,tolerance=1e-9)
    expect_true(fit$tilting$s1$regularization$used)
    conflicting <- a; conflicting$mu$s1["I(2 * x)"] <- 1.3
    expect_error(do.call(CIMA,c(conflicting,list(singularity_method=method))),"Incompatible redundant")
  }
  fit <- do.call(CIMA,c(a,list(singularity_method="minimum-norm",variance=FALSE)))
  H <- fit$tilting$s1$H_s; P <- fit$tilting$s1$H_inv
  expect_equal(H%*%P%*%H,H,tolerance=1e-8)
  expect_equal(P%*%H%*%P,P,tolerance=1e-8)
  expect_equal(H%*%P,t(H%*%P),tolerance=1e-8)
})

test_that("constant redundant tilt terms and singular H ridge handle unit changes", {
  a <- small_example(); a$formula_tilting <- ~ x*z + I(1+0*x)
  a$mu <- list(s1=c(x=.6,z=.7,"I(1 + 0 * x)"=1,"x:z"=.4))
  fit <- do.call(CIMA,c(a,list(singularity_method="minimum-norm",variance=FALSE)))
  expect_equal(coef(fit),coef(fit_small()),tolerance=1e-9)
  H <- matrix(c(2,2,2,2),2)
  ctl <- list(rank_tol=1e-9,condition_limit=1e8,singularity_method="ridge",ridge_strength=1e-8,jacobian_diagnostics=FALSE)
  base <- symmetric_inverse(H,"H",ctl)
  unit <- diag(c(1e6,1e-3)); changed <- symmetric_inverse(unit%*%H%*%unit,"H",ctl)
  expect_equal(unit%*%changed$inverse%*%unit,base$inverse,tolerance=1e-6)
  expect_equal(changed$ridge,base$ridge)
})
