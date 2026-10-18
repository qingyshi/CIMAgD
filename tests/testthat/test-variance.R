test_that("centered ratio Jacobians match numerical derivatives", {
  fit <- fit_small(saturated = FALSE)
  h <- fit$tilting$s1$h_plus; eta <- fit$tilting$s1$eta
  a <- small_example(saturated = FALSE); d <- prepare_agd(a$data_agd)
  ii <- indicator_matrix(a$Q, d); b <- model.matrix(a$formula_cate, a$Q)
  g <- as.vector(b %*% coef(fit)); eps <- 1e-5
  moment <- function(eta) {
    w <- exp(h %*% eta)[, 1]
    colSums(ii * (w * g)) / colSums(ii * w)
  }
  numeric_eta <- vapply(seq_along(eta), function(k) {
    delta <- rep(0, length(eta)); delta[k] <- eps
    (moment(eta + delta) - moment(eta - delta)) / (2 * eps)
  }, numeric(nrow(d)))
  expect_equal(unname(fit$variance$J_m_eta$s1), unname(numeric_eta), tolerance = 1e-8)
  expect_equal(fit$variance$J_m_mu$s1,
    fit$variance$J_m_eta$s1 %*% fit$tilting$s1$H_inv, tolerance = 1e-12)
  expect_lt(max(abs(fit$variance$J_m_eta$s1[, 1])), 1e-14)
})

test_that("A.6.2 keeps joint covariance, anchors margins and does not repair eigenvalues", {
  fit <- fit_small(saturated = FALSE)
  s <- fit$variance$trials$s1
  expect_equal(diag(s$Sigma)[1:5], fit$trial_info$n * fit$data$seTE^2, ignore_attr = TRUE)
  expect_true(any(abs(s$Sigma[1:5, 7:8]) > 1e-8))
  expect_equal(s$Sigma[6, ], rep(0, ncol(s$Sigma)), ignore_attr = TRUE)
  expect_gte(min(eigen(s$Sigma, symmetric = TRUE)$values), -1e-10)
  expect_equal(fit$variance$Gamma_s$s1, s$D_s %*% s$Sigma %*% t(s$D_s), ignore_attr = TRUE)
  expect_equal(vcov(fit), fit$variance$q_covariance + fit$variance$trial_covariance, ignore_attr = TRUE)
  expect_equal(fit$variance$Gamma_q, crossprod(fit$variance$xi_q) / nrow(fit$variance$xi_q))
  expect_equal(s$control$gamma1 - s$control$gamma0,
    as.vector(fit$basis$basis_q[[1]] %*% coef(fit)), tolerance = 1e-10)
  expect_true(all(s$control$gamma0 >= 0 & s$control$gamma1 <= 1))
  with_se <- fit_small(mu_se = list(s1 = c(.01, .02, .03)))
  expect_equal(diag(with_se$variance$trials$s1$Sigma)[7:9], 20000 * c(.01, .02, .03)^2, ignore_attr = TRUE)
})

test_that("complete Sigma propagation includes cross-covariance with correct signs", {
  a <- small_example(saturated = FALSE)
  base <- do.call(CIMA, a)
  sigma <- base$variance$trials$s1$Sigma
  fitted <- do.call(CIMA, c(a, list(Sigma = list(s1 = sigma))))
  expect_equal(vcov(fitted), vcov(base), tolerance = 1e-12)
  expect_null(fitted$variance$trials$s1$control)
  altered <- sigma; altered[1:5, 6:8] <- 0; altered[6:8, 1:5] <- 0
  other <- do.call(CIMA, c(a, list(Sigma = list(s1 = altered))))
  expect_gt(max(abs(vcov(other) - vcov(base))), 1e-8)
  bad <- sigma; bad[6, 6] <- 1
  expect_error(do.call(CIMA, c(a, list(Sigma = list(s1 = bad)))), "normalization")
  bad <- sigma; bad[1, 1] <- -1
  expect_error(do.call(CIMA, c(a, list(Sigma = list(s1 = bad)))), "positive semidefinite")
})

test_that("shared Q target and subgroup influence contributions are combined", {
  subgroups <- list(x1 = function(d) d$x == 1, joint = function(d) d$x == 1 & d$z == 0)
  fit <- fit_small(saturated = FALSE, target_subgroups = subgroups)
  tt <- fit$target; vv <- fit$variance
  zeta <- tt$phi_target + vv$theta_if_q %*% t(tt$J_psi_theta)
  expected <- crossprod(center_columns(zeta)) / nrow(zeta)^2 +
    tt$J_psi_theta %*% vv$trial_covariance %*% t(tt$J_psi_theta)
  expect_equal(vcov(fit, "target"), expected, ignore_attr = TRUE)
  expect_equal(tt$estimates$n, c(200L, 100L, 50L))
  independent <- fit_small(saturated = FALSE, q_target_shared = FALSE, target_subgroups = subgroups)
  expect_false(independent$target$shared_q)
  expect_equal(vcov(independent, "target"),
    Reduce(`+`, independent$target$components), ignore_attr = TRUE)
  a <- small_example(); a$target <- a$target[-1, ]
  expect_error(do.call(CIMA, c(a, list(q_target_shared = TRUE))), "row-aligned")
})

test_that("continuous reconstruction estimates residual variances and zero margins are explicit", {
  fit <- fit_small("continuous")
  control <- fit$variance$trials$s1$control
  expect_equal(control$residual_diagnostics$control$variance, 1, tolerance = 1e-9)
  expect_equal(control$residual_diagnostics$treatment$variance, 2, tolerance = 1e-9)
  anchor <- anchor_covariance(matrix(c(0,0,0,2), 2), c(3,0))
  expect_equal(anchor$Sigma, diag(c(3,0)))
  expect_equal(anchor$unmodeled_variance, 1L)
  a <- small_example(); a$data_agd$event1 <- a$data_agd$event2 <- 0L
  fit <- do.call(CIMA, a)
  expect_equal(unname(coef(fit)), rep(0,3))
  expect_equal(unname(vcov(fit)), matrix(0,3,3))
  expect_match(fit$moments$weight_note, "Identity")
})

test_that("Q influence combines contributions from all trials before covariance", {
  a <- small_example(saturated = FALSE)
  a$data_agd <- rbind(a$data_agd, transform(a$data_agd, trial = "s2"))
  a$mu$s2 <- c(.4, .3)
  fit <- do.call(CIMA, a)
  expected <- fit$variance$phi_q_m
  for (s in names(fit$tilting)) {
    t <- fit$tilting[[s]]
    expected <- expected - sweep(t$h_plus * t$weights, 2, t$mu_plus, "-") %*% t(fit$variance$J_m_mu[[s]])
  }
  expect_equal(fit$variance$xi_q, center_columns(expected), tolerance = 1e-12)
  expect_equal(fit$variance$Gamma_q, empirical_cov(expected), tolerance = 1e-12)
  expect_gt(max(abs(fit$variance$Gamma_q[1:5, 6:10])), 1e-8)
})

test_that("working joint covariance agrees with direct enumeration of (X,A,Y)", {
  fit <- fit_small()
  q <- expand.grid(x=0:1,z=0:1); prob <- c(.1,.2,.3,.4)
  indicators <- cbind(1,q$x==0,q$x==1,q$z==0,q$z==1)
  p <- colSums(indicators*prob); k <- sweep(indicators,2,p,"/")
  h <- model.matrix(~x*z,q); mu <- colSums(h*prob)
  g0 <- .3+.03*q$x+.02*q$z; g1 <- g0+.1+.05*q$x-.02*q$z
  nu0 <- colSums(k*(prob*g0)); nu1 <- colSums(k*(prob*g1))
  oracle <- matrix(0,9,9)
  for (i in 1:4) for (a in 0:1) for (y in 0:1) {
    risk <- if (a==1) g1[i] else g0[i]
    mass <- prob[i]*.5*if (y==1) risk else (1-risk)
    phi <- c(k[i,]*(a/.5*(y-nu1)-(1-a)/.5*(y-nu0)),h[i,]-mu)
    oracle <- oracle + mass*tcrossprod(phi)
  }
  expect_equal(unname(fit$variance$trials$s1$working$Sigma), oracle, tolerance=1e-10)
})

test_that("the intercept-only case reduces to the reported difference variance", {
  d <- data.frame(trial="s", subgroup="overall",stratum="TRUE",
    event1=30,n1=100,event2=20,n2=100)
  fit <- CIMA(d, mu=list(s=numeric()), Q=data.frame(x=1:20),
    formula_cate=~1,formula_tilting=~1)
  expected <- .3*.7/100+.2*.8/100
  expect_equal(unname(coef(fit)),.1,tolerance=1e-12)
  expect_equal(as.numeric(vcov(fit)),expected,tolerance=1e-12)
  expect_equal(fit$ATE$var,expected,tolerance=1e-12)
})

test_that("the summary derivative has the positive H inverse sign", {
  fit <- fit_small(saturated=FALSE); tilt <- fit$tilting$s1
  args <- small_example(saturated=FALSE)
  ii <- indicator_matrix(args$Q,prepare_agd(args$data_agd))
  g <- as.vector(model.matrix(~x+z,args$Q)%*%coef(fit))
  ctl <- list(rank_tol=1e-9,condition_limit=1e8,tilt_tol=1e-11,max_iter=200L)
  moment <- function(mu) {
    tt <- fit_tilt(tilt$h_plus,mu,"s1",ctl)
    colSums(ii*(tt$weights*g))/colSums(ii*tt$weights)
  }
  step <- 1e-5
  derivative <- vapply(2:3,function(j) {
    delta <- rep(0,3);delta[j]<-step
    (moment(tilt$mu_plus+delta)-moment(tilt$mu_plus-delta))/(2*step)
  },numeric(5))
  expect_equal(unname(fit$variance$J_m_mu$s1[,-1]),unname(derivative),tolerance=1e-7)
})

test_that("binary working regression enforces active pointwise risk constraints", {
  args <- small_example(); d <- prepare_agd(args$data_agd)
  Q <- args$Q; h <- model.matrix(~x+z,Q)
  ii <- indicator_matrix(Q,d); alpha <- sweep(ii,2,colMeans(ii),"/")
  d$nu0 <- c(.5,.05,.95,.5,.5)
  d$se_nu0 <- rep(.02,5)
  ctl <- list(rank_tol=1e-9)
  fit <- fit_working_control(h,alpha,rep(1,nrow(Q)),rep(.2,nrow(Q)),d,"binary",ctl)
  expect_true(fit$constrained)
  expect_gte(min(fit$gamma0),0)
  expect_lte(max(fit$gamma1),1)
  expect_equal(fit$gamma1-fit$gamma0,rep(.2,nrow(Q)),tolerance=1e-10)
})
