test_that("two-step uses full pilot moment covariance and recomputes final variance", {
  a <- small_example("continuous",saturated=FALSE)
  a$data_agd$mean1[1] <- a$data_agd$mean1[1]+.02
  pilot <- do.call(CIMA,a)
  fit <- do.call(CIMA,c(a,list(gmm_method="two-step")))
  C <- pilot$variance$moment_covariance
  expect_equal(fit$gmm$optimal$covariance,C,tolerance=1e-12)
  D <- fit$CATE$J_m_theta; W <- fit$moments$W
  expect_equal(unname(W %*% C),diag(nrow(C)),tolerance=1e-6)
  expected <- solve(crossprod(D,W%*%D),crossprod(D,W%*%fit$data$TE))
  expect_equal(unname(coef(fit)),as.vector(expected),tolerance=1e-8)
  expect_lt(max(abs(crossprod(D,W%*%fit$moments$residual))),1e-6)
  expect_equal(fit$gmm$initial_theta,coef(pilot))
  expect_equal(fit$gmm$steps,2L)
  expect_equal(fit$variance$moment_covariance,
    fit$variance$Gamma_q/nrow(a$Q)+fit$variance$trial_moment_covariance,tolerance=1e-12)
  J <- -solve(crossprod(D,W%*%D),crossprod(D,W))
  expect_equal(vcov(fit),J%*%fit$variance$moment_covariance%*%t(J),ignore_attr=TRUE,tolerance=1e-8)
})

test_that("two-step handles redundant moments without arbitrary covariance ridges", {
  pilot <- fit_small()
  sigma <- matrix(0,9,9)
  sigma[1:5,1:5] <- tcrossprod(pilot$CATE$J_m_theta)
  fit <- fit_small(gmm_method="two-step",Sigma=list(s1=sigma))
  expect_equal(unname(coef(fit)),c(.1,.05,-.02),tolerance=1e-9)
  expect_lt(fit$gmm$optimal$rank,nrow(fit$data))
  expect_equal(fit$gmm$optimal$ridge,0)
  expect_false(fit$regularization$used)
  expect_true(all(is.finite(vcov(fit))))
  point <- fit_small(gmm_method="two-step",variance=FALSE)
  expect_equal(coef(point),coef(fit),tolerance=1e-10)
  expect_true(all(is.na(vcov(point))))
})

test_that("informative noiseless moments are not silently discarded", {
  ctl<-list(rank_tol=1e-9,condition_limit=1e8,singularity_method="none",ridge_strength=1e-8)
  expect_error(optimal_weight(diag(c(1,0)),diag(2),ctl),"informative zero-variance")
  ctl$singularity_method<-"ridge"
  w<-optimal_weight(diag(c(1,0)),diag(2),ctl)
  expect_equal(w$ridge,1e-8)
  expect_gt(min(eigen(w$W,symmetric=TRUE)$values),0)
})

test_that("two-step copula bootstrap uses refit moments rather than integration noise", {
  a<-small_example(saturated=FALSE);a$Q<-a$target<-NULL;a$summaryX<-small_copula(700)
  fit<-do.call(CIMA,c(a,list(q_method="bootstrap",gmm_method="two-step",
    bootstrap_reps=12,bootstrap_mc_size=1800,seed=621)))
  expect_equal(fit$gmm$steps,2L)
  expect_match(fit$gmm$optimal$source,"refitted copula")
  expect_equal(fit$copula$bootstrap$successful,12L)
  expect_false(fit$copula$bootstrap$weight_refitted)
  expect_equal(vcov(fit),fit$variance$trial_covariance+fit$copula$bootstrap$theta_covariance,ignore_attr=TRUE)
  expect_equal(fit$variance$moment_covariance,fit$variance$trial_moment_covariance+
    fit$copula$bootstrap$moment_covariance,ignore_attr=TRUE)
})
