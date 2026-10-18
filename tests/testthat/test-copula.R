test_that("copula input has one source, preserves RNG and matches margins", {
  spec <- small_copula(4000)
  set.seed(4); before <- .Random.seed
  q <- XCopula(spec, seed = 22)
  expect_identical(.Random.seed, before)
  expect_equal(q, XCopula(spec, seed = 22))
  expect_equal(nrow(q), 4000L)
  expect_equal(colMeans(q), c(x=.45,z=.55), tolerance=.025)
  expect_false(identical(q, XCopula(spec, seed=23)))
  expect_error(XCopula(list(trials=list(s1=spec))), "one source")
  spec$correlation <- matrix(c(1,2,2,1),2)
  expect_error(XCopula(spec), "positive definite")
})

test_that("plugin uses original sample size and bootstrap refits both margins and dependence", {
  a <- small_example(); a$Q <- a$target <- NULL; a$summaryX <- small_copula(800)
  plugin <- do.call(CIMA, c(a, list(seed = 123)))
  expect_equal(plugin$copula$integration_n, 800L)
  expect_true(plugin$target$shared_q)
  boot <- do.call(CIMA, c(a, list(q_method="bootstrap", bootstrap_reps=12,
    bootstrap_mc_size=1800, seed=123)))
  expect_equal(boot$copula$bootstrap$successful, 12L)
  expect_equal(boot$copula$integration_n, 1800L)
  expect_gt(stats::sd(vapply(boot$copula$bootstrap$refitted_summaries,
    function(s) s$marginals$mean[1], numeric(1))), 0)
  expect_gt(stats::sd(vapply(boot$copula$bootstrap$refitted_summaries,
    function(s) s$correlation[1,2], numeric(1))), 0)
  expect_equal(vcov(boot), boot$variance$trial_covariance +
    boot$copula$bootstrap$theta_covariance, ignore_attr=TRUE)
  expect_equal(vcov(boot,"target"), Reduce(`+`, boot$target$components), ignore_attr=TRUE)
  expect_gt(boot$target$components$copula_bootstrap[1,1], 0)
})

test_that("mixed copula refit recovers latent rather than observed binary correlation", {
  spec <- list(n=15000, marginals=data.frame(variable=c("a","b","c"),
    distribution=c("normal","bernoulli","bernoulli"), mean=c(2,.3,.6), sd=c(1.5,NA,NA)),
    correlation=matrix(c(1,.4,.2,.4,1,.5,.2,.5,1),3))
  q <- XCopula(spec, seed=84); fitted <- refit_copula(q, validate_summary(spec))
  expect_equal(fitted$marginals$mean, spec$marginals$mean, tolerance=.025)
  expect_equal(unname(fitted$correlation), spec$correlation, tolerance=.04)
  expect_gt(abs(fitted$correlation[2,3] - cor(q$b,q$c)), .1)
})

test_that("bootstrap retains independent target and trial terms without adding plugin Q twice", {
  a <- small_example(saturated=FALSE)
  a$Q <- NULL; a$summaryX <- small_copula(500)
  fit <- do.call(CIMA,c(a,list(q_method="bootstrap",bootstrap_reps=8,
    bootstrap_mc_size=1500,seed=42)))
  expect_false(fit$target$shared_q)
  expect_gt(sum(diag(fit$copula$bootstrap$theta_covariance)),1e-10)
  expect_gt(sum(diag(fit$target$components$independent_target)),0)
  expect_equal(vcov(fit,"target"),Reduce(`+`,fit$target$components),ignore_attr=TRUE)
})

test_that("friendly marginal lists accept native parameters or means and SDs", {
  m<-list(n=list(distribution="normal",mean=2,sd=.5),
    g=list(distribution="gamma",mean=2,sd=1),
    b=list(distribution="beta",mean=.4,sd=.2),
    l=list(distribution="lognormal",mean=2,sd=1),
    e=list(distribution="exponential",mean=3),
    d=list(distribution="bernoulli",prob=.3))
  s<-copula_summary(12000,m)
  q<-XCopula(s,seed=415)
  expect_equal(names(q),names(m))
  expect_lt(max(abs(colMeans(q)-s$marginals$mean)),.06)
  expect_lt(max(abs(vapply(q,sd,numeric(1))-s$marginals$sd)),.07)
  expect_true(all(q$g>0 & q$l>0 & q$e>0 & q$b>0 & q$b<1))
  native<-copula_summary(12000,list(n=m$n,g=list(distribution="gamma",shape=4,rate=2),
    b=list(distribution="beta",shape1=2,shape2=3),
    l=list(distribution="lognormal",meanlog=log(2)-log(1.25)/2,sdlog=sqrt(log(1.25))),
    e=list(distribution="exponential",rate=1/3),d=m$d))
  expect_equal(XCopula(native,seed=415),q,tolerance=1e-12)
  expect_equal(XCopula(list(n=s$n,marginals=s$marginals,correlation=s$correlation),seed=415),q)
  scale<-copula_summary(12,list(g=list(distribution="gamma",shape=4,scale=.5)))
  expect_equal(scale$marginals$mean,2)
  expect_error(copula_summary(100,list(b=list(distribution="beta",mean=.4,sd=.6))),"Beta moments")
  expect_error(copula_summary(100,list(g=list(distribution="gamma",mean=-1,sd=1))),"positive mean")
  expect_error(copula_summary(100,list(g=list(distribution="gamma",shape=0,rate=1))),"positive shape")
})

test_that("non-normal bootstrap refits margins and latent dependence", {
  R<-matrix(.3,4,4);diag(R)<-1
  s<-copula_summary(10000,list(g=list(distribution="gamma",mean=2,sd=1),
    b=list(distribution="beta",mean=.4,sd=.2),
    l=list(distribution="lognormal",meanlog=.5,sdlog=.4),
    d=list(distribution="bernoulli",prob=.4)),R)
  q<-XCopula(s,seed=14);refit<-refit_copula(q,s)
  expect_lt(max(abs(refit$marginals$mean-s$marginals$mean)),.04)
  expect_lt(max(abs(refit$correlation-s$correlation)),.04)
  expect_gt(max(abs(refit$marginals$shape-s$marginals$shape),na.rm=TRUE),0)
  expect_equal(refit$marginals$parameterization[3],"log-moments")
  expect_equal(refit$marginals$meanlog[3],mean(log(q$l)))
})

test_that("new continuous margins run through plugin and bootstrap CIMA", {
  a<-small_example("continuous",saturated=FALSE)
  a$data_agd<-a$data_agd[1,,drop=FALSE]
  a$Q<-a$target<-NULL
  a$formula_cate<-~1;a$formula_tilting<-~x
  a$mu<-list(s1=c(x=2))
  a$summaryX<-copula_summary(400,list(x=list(distribution="gamma",mean=2,sd=1)))
  plugin<-do.call(CIMA,c(a,list(seed=34)))
  boot<-do.call(CIMA,c(a,list(q_method="bootstrap",bootstrap_reps=8,bootstrap_mc_size=800,seed=34)))
  expect_equal(coef(plugin),coef(boot),tolerance=1e-10)
  expect_equal(boot$copula$bootstrap$successful,8L)
  means<-vapply(boot$copula$bootstrap$refitted_summaries,function(s)s$marginals$mean,numeric(1))
  expect_gt(sd(means),0)
})
