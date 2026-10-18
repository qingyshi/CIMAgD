small_example <- function(outcome = "binary", saturated = TRUE) {
  Q <- expand.grid(x = 0:1, z = 0:1)[rep(1:4, each = 50), ]
  grid <- expand.grid(x = 0:1, z = 0:1)
  prob <- c(.1, .2, .3, .4)
  rows <- data.frame(trial = "s1", subgroup = c("overall", "x", "x", "z", "z"),
    stratum = c("TRUE", "x == 0", "x == 1", "x == 0", "x == 1"))
  indicators <- cbind(1, grid$x == 0, grid$x == 1, grid$z == 0, grid$z == 1)
  pr <- colSums(indicators * prob)
  conditional <- sweep(indicators * prob, 2, pr, "/")
  theta <- if (outcome == "binary") c(.1, .05, -.02) else c(1, .5, -.2)
  g <- as.vector(cbind(1, grid$x, grid$z) %*% theta)
  c0 <- if (outcome == "binary") .3 + .03 * grid$x + .02 * grid$z else 2 + .3 * grid$x + .2 * grid$z
  nu0 <- colSums(conditional * c0); nu1 <- colSums(conditional * (c0 + g))
  rows$n1 <- rows$n2 <- as.integer(pr * 10000)
  if (outcome == "binary") {
    rows$event1 <- round(nu1 * rows$n1); rows$event2 <- round(nu0 * rows$n2)
  } else {
    rows$mean1 <- nu1; rows$mean2 <- nu0
    rows$se1 <- sqrt((2 + colSums(conditional * sweep(matrix(c0 + g, 4, 5), 2, nu1, "-")^2)) / rows$n1)
    rows$se2 <- sqrt((1 + colSums(conditional * sweep(matrix(c0, 4, 5), 2, nu0, "-")^2)) / rows$n2)
  }
  mu <- colSums(cbind(grid$x, grid$z, grid$x * grid$z) * prob)
  if (!saturated) mu <- mu[1:2]
  list(data_agd = rows, mu = list(s1 = mu), Q = Q, target = Q,
       formula_cate = ~ x + z, formula_tilting = if (saturated) ~ x * z else ~ x + z)
}

fit_small <- function(outcome = "binary", saturated = TRUE, ...) {
  do.call(CIMA, c(small_example(outcome, saturated), list(...)))
}

small_copula <- function(n = 600L) list(n = n,
  marginals = data.frame(variable = c("x", "z"), distribution = "bernoulli", mean = c(.45, .55)),
  correlation = matrix(c(1, .2, .2, 1), 2))
