# Numerical and input-validation utilities.
`%||%` <- function(x, y) if (is.null(x)) y else x

symmetrize <- function(x) (x + t(x)) / 2

center_columns <- function(x) sweep(x, 2L, colMeans(x), "-")

empirical_cov <- function(x) crossprod(center_columns(as.matrix(x))) / nrow(x)

check_scalar <- function(x, name, lower = -Inf, integer = FALSE) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < lower ||
      (integer && x != floor(x))) stop(name, " must be a finite ",
        if (integer) "integer " else "number ", ">= ", lower, ".", call. = FALSE)
  invisible(x)
}

with_seed <- function(seed, code) {
  if (is.null(seed)) return(force(code))
  check_scalar(seed, "seed", 0, TRUE)
  had_seed <- exists(".Random.seed", .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", .GlobalEnv)
  on.exit(if (had_seed) assign(".Random.seed", old_seed, .GlobalEnv) else
    if (exists(".Random.seed", .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
  set.seed(seed)
  force(code)
}

as_rhs_formula <- function(x, name) {
  f <- if (inherits(x, "formula")) x else stats::as.formula(x)
  if (length(f) != 2L) stop(name, " must be a one-sided formula.", call. = FALSE)
  f
}

basis_matrix <- function(formula, data, label) {
  mf <- stats::model.frame(formula, data, na.action = stats::na.pass)
  b <- stats::model.matrix(formula, mf)
  if (nrow(b) != nrow(data) || any(!is.finite(b)) || !ncol(b))
    stop(label, " must produce a nonempty finite design with no missing rows.", call. = FALSE)
  b
}

# Subgroup expressions deliberately allow only comparisons and Boolean operations.
# Labels that are not expressions are matched literally (including factor levels).
stratum_indicator <- function(data, subgroup, stratum) {
  if (subgroup == "overall") return(rep(1, nrow(data)))
  if (!subgroup %in% names(data)) stop("Subgroup variable missing from Q: ", subgroup, call. = FALSE)
  x <- data[[subgroup]]
  expr <- tryCatch(parse(text = stratum), error = function(e) NULL)
  allowed <- c("==", "!=", "<", "<=", ">", ">=", "&", "|", "!", "(", "%in%", "c", "-", "+")
  valid <- function(e) {
    if (is.atomic(e)) return(TRUE)
    if (is.name(e)) return(as.character(e) %in% c("x", subgroup))
    is.call(e) && as.character(e[[1L]]) %in% allowed &&
      all(vapply(as.list(e)[-1L], valid, logical(1)))
  }
  if (!is.null(expr) && length(expr) == 1L && is.call(expr[[1L]])) {
    if (!valid(expr[[1L]])) stop("Unsupported subgroup expression: ", stratum, call. = FALSE)
    env <- list2env(stats::setNames(list(x, x), c("x", subgroup)), parent = baseenv())
    ans <- eval(expr[[1L]], envir = env)
  } else if (identical(stratum, "TRUE")) {
    ans <- rep(TRUE, nrow(data))
  } else {
    ans <- as.character(x) == stratum
  }
  if (!is.logical(ans) || length(ans) != nrow(data) || anyNA(ans))
    stop("Stratum must define one nonmissing logical value per observation: ", stratum, call. = FALSE)
  as.numeric(ans)
}

indicator_matrix <- function(data, reports) {
  cache <- new.env(parent = emptyenv())
  out <- vapply(seq_len(nrow(reports)), function(j) {
    subgroup <- reports$subgroup[j]; stratum <- reports$stratum[j]
    key <- paste0(nchar(subgroup), ":", subgroup, stratum)
    if (!exists(key, envir = cache, inherits = FALSE))
      assign(key, stratum_indicator(data, subgroup, stratum), envir = cache)
    get(key, envir = cache, inherits = FALSE)
  }, numeric(nrow(data)))
  matrix(out, nrow(data), nrow(reports), dimnames = list(NULL, reports$row_id))
}

validate_covariance <- function(x, size, name, positive = FALSE, tol = 1e-8) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !identical(dim(x), c(as.integer(size), as.integer(size))) ||
      any(!is.finite(x))) stop(name, " must be a finite ", size, " by ", size, " matrix.", call. = FALSE)
  scale <- max(abs(x), 1)
  if (max(abs(x - t(x))) > tol * scale) stop(name, " must be symmetric.", call. = FALSE)
  x <- symmetrize(x)
  d <- sqrt(pmax(diag(x), .Machine$double.eps * scale))
  ev <- eigen(x / outer(d, d), symmetric = TRUE, only.values = TRUE)$values
  if (min(ev) < -tol || (positive && min(ev) <= tol))
    stop(name, if (positive) " must be positive definite." else " must be positive semidefinite.", call. = FALSE)
  x
}

named_trial_item <- function(x, s, name) {
  if (!is.list(x) || is.null(names(x)) || !s %in% names(x))
    stop(name, " must be a list named by trial; missing '", s, "'.", call. = FALSE)
  x[[s]]
}

minimum_norm <- function(a, y, tol = 1e-10) {
  z <- svd(a)
  keep <- z$d > max(z$d, 0) * tol
  if (!any(keep)) return(rep(0, ncol(a)))
  as.vector(z$v[, keep, drop = FALSE] %*%
              (crossprod(z$u[, keep, drop = FALSE], y) / z$d[keep]))
}
