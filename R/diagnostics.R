#' Diagnose rank and conditioning of a Jacobian
#'
#' @param x Numeric matrix.
#' @param tol Relative singular-value tolerance after column scaling.
#' @param condition_limit Condition number above which a full-rank matrix is weak.
#' @return Rank, singular values, raw/scaled condition numbers, column scales,
#'   right null directions and status. Columns are scaled to unit Euclidean norm.
#' @export
diagnose_jacobian <- function(x, tol = 1e-9, condition_limit = 1e8) {
  x <- as.matrix(x)
  if (!is.numeric(x) || !length(x) || any(!is.finite(x))) stop("x must be a finite numeric matrix.", call. = FALSE)
  check_scalar(tol, "tol", .Machine$double.eps)
  check_scalar(condition_limit, "condition_limit", 1)
  scales <- sqrt(colSums(x^2)); scales[scales == 0] <- 1
  z <- sweep(x, 2, scales, "/")
  ss <- svd(z, nu = 0, nv = ncol(z))
  values <- c(ss$d, rep(0, max(0, ncol(x) - length(ss$d))))
  rank <- sum(values > max(values) * tol)
  raw <- svd(x, nu = 0, nv = 0)$d
  condition <- if (rank < ncol(x)) Inf else max(values) / min(values)
  null <- if (rank < ncol(x)) sweep(ss$v[, seq.int(rank + 1L, ncol(x)), drop = FALSE], 1, scales, "/") else matrix(numeric(), ncol(x), 0)
  rownames(null) <- colnames(x)
  list(dim = dim(x), rank = rank, full_rank = rank == ncol(x),
       singular_values = values, raw_singular_values = raw,
       condition = condition, raw_condition = if (length(raw) < ncol(x) || min(raw) == 0) Inf else max(raw) / min(raw),
       column_scales = scales, null_directions = null,
       status = if (rank < ncol(x)) "singular" else if (condition > condition_limit) "ill_conditioned" else "ok")
}

require_identified <- function(x, label, tol, condition_limit, warn = TRUE) {
  dg <- diagnose_jacobian(x, tol, condition_limit)
  if (!dg$full_rank) {
    cond <- structure(list(message = paste0(label, " is singular: rank ", dg$rank, " of ", ncol(x),
      ". Remove redundant terms or add identifying summaries/support."), call = NULL,
      matrix_name = label, diagnostics = dg), class = c("cima_singularity", "error", "condition"))
    stop(cond)
  }
  if (warn && dg$status == "ill_conditioned") warning(label, " is ill-conditioned (scaled condition ", signif(dg$condition, 4), ").", call. = FALSE)
  dg
}

#' Check the sufficient tilting span condition under Q
#'
#' @param h_plus Tilting basis matrix, including its intercept, evaluated on Q.
#' @param basis CATE basis on the same Q observations.
#' @param indicators Subgroup indicator matrix on Q.
#' @param tol Relative residual-norm tolerance.
#' @return Projection coefficients Pi_s, residual Gram matrix R_s, residuals,
#'   column diagnostics and an empirical sufficient flag. Failure is not evidence
#'   that A6 fails, and this does not check population support or dominance.
#' @export
diagnose_tilting <- function(h_plus, basis, indicators, tol = 1e-8) {
  h <- as.matrix(h_plus); b <- as.matrix(basis); ii <- as.matrix(indicators)
  if (nrow(h) != nrow(b) || nrow(h) != nrow(ii) || any(!is.finite(c(h, b, ii))))
    stop("Diagnostic matrices must have equal rows and finite entries.", call. = FALSE)
  colnames(b) <- colnames(b) %||% paste0("b", seq_len(ncol(b)))
  colnames(ii) <- colnames(ii) %||% paste0("I", seq_len(ncol(ii)))
  agg <- do.call(cbind, lapply(seq_len(ncol(ii)), function(j) {
    v <- cbind(ii[, j], b * ii[, j])
    colnames(v) <- c(paste0(colnames(ii)[j], ":I"), paste0(colnames(ii)[j], ":I*", colnames(b)))
    v
  }))
  scales <- sqrt(colMeans(h^2)); scales[scales == 0] <- 1
  hs <- sweep(h, 2, scales, "/")
  z <- svd(hs); keep <- z$d > max(z$d) * tol
  pi <- matrix(0, ncol(h), ncol(agg))
  if (any(keep)) pi <- sweep(z$v[, keep, drop = FALSE] %*%
    (crossprod(z$u[, keep, drop = FALSE], agg) / z$d[keep]), 1, scales, "/")
  dimnames(pi) <- list(colnames(h), colnames(agg))
  residual <- agg - h %*% pi
  R <- crossprod(residual) / nrow(h)
  rms <- sqrt(pmax(diag(R), 0)); norm <- sqrt(colMeans(agg^2))
  relative <- ifelse(norm > 0, rms / pmax(norm, .Machine$double.eps), rms)
  columns <- data.frame(term = colnames(agg), residual_norm = rms, relative_residual = relative,
                        in_span = relative <= tol, row.names = NULL)
  list(Pi_s = pi, R_s = R, residuals = residual, columns = columns,
       sufficient = all(columns$in_span), unresolved = columns$term[!columns$in_span],
       interpretation = "Empirical sufficient span check under Q, not a test of A6.")
}
