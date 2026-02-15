
utils::globalVariables(
  )

# internal helper: quantile functions
qnormal <- function(p, mean, sd, min, max) stats::qnorm(p, mean, sd)
qbinomial <- function(p, mean, n, min, max) stats::qbinom(p, n, mean)
qtruncnormal <- function(p, mean, sd, min, max) {
  mean + sd * stats::qnorm(
    stats::pnorm(min, mean, sd) + p * (
      stats::pnorm(max, mean, sd) - stats::pnorm(min, mean, sd)
    )
  )
}


# internal helper: indicator functions for subgroup and stratum
IFun <- function(x, s, dataAgD) {
  dat_s <- dataAgD[dataAgD$trial == s, , drop = FALSE]
  if (nrow(dat_s) == 0) {
    return(numeric(0))
  }

  out <- rep(1, nrow(dat_s))
  nms <- names(x)

  for (j in seq_len(nrow(dat_s))) {
    subgroup_j <- as.character(dat_s$subgroup[j])
    stratum_j <- as.character(dat_s$stratum[j])

    if (is.na(subgroup_j) || subgroup_j == "overall" || is.na(stratum_j) ||
      stratum_j == "TRUE") {
      out[j] <- 1
      next
    }

    if (!(subgroup_j %in% nms)) {
      out[j] <- 0
      next
    }

    val <- x[[subgroup_j]]
    env <- new.env(parent = baseenv())
    env$x <- val
    env[[subgroup_j]] <- val

    test <- try(eval(parse(text = stratum_j), envir = env), silent = TRUE)
    if (inherits(test, "try-error") || length(test) != 1L || is.na(test)) {
      out[j] <- as.numeric(as.character(val) == stratum_j)
    } else {
      out[j] <- as.numeric(isTRUE(test))
    }
  }

  out
}
IMatFun <- function(Q, dataAgD) {
  trials <- unique(as.character(dataAgD$trial))
  Q_df <- as.data.frame(Q)

  out <- lapply(trials, function(s) {
    mat <- t(vapply(
      seq_len(nrow(Q_df)),
      function(i) IFun(Q_df[i, , drop = FALSE], s, dataAgD),
      FUN.VALUE = numeric(sum(dataAgD$trial == s))
    ))
    if (is.null(dim(mat))) {
      mat <- matrix(mat, ncol = sum(dataAgD$trial == s))
    }
    colnames(mat) <- rownames(dataAgD[dataAgD$trial == s, , drop = FALSE])
    mat
  })

  names(out) <- trials
  out
}
