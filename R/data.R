#' Prepare contrast-based aggregate reports
#'
#' @param data_agd Data frame with trial, subgroup, stratum and both arms' reports.
#'   Binary reports require event1, n1, event2, n2. Continuous reports require
#'   mean1, se1, n1, mean2, se2, n2; se1/se2 are standard errors of means.
#' @param outcome One of auto, binary, continuous.
#' @return A validated data frame with missing TE, seTE, n, comparison and Time
#'   filled. Arm 1 is treatment; arm 2 is control. No continuity correction is used.
#' @export
prepare_agd <- function(data_agd, outcome = c("auto", "binary", "continuous")) {
  outcome <- match.arg(outcome)
  d <- as.data.frame(data_agd, stringsAsFactors = FALSE)
  required <- c("trial", "subgroup", "stratum", "n1", "n2")
  missing <- setdiff(required, names(d))
  if (length(missing)) stop("Missing aggregate columns: ", paste(missing, collapse = ", "), call. = FALSE)
  if (!nrow(d)) stop("data_agd has no rows.", call. = FALSE)
  for (nm in c("trial", "subgroup", "stratum")) {
    d[[nm]] <- as.character(d[[nm]])
    if (anyNA(d[[nm]]) || any(!nzchar(d[[nm]]))) stop("Missing values in ", nm, ".", call. = FALSE)
  }
  has_binary <- all(c("event1", "event2") %in% names(d))
  has_continuous <- all(c("mean1", "se1", "mean2", "se2") %in% names(d))
  if (outcome == "auto") {
    if (has_binary && has_continuous) stop("Both reporting formats present; specify outcome.", call. = FALSE)
    outcome <- if (has_binary) "binary" else "continuous"
  }
  cols <- c("n1", "n2", if (outcome == "binary") c("event1", "event2") else c("mean1", "se1", "mean2", "se2"))
  if (!all(cols %in% names(d))) stop("Incomplete ", outcome, " arm reports; require ", paste(cols, collapse = ", "), ".", call. = FALSE)
  for (nm in cols) if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]])))
    stop(nm, " must be numeric and finite.", call. = FALSE)
  if (any(d$n1 <= 0 | d$n2 <= 0 | d$n1 != floor(d$n1) | d$n2 != floor(d$n2)))
    stop("Arm sample sizes must be positive integers.", call. = FALSE)
  if (outcome == "binary") {
    for (a in 1:2) if (any(d[[paste0("event", a)]] < 0 |
      d[[paste0("event", a)]] > d[[paste0("n", a)]] |
      d[[paste0("event", a)]] != floor(d[[paste0("event", a)]])))
      stop("Event counts must be integers between zero and arm sample size.", call. = FALSE)
    nu1 <- d$event1 / d$n1
    nu0 <- d$event2 / d$n2
    se1 <- sqrt(nu1 * (1 - nu1) / d$n1)
    se0 <- sqrt(nu0 * (1 - nu0) / d$n2)
  } else {
    if (any(d$se1 < 0 | d$se2 < 0)) stop("Arm standard errors cannot be negative.", call. = FALSE)
    nu1 <- d$mean1; nu0 <- d$mean2; se1 <- d$se1; se0 <- d$se2
  }
  fill <- function(name, value) {
    if (!name %in% names(d)) d[[name]] <<- value else {
      miss <- is.na(d[[name]])
      d[[name]][miss] <<- rep_len(value, nrow(d))[miss]
    }
  }
  if ("comparison" %in% names(d)) d$comparison <- as.character(d$comparison)
  fill("comparison", "A vs B")
  d$comparison <- as.character(d$comparison)
  if (length(unique(d$comparison)) != 1L || any(!nzchar(d$comparison)))
    stop("CIMA estimates one consistently oriented comparison at a time.", call. = FALSE)
  fill("TE", nu1 - nu0); fill("seTE", sqrt(se1^2 + se0^2))
  fill("n", d$n1 + d$n2); fill("Time", 1)
  for (nm in c("TE", "seTE", "n", "Time"))
    if (!is.numeric(d[[nm]]) || any(!is.finite(d[[nm]]))) stop(nm, " must be finite and numeric.", call. = FALSE)
  if (any(d$seTE < 0)) stop("seTE cannot be negative.", call. = FALSE)
  if (any(d$n != d$n1 + d$n2)) stop("n must equal n1 + n2 for each report.", call. = FALSE)
  if (outcome == "binary" && any(abs(d$TE) > 1)) stop("Binary TE must be a risk difference in [-1, 1].", call. = FALSE)
  if (anyDuplicated(d[c("trial", "comparison", "subgroup", "stratum", "Time")]))
    stop("Duplicate trial/subgroup/stratum/Time report.", call. = FALSE)
  # These internal columns keep the reported format separate from estimator notation.
  d$nu1 <- nu1; d$nu0 <- nu0; d$se_nu1 <- se1; d$se_nu0 <- se0
  d$row_id <- paste(d$trial, d$subgroup, d$stratum, d$Time, sep = ":")
  rownames(d) <- NULL
  attr(d, "outcome") <- outcome
  d
}

resolve_trial_info <- function(d, indicators, trial_info = NULL) {
  trials <- unique(d$trial)
  if (!is.null(trial_info)) {
    trial_info <- as.data.frame(trial_info)
    if (!all(c("trial", "n", "allocation") %in% names(trial_info)) ||
        anyDuplicated(trial_info$trial)) stop("trial_info requires unique trial, n, allocation columns.", call. = FALSE)
  }
  out <- lapply(trials, function(s) {
    ds <- d[d$trial == s, , drop = FALSE]
    if (!is.null(trial_info) && s %in% trial_info$trial) {
      z <- trial_info[match(s, trial_info$trial), ]
      n_s <- z$n; allocation <- z$allocation
    } else {
      overall <- which(ds$subgroup == "overall")
      if (length(overall)) {
        if (length(unique(ds$n[overall])) != 1L || length(unique(ds$n1[overall])) != 1L)
          stop("Different overall samples within trial ", s, "; provide trial_info and a joint Sigma.", call. = FALSE)
        n_s <- ds$n[overall[1L]]; allocation <- ds$n1[overall[1L]] / n_s
      } else {
        candidates <- list()
        for (sg in unique(ds$subgroup)) for (tt in unique(ds$Time)) {
          ii <- which(d$trial == s & d$subgroup == sg & d$Time == tt)
          if (length(ii) > 1L && all(rowSums(indicators[, ii, drop = FALSE]) == 1))
            candidates[[length(candidates) + 1L]] <- c(sum(d$n[ii]), sum(d$n1[ii]))
        }
        if (!length(candidates)) stop("Cannot identify total n/allocation for trial ", s,
          "; supply trial_info or an overall/exhaustive partition report.", call. = FALSE)
        cc <- do.call(rbind, candidates)
        if (length(unique(cc[, 1])) != 1L || length(unique(cc[, 2])) != 1L)
          stop("Inconsistent partition totals for trial ", s, "; supply trial_info.", call. = FALSE)
        n_s <- cc[1, 1]; allocation <- cc[1, 2] / n_s
      }
    }
    check_scalar(n_s, paste0("trial_info n for ", s), 2, TRUE)
    check_scalar(allocation, paste0("allocation for ", s), 0)
    if (allocation <= 0 || allocation >= 1 || any(ds$n > n_s))
      stop("Invalid allocation or subgroup n larger than trial n for ", s, ".", call. = FALSE)
    data.frame(trial = s, n = n_s, allocation = allocation)
  })
  do.call(rbind, out)
}
