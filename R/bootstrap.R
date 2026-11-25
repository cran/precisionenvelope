envelope_bootstrap <- function(fit, R = 999L, alpha = 0.05,
                               failure = c("abort", "discard"),
                               min_success_fraction = 0.75,
                               min_success_count = 40L) {
  if (!inherits(fit, "envelope_fit")) {
    .stop("`fit` must be an object returned by envelope_fit().")
  }
  if (is.null(fit$requirement)) {
    .stop("the fit has no requirement; refit with `target` or `variance_limit`.")
  }
  failure <- match.arg(failure)
  if (!.is_number(R) || R < 1 || R != round(R)) {
    .stop("`R` must be a positive whole number.")
  }
  if (!.is_number(alpha) || alpha <= 0 || alpha >= 1) {
    .stop("`alpha` must lie strictly between 0 and 1.")
  }
  R <- as.integer(R)
  required <- if (failure == "abort") {
    R
  } else {
    if (!.is_number(min_success_fraction) || min_success_fraction <= 0 ||
        min_success_fraction > 1) {
      .stop("`min_success_fraction` must lie in (0, 1].")
    }
    max(as.integer(min_success_count), ceiling(min_success_fraction * R))
  }
  if (required > R) {
    .stop("`R` is smaller than the required number of successful refits.")
  }

  theta <- fit$theta
  x <- fit$data$x
  r <- fit$data$r
  y <- fit$data$y
  n <- length(r)
  draws <- matrix(NA_real_, R, length(theta))
  ok <- logical(R)
  for (i in seq_len(R)) {
    idx <- sample.int(n, n, replace = TRUE)
    value <- tryCatch(
      .envelope_point(x[idx, , drop = FALSE], r[idx], y[idx], fit$spec, fit$learners)$theta,
      error = function(err) NULL
    )
    if (!is.null(value) && length(value) == length(theta) && all(is.finite(value))) {
      draws[i, ] <- value
      ok[i] <- TRUE
    }
  }
  successes <- sum(ok)
  if (successes < required) {
    .stop("only ", successes, " of ", R, " bootstrap refits succeeded; ",
          required, " were required under failure = \"", failure, "\".")
  }
  draws <- draws[ok, , drop = FALSE]
  scale <- apply(draws, 2L, stats::sd)
  if (any(!is.finite(scale)) || any(scale <= 1e-12)) {
    .stop("the bootstrap produced a coordinate with no variation.")
  }
  maxima <- apply(abs(sweep(sweep(draws, 2L, theta), 2L, scale, "/")), 1L, max)
  rank <- ceiling((successes + 1) * (1 - alpha))
  if (rank > successes) {
    .stop("at least ", ceiling((1 - alpha) / alpha),
          " successful refits are needed for alpha = ", alpha, ".")
  }
  critical <- sort(maxima)[rank]
  lower_band <- theta - critical * scale
  upper_band <- theta + critical * scale

  k <- length(fit$spec$p)
  g_lo <- lower_band[seq_len(k)]
  g_hi <- upper_band[seq_len(k)]
  p_lo <- lower_band[k + seq_len(k)]
  p_hi <- upper_band[k + seq_len(k)]

  feasibility <- "not applicable"
  allow_phase <- TRUE
  if (fit$requirement$type == "variance") {
    f_lo <- lower_band[length(theta)]
    f_hi <- upper_band[length(theta)]
    if (f_lo > 0) {
      feasibility <- "feasible with complete data"
    } else if (f_hi < 0) {
      feasibility <- "infeasible even with complete data"
      allow_phase <- FALSE
    } else {
      feasibility <- "indeterminate"
      allow_phase <- FALSE
    }
  }
  phase <- rep("indeterminate", k)
  thresholds <- NULL
  if (allow_phase) {
    phase[g_lo > 0] <- "guaranteed"
    phase[p_hi < 0] <- "inadmissible"
    phase[g_hi < 0 & p_lo > 0] <- "allocation-sensitive"
    q <- 1 - fit$spec$p
    o <- order(q)
    thresholds <- .threshold_set(q[o], g_lo[o], g_hi[o], p_lo[o], p_hi[o],
                                 q_max = 1 - fit$spec$lower)
  }

  band <- data.frame(p = fit$spec$p, q = 1 - fit$spec$p,
                     margin_guaranteed = theta[seq_len(k)],
                     guaranteed_lower = g_lo, guaranteed_upper = g_hi,
                     margin_possible = theta[k + seq_len(k)],
                     possible_lower = p_lo, possible_upper = p_hi,
                     phase = phase)
  structure(
    list(band = band, critical_value = critical, rank = rank, alpha = alpha,
         attempted = R, successful = successes, failure = failure,
         feasibility = feasibility, thresholds = thresholds, fit = fit),
    class = "envelope_band"
  )
}

# Tighten a band for a nonincreasing function of q (values given in
# increasing q order).
.monotone_tighten <- function(lo, hi) {
  lo2 <- rev(cummax(rev(lo)))
  hi2 <- cummin(hi)
  if (any(lo2 > hi2 + 1e-12)) {
    .stop("the monotone confidence band is empty.")
  }
  list(lo = lo2, hi = hi2)
}

# Invert a band for a nonincreasing margin into an interval for the
# missingness rate at which the margin changes sign.
.invert_margin <- function(q, lo, hi, q_max) {
  tight <- .monotone_tighten(lo, hi)
  accepted <- q[tight$lo > 0]
  rejected <- q[tight$hi < 0]
  lower_q <- if (length(accepted)) max(accepted) else 0
  upper_q <- if (length(rejected)) min(rejected) else q_max
  if (lower_q > upper_q) {
    .stop("band inversion produced an empty threshold interval.")
  }
  c(lower = lower_q, upper = upper_q)
}

.threshold_set <- function(q, g_lo, g_hi, p_lo, p_hi, q_max) {
  g <- .invert_margin(q, g_lo, g_hi, q_max)
  pp <- .invert_margin(q, p_lo, p_hi, q_max)
  if (g[["lower"]] > pp[["upper"]]) {
    .stop("the threshold intervals are incompatible with qG <= qP.")
  }
  rbind(
    guaranteed = c(lower = g[["lower"]], upper = min(g[["upper"]], pp[["upper"]])),
    possible = c(lower = max(pp[["lower"]], g[["lower"]]), upper = pp[["upper"]])
  )
}

as.data.frame.envelope_band <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$band
}

print.envelope_band <- function(x, digits = 3, ...) {
  cat("Simultaneous refitted-bootstrap band for the precision envelope\n")
  cat("  Successful refits:", x$successful, "of", x$attempted,
      " (failure policy:", x$failure, ")\n")
  cat("  Critical value:", .fmt(x$critical_value, digits),
      " (rank", x$rank, ", level", .fmt(1 - x$alpha, 2), ")\n")
  if (x$feasibility != "not applicable") {
    cat("  Requirement with complete data:", x$feasibility, "\n")
  }
  cat("\nCertified phase by overall missingness rate:\n")
  b <- x$band
  shown <- data.frame(q = .fmt(b$q, digits), phase = b$phase)
  print(shown, row.names = FALSE, right = TRUE)
  if (!is.null(x$thresholds)) {
    th <- x$thresholds
    cat("\nConfidence intervals for the threshold missingness rates:\n")
    cat("  qG in [", .fmt(th["guaranteed", "lower"], digits), ", ",
        .fmt(th["guaranteed", "upper"], digits), "]\n", sep = "")
    cat("  qP in [", .fmt(th["possible", "lower"], digits), ", ",
        .fmt(th["possible", "upper"], digits), "]\n", sep = "")
  }
  invisible(x)
}
