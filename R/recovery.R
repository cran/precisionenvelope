recovery_plan <- function(b, baseline, A = NULL, target = NULL,
                          variance_limit = NULL, residual_target = NULL,
                          cost = 1, weights = NULL) {
  .check_residual_law(b)
  n <- length(b)
  if (!is.numeric(baseline) || length(baseline) != n || any(!is.finite(baseline)) ||
      any(baseline <= 0) || any(baseline > 1)) {
    .stop("`baseline` must give one observation probability in (0, 1] per value of `b`.")
  }
  if (length(cost) == 1L) {
    cost <- rep(cost, n)
  }
  if (!is.numeric(cost) || length(cost) != n || any(!is.finite(cost)) || any(cost <= 0)) {
    .stop("`cost` must be positive, either a single value or one value per value of `b`.")
  }
  w <- .normalize_weights(n, weights)
  B <- sum(w * b)

  given <- !c(is.null(target), is.null(variance_limit), is.null(residual_target))
  if (sum(given) != 1L) {
    .stop("supply exactly one of `target`, `variance_limit` and `residual_target`.")
  }
  if (!is.null(residual_target)) {
    if (!.is_number(residual_target) || residual_target < 0) {
      .stop("`residual_target` must be a single non-negative number.")
    }
    c_target <- residual_target
  } else {
    if (!.is_number(A) || A < 0) {
      .stop("`A` must be supplied as a single non-negative number with `target` or `variance_limit`.")
    }
    c_target <- .requirement(target, variance_limit, A, B)$C_required
  }

  residual <- function(e) sum(w * b / e)
  c_complete <- B
  c_baseline <- residual(baseline)
  priority <- b / (cost * baseline^2)

  build <- function(status, final, eta) {
    delta <- final - baseline
    recovery <- numeric(n)
    open <- baseline < 1
    recovery[open] <- delta[open] / (1 - baseline[open])
    units <- data.frame(b = b, baseline = baseline, cost = cost,
                        final = final, recovery_among_missing = recovery,
                        priority = priority)
    structure(
      list(status = status, units = units,
           baseline_rate = sum(w * baseline), final_rate = sum(w * final),
           added_rate = sum(w * delta), added_cost = sum(w * cost * delta),
           C_complete = c_complete, C_baseline = c_baseline,
           C_target = c_target, C_achieved = residual(final), eta = eta),
      class = "recovery_plan"
    )
  }

  tol <- 1e-12
  atol <- tol * max(1, abs(c_complete), abs(c_baseline), abs(c_target))
  if (c_target >= c_baseline - atol) {
    return(build("already meets the requirement", baseline, 0))
  }
  if (c_target < c_complete - atol) {
    out <- build("infeasible", baseline, NA_real_)
    out$C_achieved <- NA_real_
    return(out)
  }
  if (c_target <= c_complete + atol) {
    final <- baseline
    final[b > 0] <- 1
    return(build("requires every informative outcome", final, NA_real_))
  }

  allocation <- function(eta) {
    raw <- ifelse(b > 0, sqrt(eta * b / cost), 0)
    pmin(1, pmax(baseline, raw))
  }
  eta_lo <- 0
  eta_hi <- 1
  while (residual(allocation(eta_hi)) > c_target) {
    eta_hi <- 4 * eta_hi
    if (!is.finite(eta_hi)) {
      .stop("failed to bracket the recovery multiplier.")
    }
  }
  for (i in seq_len(240L)) {
    eta_mid <- 0.5 * (eta_lo + eta_hi)
    if (residual(allocation(eta_mid)) <= c_target) {
      eta_hi <- eta_mid
    } else {
      eta_lo <- eta_mid
    }
    if (eta_hi - eta_lo <= tol * max(1, eta_hi)) {
      break
    }
  }
  build("planned", allocation(eta_hi), eta_hi)
}

print.recovery_plan <- function(x, digits = 3, ...) {
  cat("Minimum-cost targeted recovery plan\n")
  cat("  Status:", x$status, "\n")
  cat("  Residual variance component: baseline", format(signif(x$C_baseline, 4)),
      ", required", format(signif(x$C_target, 4)),
      ", complete data", format(signif(x$C_complete, 4)), "\n")
  if (x$status %in% c("planned", "requires every informative outcome")) {
    cat("  Observation rate:", .fmt(x$baseline_rate, digits), "->",
        .fmt(x$final_rate, digits), "\n")
    missing_rate <- 1 - x$baseline_rate
    share <- if (missing_rate > 0) x$added_rate / missing_rate else 0
    cat("  Share of currently missing outcomes to recover:", .fmt(share, digits), "\n")
    cat("  Expected added cost per unit:", format(signif(x$added_cost, 4)), "\n")
  }
  invisible(x)
}
