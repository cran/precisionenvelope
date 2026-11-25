precision_envelope <- function(p, lower, A = NULL, r2 = NULL, b = 1,
                               weights = NULL,
                               covariate = c("continuous", "strata"),
                               target = NULL, variance_limit = NULL) {
  covariate <- match.arg(covariate)
  .check_lower(lower)
  p <- .check_rates(p, lower)
  .check_residual_law(b)
  w <- .normalize_weights(length(b), weights)
  B <- sum(w * b)
  if (B <= 0) {
    .stop("the residual variance law `b` must have a positive mean.")
  }
  if (is.null(A) == is.null(r2)) {
    .stop("supply exactly one of `A` and `r2`.")
  }
  if (!is.null(r2)) {
    if (!.is_number(r2) || r2 < 0 || r2 >= 1) {
      .stop("`r2` must be a single number in [0, 1).")
    }
    A <- r2 / (1 - r2) * B
  } else if (!.is_number(A) || A < 0) {
    .stop("`A` must be a single non-negative number.")
  }
  if (covariate == "strata" && length(b) > .max_strata) {
    .stop("`covariate = \"strata\"` allows at most ", .max_strata,
          " strata; use `covariate = \"continuous\"` for a sample of units.")
  }

  V_full <- A + B
  c_fav <- vapply(p, function(pp) .c_favorable(b, w, pp, lower), numeric(1))
  c_adv <- vapply(p, function(pp) .c_adverse(b, w, pp, lower, covariate), numeric(1))
  c_mcar <- B / p

  tab <- data.frame(
    p = p,
    q = 1 - p,
    C_favorable = c_fav,
    C_mcar = c_mcar,
    C_adverse = c_adv,
    rho_favorable = V_full / (A + c_fav),
    rho_mcar = V_full / (A + c_mcar),
    rho_adverse = V_full / (A + c_adv)
  )
  tab$se_inflation_favorable <- 1 / sqrt(tab$rho_favorable) - 1
  tab$se_inflation_mcar <- 1 / sqrt(tab$rho_mcar) - 1
  tab$se_inflation_adverse <- 1 / sqrt(tab$rho_adverse) - 1
  tab$q_equivalent_favorable <- 1 - B / c_fav
  tab$q_equivalent_adverse <- 1 - B / c_adv

  req <- .requirement(target, variance_limit, A, B)
  thresholds <- NULL
  if (!is.null(req)) {
    creq <- req$C_required
    tab$margin_guaranteed <- creq - c_adv
    tab$margin_possible <- creq - c_fav
    feasible <- creq >= B
    tab$phase <- if (feasible) {
      .phase(tab$margin_guaranteed, tab$margin_possible)
    } else {
      rep("infeasible", nrow(tab))
    }
    thresholds <- .envelope_thresholds(b, w, A, lower, covariate, creq)
    req$rho_required <- V_full / (A + creq)
    req$complete_data_feasible <- feasible
  }

  structure(
    list(table = tab, A = A, B = B, V_full = V_full, lower = lower,
         q_max = 1 - lower, covariate = covariate, requirement = req,
         thresholds = thresholds, b = b, weights = w),
    class = "precision_envelope"
  )
}

# Missingness rates at which the guaranteed, MCAR and possible margins reach
# zero. Each margin is nonincreasing in the missingness rate q.
.envelope_thresholds <- function(b, w, A, lower, covariate, creq) {
  B <- sum(w * b)
  q_max <- 1 - lower
  margin <- list(
    guaranteed = function(q) creq - .c_adverse(b, w, 1 - q, lower, covariate),
    mcar = function(q) creq - B / (1 - q),
    possible = function(q) creq - .c_favorable(b, w, 1 - q, lower)
  )
  vapply(margin, function(f) {
    if (f(0) < 0) {
      return(NA_real_)
    }
    if (f(q_max) >= 0) {
      return(q_max)
    }
    stats::uniroot(f, c(0, q_max), tol = 1e-10)$root
  }, numeric(1))
}

residual_law_gamma <- function(cv, n = 1000L, mean = 1) {
  if (!.is_number(cv) || cv < 0) {
    .stop("`cv` must be a single non-negative number.")
  }
  if (!.is_number(n) || n < 1 || n != round(n)) {
    .stop("`n` must be a positive whole number.")
  }
  if (!.is_number(mean) || mean <= 0) {
    .stop("`mean` must be a single positive number.")
  }
  if (cv == 0) {
    return(rep(mean, n))
  }
  u <- (seq_len(n) - 0.5) / n
  x <- stats::qgamma(u, shape = 1 / cv^2, scale = cv^2)
  x * mean / base::mean(x)
}

as.data.frame.precision_envelope <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
}

print.precision_envelope <- function(x, digits = 3, ...) {
  cat("Sharp precision envelope\n")
  cat("  Overlap lower bound:", .fmt(x$lower, digits),
      "  (missingness rates up to", .fmt(x$q_max, digits), "are admissible)\n")
  cat("  Variance components: A =", format(signif(x$A, 4)),
      " B =", format(signif(x$B, 4)),
      "  R2 =", .fmt(x$A / x$V_full, digits), "\n")
  cat("  Covariate space:", x$covariate, "\n")
  req <- x$requirement
  if (!is.null(req)) {
    if (req$type == "retention") {
      cat("  Requirement: information retention of at least", .fmt(req$value, digits), "\n")
    } else {
      cat("  Requirement: efficient variance of at most", format(signif(req$value, 4)),
          "(retention of at least", .fmt(req$rho_required, digits), ")\n")
    }
  }
  cat("\n")
  tab <- x$table
  shown <- data.frame(
    q = .fmt(tab$q, digits),
    favorable = .fmt(tab$rho_favorable, digits),
    MCAR = .fmt(tab$rho_mcar, digits),
    adverse = .fmt(tab$rho_adverse, digits),
    "SE increase" = paste0(.fmt(100 * tab$se_inflation_favorable, 1), "% to ",
                           .fmt(100 * tab$se_inflation_adverse, 1), "%"),
    check.names = FALSE
  )
  if (!is.null(tab$phase)) {
    shown$phase <- tab$phase
  }
  cat("Information retention by overall missingness rate q:\n")
  print(shown, row.names = FALSE, right = TRUE)
  if (!is.null(x$thresholds)) {
    th <- x$thresholds
    cat("\nThreshold missingness rates for the requirement:\n")
    if (all(is.na(th))) {
      cat("  The requirement cannot be met even without missing data.\n")
    } else {
      cat("  guaranteed up to q =", .fmt(th[["guaranteed"]], digits), "\n")
      cat("  met under MCAR up to q =", .fmt(th[["mcar"]], digits), "\n")
      cat("  attainable by some admissible mechanism up to q =",
          .fmt(th[["possible"]], digits), "\n")
    }
  }
  invisible(x)
}

plot.precision_envelope <- function(x, y, n = 201L, ...) {
  q <- seq(0, x$q_max, length.out = n)
  p <- 1 - q
  b <- x$b
  w <- x$weights
  A <- x$A
  V <- x$V_full
  c_fav <- vapply(p, function(pp) .c_favorable(b, w, pp, x$lower), numeric(1))
  c_adv <- vapply(p, function(pp) .c_adverse(b, w, pp, x$lower, x$covariate), numeric(1))
  rho_fav <- V / (A + c_fav)
  rho_adv <- V / (A + c_adv)
  rho_mcar <- V / (A + x$B / p)

  old <- graphics::par(mar = c(4.5, 4.5, 6, 1))
  on.exit(graphics::par(old), add = TRUE)
  y_low <- min(rho_adv, if (!is.null(x$requirement)) x$requirement$rho_required)
  graphics::plot(q, rho_fav, type = "n", ylim = c(max(0, y_low - 0.05), 1),
                 xlab = "Overall missingness rate", ylab = "Information retention",
                 las = 1, ...)
  graphics::polygon(c(q, rev(q)), c(rho_fav, rev(rho_adv)),
                    col = grDevices::gray(0.92), border = NA)
  graphics::lines(q, rho_adv, lwd = 2, col = "black")
  graphics::lines(q, rho_mcar, lwd = 2, lty = 2, col = grDevices::gray(0.45))
  graphics::lines(q, rho_fav, lwd = 2, lty = 3, col = "black")
  labels <- c("Favorable mechanism", "MCAR", "Adverse mechanism")
  ltys <- c(3, 2, 1)
  cols <- c("black", grDevices::gray(0.45), "black")
  if (!is.null(x$requirement)) {
    graphics::abline(h = x$requirement$rho_required, lty = 4, col = grDevices::gray(0.25))
    labels <- c(labels, "Requirement")
    ltys <- c(ltys, 4)
    cols <- c(cols, grDevices::gray(0.25))
    th <- x$thresholds
    th <- th[is.finite(th) & th < x$q_max]
    if (length(th) > 0L) {
      graphics::abline(v = th, lty = 3, col = grDevices::gray(0.55))
      names_map <- c(guaranteed = "qG", mcar = "qMCAR", possible = "qP")
      graphics::axis(3, at = th, labels = names_map[names(th)],
                     tick = TRUE, cex.axis = 0.8)
    }
  }
  # Place the legend in the top margin, about three lines above the plot
  # region, so that it never overlaps a curve or the threshold labels.
  usr <- graphics::par("usr")
  inch_per_line <- graphics::par("mai")[3] / graphics::par("mar")[3]
  y_per_inch <- diff(usr[3:4]) / graphics::par("pin")[2]
  graphics::legend(x = mean(usr[1:2]), y = usr[4] + 3.4 * inch_per_line * y_per_inch,
                   legend = labels, lty = ltys, col = cols, lwd = 2, ncol = 2,
                   xjust = 0.5, yjust = 0.5, bty = "n", cex = 0.8, xpd = NA)
  invisible(data.frame(q = q, rho_favorable = rho_fav, rho_mcar = rho_mcar,
                       rho_adverse = rho_adv))
}
