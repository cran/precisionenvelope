envelope_fit <- function(data, outcome, covariates, observed = NULL, p, lower,
                         target = NULL, variance_limit = NULL, folds = 2L,
                         learners = envelope_learners(), variance_floor = NULL,
                         propensity_clip = c(0.03, 0.995),
                         min_observed_training = 20L) {
  if (!is.data.frame(data)) {
    .stop("`data` must be a data frame.")
  }
  if (!is.character(outcome) || length(outcome) != 1L || !outcome %in% names(data)) {
    .stop("`outcome` must name one column of `data`.")
  }
  if (!is.character(covariates) || length(covariates) == 0L ||
      !all(covariates %in% names(data))) {
    .stop("`covariates` must name one or more columns of `data`.")
  }
  if (!is.null(observed) &&
      (!is.character(observed) || length(observed) != 1L || !observed %in% names(data))) {
    .stop("`observed` must be NULL or name one column of `data`.")
  }
  .check_lower(lower)
  p <- .check_rates(p, lower)
  if (is.unsorted(p, strictly = TRUE)) {
    .stop("`p` must be strictly increasing.")
  }
  if (!.is_number(folds) || folds < 2 || folds != round(folds)) {
    .stop("`folds` must be a whole number of at least 2.")
  }
  if (!is.numeric(propensity_clip) || length(propensity_clip) != 2L ||
      !(0 < propensity_clip[1] && propensity_clip[1] < propensity_clip[2] &&
        propensity_clip[2] <= 1)) {
    .stop("`propensity_clip` must be two increasing numbers in (0, 1].")
  }
  if (!.is_number(min_observed_training) || min_observed_training < 4) {
    .stop("`min_observed_training` must be at least 4.")
  }
  learners <- do.call(envelope_learners, as.list(learners))

  x <- as.matrix(data[, covariates, drop = FALSE])
  if (!is.numeric(x) || any(!is.finite(x))) {
    .stop("covariates must be numeric and completely observed.")
  }
  y <- data[[outcome]]
  if (!is.numeric(y)) {
    .stop("the outcome must be numeric.")
  }
  if (is.null(observed)) {
    r <- as.integer(!is.na(y))
  } else {
    r_raw <- data[[observed]]
    if (is.logical(r_raw)) {
      r_raw <- as.integer(r_raw)
    }
    if (!is.numeric(r_raw) || any(!is.finite(r_raw)) || !all(r_raw %in% c(0, 1))) {
      .stop("the observation indicator must contain only 0 and 1.")
    }
    r <- as.integer(r_raw)
  }
  if (any(r == 1L & !is.finite(y))) {
    .stop("the outcome must be finite wherever it is marked as observed.")
  }
  if (length(unique(r)) < 2L) {
    .stop("the data must contain both observed and missing outcomes.")
  }
  y[r == 0L] <- 0

  if (is.null(variance_floor)) {
    variance_floor <- 1e-3 * stats::var(y[r == 1L])
  }
  if (!.is_number(variance_floor) || variance_floor <= 0) {
    .stop("`variance_floor` must be a single positive number.")
  }

  spec <- list(p = p, lower = lower, target = target,
               variance_limit = variance_limit, folds = as.integer(folds),
               variance_floor = variance_floor, propensity_clip = propensity_clip,
               min_observed_training = min_observed_training)
  point <- .envelope_point(x, r, y, spec, learners)

  structure(
    c(point,
      list(n = length(r), observed_count = sum(r), observation_rate = mean(r),
           share_below_lower = mean(point$nuisance$e < lower),
           spec = spec, learners = learners,
           data = list(x = x, r = r, y = y),
           variables = list(outcome = outcome, covariates = covariates))),
    class = "envelope_fit"
  )
}

# Cross-fitted nuisance predictions: outcome mean m(W), residual variance
# v(W) = Var(Y | W) and observation probability e(W).
.crossfit <- function(x, r, y, spec, learners) {
  n <- length(r)
  fold <- integer(n)
  fold[sample.int(n)] <- .chunks(n, spec$folds)
  m <- v <- e <- numeric(n)
  for (k in seq_len(spec$folds)) {
    train <- which(fold != k)
    test <- which(fold == k)
    obs <- train[r[train] == 1L]
    if (length(obs) < spec$min_observed_training) {
      .stop("too few observed outcomes in a cross-fitting training fold.")
    }
    if (length(unique(r[train])) < 2L) {
      .stop("the observation model has only one class in a training fold.")
    }
    m[test] <- learners$mean(x[obs, , drop = FALSE], y[obs], x[test, , drop = FALSE])

    # Squared residuals for the variance model come from an inner two-way
    # split, so that no residual is computed from its own fitted mean.
    inner <- integer(length(obs))
    inner[sample.int(length(obs))] <- .chunks(length(obs), 2L)
    sq <- numeric(length(obs))
    for (j in 1:2) {
      ev <- which(inner == j)
      tr <- which(inner != j)
      if (length(tr) < 2L || length(ev) < 1L) {
        .stop("too few observed outcomes for the inner variance split.")
      }
      fitted <- learners$mean(x[obs[tr], , drop = FALSE], y[obs[tr]],
                              x[obs[ev], , drop = FALSE])
      sq[ev] <- (y[obs[ev]] - fitted)^2
    }
    v[test] <- pmax(learners$variance(x[obs, , drop = FALSE], sq,
                                      x[test, , drop = FALSE]),
                    spec$variance_floor)
    e[test] <- pmin(pmax(learners$propensity(x[train, , drop = FALSE], r[train],
                                             x[test, , drop = FALSE]),
                         spec$propensity_clip[1]),
                    spec$propensity_clip[2])
  }
  if (any(!is.finite(c(m, v, e)))) {
    .stop("a nuisance learner returned a non-finite prediction.")
  }
  list(m = m, v = v, e = e)
}

.envelope_point <- function(x, r, y, spec, learners) {
  nuisance <- .crossfit(x, r, y, spec, learners)
  .envelope_from_nuisance(r, y, nuisance, spec)
}

# One-step estimates and influence values given cross-fitted nuisances.
.envelope_from_nuisance <- function(r, y, nuisance, spec) {
  m <- nuisance$m
  v <- nuisance$v
  e <- nuisance$e
  n <- length(r)
  lower <- spec$lower
  res <- ifelse(r == 1L, y - m, 0)

  psi_mu <- m + r / e * res
  mu <- mean(psi_mu)
  psi_a <- (m - mu)^2 + 2 * r / e * (m - mu) * res
  psi_b <- v + r / e * (res^2 - v)
  psi_c <- v / e + r / e^2 * (res^2 - v) - v / e^2 * (r - e)
  A <- mean(psi_a)
  B <- mean(psi_b)
  C <- mean(psi_c)
  V_full <- A + B
  V_obs <- A + C
  components <- c(A = A, B = B, C = C)
  if (any(components <= 0)) {
    bad <- names(components)[components <= 0]
    .stop("the one-step estimate of ", paste(bad, collapse = " and "),
          " is not positive (", paste(sprintf("%s = %.4g", names(components), components),
                                      collapse = ", "), "). ",
          if ("A" %in% bad) {
            paste("The covariates may explain almost none of the outcome variance;",
                  "check the mean learner, or plan with precision_envelope(A = 0).")
          })
  }
  da <- psi_a - A
  db <- psi_b - B
  dc <- psi_c - C
  rho <- V_full / V_obs
  drho <- (da + db) / V_obs - V_full * (da + dc) / V_obs^2
  rho_se <- stats::sd(drho) / sqrt(n)

  req <- .requirement(spec$target, spec$variance_limit, A, B)

  w <- rep(1 / n, n)
  k <- length(spec$p)
  c_fav <- c_adv <- numeric(k)
  for (i in seq_len(k)) {
    pp <- spec$p[i]
    alloc <- .favorable_allocation(v, w, pp, lower)
    psi_minus <- switch(
      alloc$branch,
      lower_endpoint = v / lower + r / (e * lower) * (res^2 - v),
      saturated = psi_b,
      zero_information = psi_b,
      v / alloc$pi + alloc$lambda * (alloc$pi - pp) + r / (e * alloc$pi) * (res^2 - v)
    )
    a <- (1 - pp) / (1 - lower)
    threshold <- stats::quantile(v, 1 - a, type = 7, names = FALSE)
    mult <- 1 / lower - 1
    psi_plus <- v + mult * (pmax(v - threshold, 0) + a * threshold) +
      r / e * (1 + mult * (v > threshold)) * (res^2 - v)
    c_fav[i] <- mean(psi_minus)
    c_adv[i] <- mean(psi_plus)
  }

  tab <- data.frame(
    p = spec$p,
    q = 1 - spec$p,
    C_favorable = c_fav,
    C_adverse = c_adv,
    rho_favorable = V_full / (A + c_fav),
    rho_adverse = V_full / (A + c_adv)
  )
  theta <- numeric(0)
  if (!is.null(req)) {
    creq <- req$C_required
    tab$margin_guaranteed <- creq - c_adv
    tab$margin_possible <- creq - c_fav
    tab$phase <- .phase(tab$margin_guaranteed, tab$margin_possible)
    theta <- c(tab$margin_guaranteed, tab$margin_possible)
    names(theta) <- c(paste0("G_", spec$p), paste0("P_", spec$p))
    if (req$type == "variance") {
      # The requirement can hold at some rate only if it holds with complete data.
      feasibility <- req$value - V_full
      if (feasibility < 0) {
        tab$phase <- "infeasible"
      }
      theta <- c(theta, complete_data_feasibility = feasibility)
    }
  }

  list(mu = mu, A = A, B = B, C_realized = C, V_full = V_full, V_observed = V_obs,
       rho_realized = rho, rho_se = rho_se,
       rho_ci = rho + c(-1, 1) * stats::qnorm(0.975) * rho_se,
       requirement = req, table = tab, theta = theta,
       nuisance = list(m = m, v = v, e = e))
}

print.envelope_fit <- function(x, digits = 3, ...) {
  cat("Estimated precision envelope for the mean of", x$variables$outcome, "\n")
  cat("  n =", x$n, " observed =", x$observed_count,
      " observation rate =", .fmt(x$observation_rate, digits), "\n")
  cat("  Overlap lower bound:", .fmt(x$spec$lower, digits), "\n")
  cat("  Estimated mean:", format(signif(x$mu, 6)), "\n")
  cat("  Variance components: A =", format(signif(x$A, 4)), " B =", format(signif(x$B, 4)),
      " C (realized) =", format(signif(x$C_realized, 4)), "\n")
  cat("  Realized information retention:", .fmt(x$rho_realized, digits),
      " (95% Wald interval", .fmt(x$rho_ci[1], digits), "to", .fmt(x$rho_ci[2], digits), ")\n")
  if (x$share_below_lower > 0) {
    cat("  Note: a share of", .fmt(x$share_below_lower, digits),
        "of fitted observation probabilities lies below the overlap lower bound.\n")
  }
  if (!is.null(x$requirement)) {
    if (x$requirement$type == "retention") {
      cat("  Requirement: information retention of at least", .fmt(x$requirement$value, digits), "\n")
    } else {
      cat("  Requirement: efficient variance of at most", format(signif(x$requirement$value, 4)), "\n")
    }
  }
  cat("\nEstimated envelope (point estimates):\n")
  tab <- x$table
  shown <- data.frame(p = .fmt(tab$p, digits), q = .fmt(tab$q, digits),
                      favorable = .fmt(tab$rho_favorable, digits),
                      adverse = .fmt(tab$rho_adverse, digits))
  if (!is.null(tab$phase)) {
    shown$phase <- tab$phase
  }
  print(shown, row.names = FALSE, right = TRUE)
  cat("\nUse envelope_bootstrap() for a simultaneous band and certified phases.\n")
  invisible(x)
}

as.data.frame.envelope_fit <- function(x, row.names = NULL, optional = FALSE, ...) {
  x$table
}

plot.envelope_fit <- function(x, y, ...) {
  tab <- x$table
  old <- graphics::par(mar = c(4.5, 4.5, 6, 1))
  on.exit(graphics::par(old), add = TRUE)
  y_low <- min(tab$rho_adverse, x$rho_realized)
  graphics::plot(tab$q, tab$rho_favorable, type = "n",
                 xlim = range(c(0, tab$q, 1 - x$observation_rate)),
                 ylim = c(max(0, y_low - 0.05), 1),
                 xlab = "Overall missingness rate", ylab = "Information retention",
                 las = 1, ...)
  graphics::polygon(c(tab$q, rev(tab$q)), c(tab$rho_favorable, rev(tab$rho_adverse)),
                    col = grDevices::gray(0.92), border = NA)
  graphics::lines(tab$q, tab$rho_adverse, type = "b", pch = 16, lwd = 2)
  graphics::lines(tab$q, tab$rho_favorable, type = "b", pch = 1, lwd = 2, lty = 3)
  graphics::points(1 - x$observation_rate, x$rho_realized, pch = 18, cex = 1.6)
  labels <- c("Favorable endpoint", "Adverse endpoint", "Realized retention")
  ltys <- c(3, 1, NA)
  pchs <- c(1, 16, 18)
  if (!is.null(x$requirement)) {
    rho_req <- x$V_full / (x$A + x$requirement$C_required)
    graphics::abline(h = rho_req, lty = 4, col = grDevices::gray(0.25))
    labels <- c(labels, "Requirement")
    ltys <- c(ltys, 4)
    pchs <- c(pchs, NA)
  }
  usr <- graphics::par("usr")
  inch_per_line <- graphics::par("mai")[3] / graphics::par("mar")[3]
  y_per_inch <- diff(usr[3:4]) / graphics::par("pin")[2]
  graphics::legend(x = mean(usr[1:2]), y = usr[4] + 3.4 * inch_per_line * y_per_inch,
                   legend = labels, lty = ltys, pch = pchs, lwd = 2, ncol = 2,
                   xjust = 0.5, yjust = 0.5, bty = "n", cex = 0.8, xpd = NA)
  invisible(tab)
}
