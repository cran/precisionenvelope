# Envelope endpoints for a fixed law of the residual variance b(W).
#
# The favorable endpoint minimizes E{b(W) / pi(W)} over observation laws with
# lower <= pi(W) <= 1 and E{pi(W)} = p. The minimizer is the clipped
# square-root rule pi(W) = min{1, max{lower, sqrt(b(W) / lambda)}}.
#
# The adverse endpoint maximizes the same functional. On a covariate space
# without atoms the maximizer places pi = lower on the upper tail of b(W) of
# mass a = (1 - p) / (1 - lower) and pi = 1 elsewhere, so that
# C_adverse = B + (1 / lower - 1) * (integral of the upper a-tail of b).
# When W takes finitely many values that a mechanism cannot split, the
# maximizer is an extreme point of the feasible set and is found by
# enumeration.

.favorable_allocation <- function(b, w, p, lower, tol = 1e-12, maxit = 200L) {
  n <- length(b)
  positive <- b > 0
  positive_mass <- sum(w[positive])
  zero_mass <- 1 - positive_mass
  p_saturation <- lower + (1 - lower) * positive_mass
  p <- min(max(p, lower), 1)

  if (p <= lower + tol) {
    return(list(pi = rep(lower, n), lambda = Inf, branch = "lower_endpoint",
                p_saturation = p_saturation))
  }
  if (positive_mass == 0) {
    return(list(pi = rep(p, n), lambda = 0, branch = "zero_information",
                p_saturation = lower))
  }
  if (p >= p_saturation - tol) {
    pi <- rep(1, n)
    if (zero_mass > 0) {
      z <- (p - positive_mass) / zero_mass
      pi[!positive] <- min(max(z, lower), 1)
    }
    return(list(pi = pi, lambda = 0, branch = "saturated",
                p_saturation = p_saturation))
  }

  allocation_at <- function(lambda) {
    pmin(pmax(sqrt(b / lambda), lower), 1)
  }
  # At lambda = min positive b every positive value is clipped at one, so the
  # mean exceeds p; at lambda = max b / lower^2 every value is clipped at
  # lower, so the mean falls below p.
  log_lo <- log(min(b[positive]))
  log_hi <- log(max(b[positive]) / lower^2)
  final <- 0.5 * (log_lo + log_hi)
  for (i in seq_len(maxit)) {
    mid <- 0.5 * (log_lo + log_hi)
    h <- sum(w * allocation_at(exp(mid))) - p
    if (abs(h) <= tol) {
      final <- mid
      break
    }
    if (h > 0) {
      log_lo <- mid
    } else {
      log_hi <- mid
    }
    final <- 0.5 * (log_lo + log_hi)
  }
  lambda <- exp(final)
  list(pi = allocation_at(lambda), lambda = lambda, branch = "interior",
       p_saturation = p_saturation)
}

.c_favorable <- function(b, w, p, lower) {
  alloc <- .favorable_allocation(b, w, p, lower)
  sum(w * b / alloc$pi)
}

# Integral of the upper tail of mass a of the weighted law of b.
.upper_tail_integral <- function(b, w, a) {
  if (a <= 0) {
    return(0)
  }
  o <- order(b, decreasing = TRUE)
  bs <- b[o]
  ws <- w[o]
  before <- c(0, cumsum(ws)[-length(ws)])
  take <- pmax(0, pmin(ws, a - before))
  sum(bs * take)
}

.c_adverse_continuous <- function(b, w, p, lower) {
  a <- min(max((1 - p) / (1 - lower), 0), 1)
  sum(w * b) + (1 / lower - 1) * .upper_tail_integral(b, w, a)
}

.max_strata <- 15L

.c_adverse_strata <- function(b, w, p, lower, tol = 1e-12) {
  keep <- w > 0
  b <- b[keep]
  w <- w[keep]
  k <- length(b)
  if (k > .max_strata) {
    .stop("exact enumeration for `covariate = \"strata\"` is limited to ",
          .max_strata, " strata.")
  }
  a <- min(max((1 - p) / (1 - lower), 0), 1)
  wb <- w * b
  best <- -Inf
  # An extreme point has every stratum at pi = lower or pi = 1 except at most
  # one stratum f, whose value is fixed by the mean constraint.
  for (f in seq_len(k)) {
    others <- seq_len(k)[-f]
    m <- length(others)
    if (m > 0L) {
      codes <- seq_len(2L^m) - 1L
      at_lower <- vapply(seq_len(m), function(j) bitwAnd(codes, 2L^(j - 1L)) > 0L,
                         logical(length(codes)))
      at_lower <- matrix(at_lower, nrow = length(codes))
      mass <- drop(at_lower %*% w[others])
      base <- drop(at_lower %*% (wb[others] / lower)) +
        drop((!at_lower) %*% wb[others])
    } else {
      mass <- 0
      base <- 0
    }
    x_f <- (a - mass) / w[f]
    ok <- x_f >= -tol & x_f <= 1 + tol
    if (!any(ok)) {
      next
    }
    x_f <- pmin(pmax(x_f[ok], 0), 1)
    value <- base[ok] + wb[f] / (1 - (1 - lower) * x_f)
    best <- max(best, value)
  }
  if (!is.finite(best)) {
    .stop("no feasible observation law was found for the strata.")
  }
  best
}

.c_adverse <- function(b, w, p, lower, covariate) {
  if (covariate == "strata") {
    .c_adverse_strata(b, w, p, lower)
  } else {
    .c_adverse_continuous(b, w, p, lower)
  }
}
