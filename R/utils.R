# Internal helpers for argument checking and small numerical tasks.

.stop <- function(...) stop(..., call. = FALSE)

.is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

.check_lower <- function(lower) {
  if (!.is_number(lower) || lower <= 0 || lower >= 1) {
    .stop("`lower` must be a single number strictly between 0 and 1.")
  }
  invisible(lower)
}

.check_rates <- function(p, lower, tol = 1e-10) {
  if (!is.numeric(p) || length(p) == 0L || any(!is.finite(p))) {
    .stop("`p` must be a non-empty numeric vector of finite observation rates.")
  }
  if (any(p < lower - tol) || any(p > 1 + tol)) {
    .stop("every observation rate in `p` must lie between `lower` and 1.")
  }
  # Values such as 1 - 0.8 fall a rounding error below 0.2; clamp them.
  pmin(pmax(p, lower), 1)
}

.normalize_weights <- function(n, weights) {
  if (is.null(weights)) {
    return(rep(1 / n, n))
  }
  if (!is.numeric(weights) || length(weights) != n || any(!is.finite(weights))) {
    .stop("`weights` must be a finite numeric vector with one entry per value of `b`.")
  }
  if (any(weights < 0)) {
    .stop("`weights` must be non-negative.")
  }
  total <- sum(weights)
  if (total <= 0) {
    .stop("`weights` must have a positive sum.")
  }
  weights / total
}

.check_residual_law <- function(b) {
  if (!is.numeric(b) || length(b) == 0L || any(!is.finite(b))) {
    .stop("`b` must be a non-empty numeric vector of finite residual variances.")
  }
  if (any(b < 0)) {
    .stop("`b` must be non-negative.")
  }
  invisible(b)
}

# Translate a precision requirement into the residual budget C_required.
# A retention requirement rho0 asks for V_full / V_observed >= rho0; an
# absolute requirement asks for V_observed <= variance_limit.
.requirement <- function(target, variance_limit, A, B) {
  if (!is.null(target) && !is.null(variance_limit)) {
    .stop("supply at most one of `target` and `variance_limit`.")
  }
  if (!is.null(target)) {
    if (!.is_number(target) || target <= 0 || target > 1) {
      .stop("`target` must be a single retention requirement in (0, 1].")
    }
    return(list(type = "retention", value = target,
                C_required = (A + B) / target - A))
  }
  if (!is.null(variance_limit)) {
    if (!.is_number(variance_limit) || variance_limit <= 0) {
      .stop("`variance_limit` must be a single positive number.")
    }
    return(list(type = "variance", value = variance_limit,
                C_required = variance_limit - A))
  }
  NULL
}

.phase <- function(margin_guaranteed, margin_possible) {
  ifelse(margin_guaranteed >= 0, "guaranteed",
         ifelse(margin_possible >= 0, "allocation-sensitive", "inadmissible"))
}

# Split 1..n into k consecutive chunks of nearly equal size; the first
# n %% k chunks receive one extra element.
.chunks <- function(n, k) {
  sizes <- rep(n %/% k, k) + c(rep(1L, n %% k), rep(0L, k - n %% k))
  rep(seq_len(k), times = sizes)
}

.fmt <- function(x, digits = 3) {
  formatC(x, digits = digits, format = "f")
}
