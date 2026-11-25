# Nuisance learners. A learner is a function(x, y, newx) that fits a model
# to the numeric matrix x and response y and returns predictions for newx.

.scaled_design <- function(x, newx) {
  center <- colMeans(x)
  scale <- sqrt(colMeans(sweep(x, 2L, center)^2))
  scale[!(scale > 1e-12)] <- 1
  train <- cbind(1, sweep(sweep(x, 2L, center), 2L, scale, "/"))
  evaluate <- cbind(1, sweep(sweep(newx, 2L, center), 2L, scale, "/"))
  list(train = train, evaluate = evaluate)
}

learner_linear <- function() {
  function(x, y, newx) {
    d <- .scaled_design(x, newx)
    fit <- stats::lm.fit(d$train, y)
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    drop(d$evaluate %*% coef)
  }
}

learner_logistic <- function() {
  function(x, y, newx) {
    d <- .scaled_design(x, newx)
    fit <- suppressWarnings(stats::glm.fit(d$train, y, family = stats::binomial()))
    coef <- fit$coefficients
    coef[is.na(coef)] <- 0
    eta <- drop(d$evaluate %*% coef)
    stats::plogis(pmin(pmax(eta, -30), 30))
  }
}

envelope_learners <- function(mean = learner_linear(), variance = learner_linear(),
                              propensity = learner_logistic()) {
  learners <- list(mean = mean, variance = variance, propensity = propensity)
  ok <- vapply(learners, is.function, logical(1))
  if (!all(ok)) {
    .stop("each learner must be a function(x, y, newx) returning predictions.")
  }
  learners
}
