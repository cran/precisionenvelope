# Cross-implementation check: given fixed nuisance predictions, the one-step
# estimates must agree with the article's reference implementation.

fixture <- utils::read.csv(test_path("fixtures", "nuisance_fixture.csv"))
grid <- c(0.1, 0.6, 0.7, 0.8, 0.9, 1.0)

from_fixture <- function(target = NULL, variance_limit = NULL) {
  r <- fixture$r
  y <- ifelse(r == 1, fixture$y, 0)
  spec <- list(p = grid, lower = 0.1, target = target, variance_limit = variance_limit)
  .envelope_from_nuisance(r, y, list(m = fixture$m, v = fixture$v, e = fixture$e), spec)
}

test_that("one-step estimates agree with the reference implementation", {
  est <- from_fixture(target = 0.8)
  expect_equal(est$mu, 10.812489952161936, tolerance = 1e-10)
  expect_equal(est$A, 5.31868350984671, tolerance = 1e-10)
  expect_equal(est$B, 1.658303859844213, tolerance = 1e-10)
  expect_equal(est$C_realized, 2.136324995969498, tolerance = 1e-10)
  expect_equal(est$rho_realized, 0.9358791964150885, tolerance = 1e-10)
  expect_equal(est$rho_se, 0.012134340240425796, tolerance = 1e-8)
  expect_equal(est$requirement$C_required, 3.402550702266942, tolerance = 1e-10)
  expect_equal(est$table$C_favorable,
               c(16.583038598442126, 2.53174567501295, 2.174331226544248,
                 1.9151348615760133, 1.7362898180879212, 1.658303859844213),
               tolerance = 1e-8)
  expect_equal(est$table$C_adverse,
               c(16.58303859844213, 11.442373237338424, 9.725125692479608,
                 7.38262462302, 5.1374227742596705, 1.658303859844213),
               tolerance = 1e-10)
  expect_equal(unname(est$theta),
               c(-13.180487896175187, -8.039822535071483, -6.322574990212666,
                 -3.9800739207530578, -1.7348720719927284, 1.744246842422729,
                 -13.180487896175183, 0.8708050272539922, 1.2282194757226943,
                 1.4874158406909288, 1.666260884179021, 1.744246842422729),
               tolerance = 1e-8)
  expect_identical(est$table$phase,
                   c("inadmissible", rep("allocation-sensitive", 4), "guaranteed"))
})

test_that("an absolute requirement adds the complete-data feasibility coordinate", {
  est <- from_fixture(variance_limit = 8.5)
  expect_equal(est$requirement$C_required, 3.18131649015329, tolerance = 1e-10)
  expect_length(est$theta, 13L)
  expect_equal(unname(est$theta[13]), 1.5230126303090774, tolerance = 1e-10)
  expect_equal(unname(est$theta[1:6]),
               c(-13.40172210828884, -8.261056747185133, -6.543809202326318,
                 -4.20130813286671, -1.9561062841063803, 1.5230126303090772),
               tolerance = 1e-8)
})

simulate_data <- function(n = 400) {
  w1 <- stats::rnorm(n)
  w2 <- stats::rbinom(n, 1, 0.4)
  y <- 10 + 2 * w1 + 1.5 * w2 + stats::rnorm(n, sd = exp(0.3 * w1))
  r <- stats::rbinom(n, 1, stats::plogis(1 + 0.8 * w1))
  data.frame(y = ifelse(r == 1, y, NA), w1 = w1, w2 = w2)
}

test_that("envelope_fit runs end to end and is reproducible with a seed", {
  set.seed(11)
  dat <- simulate_data()
  set.seed(1)
  fit1 <- envelope_fit(dat, outcome = "y", covariates = c("w1", "w2"),
                       p = c(0.6, 0.8, 0.9), lower = 0.1, target = 0.8)
  set.seed(1)
  fit2 <- envelope_fit(dat, outcome = "y", covariates = c("w1", "w2"),
                       p = c(0.6, 0.8, 0.9), lower = 0.1, target = 0.8)
  expect_s3_class(fit1, "envelope_fit")
  expect_identical(fit1$table, fit2$table)
  expect_true(fit1$rho_realized > 0 && fit1$rho_realized <= 1.2)
  expect_true(all(fit1$table$C_favorable <= fit1$table$C_adverse + 1e-8))
  expect_equal(fit1$observed_count, sum(!is.na(dat$y)))
  expect_output(print(fit1), "Estimated precision envelope")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(fit1), "data.frame")
})

test_that("an explicit observation indicator and custom learners are accepted", {
  set.seed(12)
  dat <- simulate_data(300)
  dat$seen <- as.integer(!is.na(dat$y))
  quadratic <- function(x, y, newx) {
    learner_linear()(cbind(x, x^2), y, cbind(newx, newx^2))
  }
  set.seed(2)
  fit <- envelope_fit(dat, outcome = "y", covariates = c("w1", "w2"), observed = "seen",
                      p = c(0.7, 0.9), lower = 0.1,
                      learners = list(mean = quadratic, variance = quadratic))
  expect_s3_class(fit, "envelope_fit")
  expect_null(fit$requirement)
  expect_true(fit$A > 0)
})

test_that("the refitted bootstrap certifies phases and threshold intervals", {
  set.seed(13)
  dat <- simulate_data(400)
  set.seed(3)
  fit <- envelope_fit(dat, outcome = "y", covariates = c("w1", "w2"),
                      p = c(0.6, 0.7, 0.8, 0.9, 1.0), lower = 0.1, target = 0.8)
  band <- envelope_bootstrap(fit, R = 39)
  expect_s3_class(band, "envelope_band")
  expect_identical(band$successful, 39L)
  expect_identical(band$rank, 38)
  b <- band$band
  expect_true(all(b$guaranteed_lower <= b$margin_guaranteed))
  expect_true(all(b$possible_upper >= b$margin_possible))
  expect_true(all(b$phase %in% c("guaranteed", "allocation-sensitive",
                                 "inadmissible", "indeterminate")))
  th <- band$thresholds
  expect_true(th["guaranteed", "lower"] <= th["guaranteed", "upper"])
  expect_true(th["possible", "lower"] <= th["possible", "upper"])
  expect_true(th["guaranteed", "lower"] <= th["possible", "upper"])
  expect_output(print(band), "Simultaneous refitted-bootstrap band")
})

test_that("envelope_fit input checks", {
  dat <- data.frame(y = c(1, NA, 3), w = c(1, 2, 3))
  expect_error(envelope_fit(dat, "z", "w", p = 0.8, lower = 0.1), "outcome")
  expect_error(envelope_fit(dat, "y", "w", p = c(0.9, 0.8), lower = 0.1), "increasing")
  dat_full <- data.frame(y = c(1, 2, 3), w = c(1, 2, 3))
  expect_error(envelope_fit(dat_full, "y", "w", p = 0.8, lower = 0.1), "both observed and missing")
})

test_that("the public MASLD illustration is reproduced approximately", {
  skip_on_cran()
  skip_if_not_installed("survival")
  nafld1 <- survival::nafld1
  cases <- nafld1[!is.na(nafld1$case.id) & nafld1$id == nafld1$case.id, ]
  dat <- data.frame(bmi = as.numeric(cases$bmi), age = as.numeric(cases$age),
                    male = as.integer(cases$male == 1))
  dat$age_centered_sq <- (dat$age - mean(dat$age))^2
  expect_identical(nrow(dat), 3514L)
  expect_identical(sum(!is.na(dat$bmi)), 2989L)
  set.seed(20260901)
  fit <- envelope_fit(dat, outcome = "bmi", covariates = c("age", "male", "age_centered_sq"),
                      p = c(0.70, 0.75, 0.80, 0.85, 0.90, 0.97), lower = 0.05,
                      target = 0.80, variance_floor = 0.05)
  # Published values: mean 34.24, B = 55.38, C = 65.00, retention 0.858,
  # and at p = 0.85 favorable retention 0.865 and adverse retention 0.177.
  expect_equal(fit$mu, 34.24, tolerance = 0.002)
  expect_equal(fit$B, 55.38, tolerance = 0.03)
  expect_equal(fit$C_realized, 65.00, tolerance = 0.03)
  expect_equal(fit$rho_realized, 0.858, tolerance = 0.01)
  at <- fit$table[fit$table$p == 0.85, ]
  expect_equal(at$rho_favorable, 0.865, tolerance = 0.02)
  # With lower = 0.05 the adverse endpoint weights the upper tail of the fitted
  # residual variance by 1 / 0.05 - 1 = 19, so it moves with the random
  # cross-fitting split more than the other quantities do.
  expect_gt(at$rho_adverse, 0.12)
  expect_lt(at$rho_adverse, 0.26)
  expect_identical(at$phase, "allocation-sensitive")
  expect_identical(fit$table$phase[fit$table$p == 0.70], "inadmissible")
})
