# The analytic illustration of the article: a residual variance law with two
# equally likely levels 0.1 and 1.9 (so B = 1), A = 1, overlap lower bound
# 0.2 and a retention requirement of 0.80. Reference values were computed
# with the article's reference implementation.

ill <- function(q, ...) {
  precision_envelope(p = 1 - q, lower = 0.2, A = 1, b = c(0.1, 1.9), ...)
}

test_that("favorable, MCAR and adverse endpoints match the analytic illustration", {
  env <- ill(c(0, 0.10, 0.25, 0.40, 0.60, 0.80))
  tab <- env$table
  expect_equal(tab$C_favorable,
               c(1, 1.0125, 1.05, 1.196574912, 1.8333333333, 5), tolerance = 1e-8)
  expect_equal(tab$C_adverse, c(1, 1.95, 3.375, 4.8, 4.9, 5), tolerance = 1e-10)
  expect_equal(tab$C_mcar, 1 / (1 - c(0, 0.10, 0.25, 0.40, 0.60, 0.80)), tolerance = 1e-12)
  expect_equal(tab$rho_favorable[3], 0.9756097561, tolerance = 1e-9)
  expect_equal(tab$rho_mcar[3], 0.8571428571, tolerance = 1e-9)
  expect_equal(tab$rho_adverse[3], 0.4571428571, tolerance = 1e-9)
})

test_that("the 25% example raises the standard error by 1.2% at best", {
  env <- ill(0.25)
  expect_equal(env$table$se_inflation_favorable, 0.012422836565834539, tolerance = 1e-8)
  expect_equal(env$table$se_inflation_mcar, 0.08012344973464325, tolerance = 1e-10)
  expect_equal(env$table$se_inflation_adverse, 0.47901994577490403, tolerance = 1e-10)
})

test_that("threshold missingness rates match the analytic illustration", {
  env <- ill(c(0.05, 0.25, 0.60), target = 0.80)
  expect_equal(env$requirement$C_required, 1.5)
  expect_equal(unname(env$thresholds["guaranteed"]), 1 / 19, tolerance = 1e-8)
  expect_equal(unname(env$thresholds["mcar"]), 1 / 3, tolerance = 1e-8)
  expect_equal(unname(env$thresholds["possible"]), 0.52, tolerance = 1e-8)
  expect_identical(env$table$phase,
                   c("guaranteed", "allocation-sensitive", "inadmissible"))
})

test_that("the envelope always contains MCAR and orders the endpoints", {
  b <- residual_law_gamma(cv = 0.8, n = 400)
  env <- precision_envelope(p = seq(0.3, 1, by = 0.05), lower = 0.1, r2 = 0.3, b = b)
  tab <- env$table
  expect_true(all(tab$C_favorable <= tab$C_mcar + 1e-10))
  expect_true(all(tab$C_mcar <= tab$C_adverse + 1e-10))
  expect_true(all(diff(tab$rho_adverse) >= -1e-12))
  expect_true(all(tab$q_equivalent_favorable <= tab$q + 1e-10))
  expect_true(all(tab$q_equivalent_adverse >= tab$q - 1e-10))
})

test_that("a homoskedastic law makes MCAR the favorable mechanism", {
  env <- precision_envelope(p = c(0.5, 0.8), lower = 0.1, r2 = 0.4)
  expect_equal(env$table$C_favorable, env$table$C_mcar, tolerance = 1e-10)
  a <- (1 - c(0.5, 0.8)) / 0.9
  expect_equal(env$table$C_adverse, 1 + 9 * a, tolerance = 1e-12)
})

test_that("r2 and A give the same envelope", {
  b <- c(0.5, 1, 2)
  e1 <- precision_envelope(p = 0.7, lower = 0.1, r2 = 0.25, b = b)
  e2 <- precision_envelope(p = 0.7, lower = 0.1, A = 0.25 / 0.75 * mean(b), b = b)
  expect_equal(e1$table, e2$table)
})

test_that("indivisible strata give the exact extreme-point bound", {
  b <- c(0.5, 1.0, 3.0)
  w <- c(0.5, 0.3, 0.2)
  strata <- precision_envelope(p = c(0.9, 0.75, 0.6), lower = 0.1, A = 1, b = b,
                               weights = w, covariate = "strata")
  continuous <- precision_envelope(p = c(0.9, 0.75, 0.6), lower = 0.1, A = 1, b = b,
                                   weights = w)
  expect_equal(strata$table$C_adverse, c(1.75, 6.641304347826089, 7.375), tolerance = 1e-10)
  expect_equal(continuous$table$C_adverse, c(4.15, 7.25, 8.75), tolerance = 1e-10)
})

test_that("an absolute requirement that fails with complete data is infeasible", {
  env <- precision_envelope(p = c(0.8, 0.9), lower = 0.1, A = 1, b = 1,
                            variance_limit = 1.5)
  expect_identical(env$table$phase, c("infeasible", "infeasible"))
  expect_true(all(is.na(env$thresholds)))
  expect_false(env$requirement$complete_data_feasible)
})

test_that("invalid inputs are rejected", {
  expect_error(precision_envelope(p = 0.8, lower = 0.1), "exactly one")
  expect_error(precision_envelope(p = 0.05, lower = 0.1, A = 1), "between")
  expect_error(precision_envelope(p = 0.8, lower = 1, A = 1), "lower")
  expect_error(precision_envelope(p = 0.8, lower = 0.1, A = 1, b = -1), "non-negative")
  expect_error(precision_envelope(p = 0.8, lower = 0.1, A = 1, target = 0.8,
                                  variance_limit = 2), "at most one")
})

test_that("print and plot methods run", {
  env <- ill(c(0.1, 0.3), target = 0.8)
  expect_output(print(env), "Sharp precision envelope")
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  curves <- plot(env, n = 21)
  expect_s3_class(curves, "data.frame")
  expect_equal(nrow(curves), 21L)
  expect_s3_class(as.data.frame(env), "data.frame")
})

test_that("residual_law_gamma has the requested mean and spread", {
  b <- residual_law_gamma(cv = 0.5, n = 2000)
  expect_equal(mean(b), 1, tolerance = 1e-12)
  expect_equal(stats::sd(b) / mean(b), 0.5, tolerance = 0.02)
  expect_identical(residual_law_gamma(cv = 0, n = 3, mean = 2), c(2, 2, 2))
})

test_that("allocation by the two region labels reproduces Section 5.1 of the article", {
  env <- precision_envelope(p = 0.75, lower = 0.2, A = 1, b = c(0.1, 1.9),
                            covariate = "strata", target = 0.80)
  expect_equal(env$table$C_adverse, 1.95, tolerance = 1e-10)
  expect_equal(env$table$rho_adverse, 2 / 2.95, tolerance = 1e-10)
  expect_equal(round(100 * env$table$se_inflation_adverse, 1), 21.4)
  expect_equal(unname(env$thresholds["guaranteed"]), (1 - 0.95 / 1.45) / 2, tolerance = 1e-8)
  expect_equal(round(100 * unname(env$thresholds["guaranteed"]), 1), 17.2)
})

test_that("an interval half-width requirement maps to the retention requirement", {
  # Section 5.1: rho0 = 0.80 is equivalent to Vlim = 2.5 and, for n = 1000, to
  # a 95% half-width limit of 0.098.
  vlim <- 1000 * (1.959964 * sqrt(2.5 / 1000) / stats::qnorm(0.975))^2
  a <- precision_envelope(p = 0.75, lower = 0.2, A = 1, b = c(0.1, 1.9), target = 0.8)
  b <- precision_envelope(p = 0.75, lower = 0.2, A = 1, b = c(0.1, 1.9),
                          variance_limit = vlim)
  expect_equal(b$requirement$rho_required, 0.8, tolerance = 1e-6)
  expect_identical(a$table$phase, b$table$phase)
})
