## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", fig.width = 6.5,
                      fig.height = 4.8)
has_survival <- requireNamespace("survival", quietly = TRUE)

## ----setup--------------------------------------------------------------------
library(precisionenvelope)

## -----------------------------------------------------------------------------
env <- precision_envelope(p = c(0.95, 0.85, 0.75, 0.65), lower = 0.10,
                          r2 = 0.30, target = 0.80)
env

## -----------------------------------------------------------------------------
het <- precision_envelope(p = c(0.95, 0.85, 0.75, 0.65), lower = 0.10,
                          r2 = 0.30, b = residual_law_gamma(cv = 0.8),
                          target = 0.80)
het
plot(het)

## -----------------------------------------------------------------------------
ill <- precision_envelope(p = 0.75, lower = 0.2, A = 1, b = c(0.1, 1.9),
                          target = 0.80)
ill$table[, c("q", "se_inflation_favorable", "se_inflation_mcar",
              "se_inflation_adverse", "phase")]
ill$thresholds

## -----------------------------------------------------------------------------
strata <- precision_envelope(p = 0.75, lower = 0.2, A = 1, b = c(0.1, 1.9),
                             covariate = "strata", target = 0.80)
strata$table[, c("q", "rho_adverse", "se_inflation_adverse")]
strata$thresholds

## -----------------------------------------------------------------------------
plan <- recovery_plan(b = c(0.1, 1.9), baseline = c(1, 0.5), A = 1,
                      target = 0.80)
plan

## ----eval = has_survival------------------------------------------------------
nafld1 <- survival::nafld1
cases <- nafld1[!is.na(nafld1$case.id) & nafld1$id == nafld1$case.id, ]
dat <- data.frame(bmi = as.numeric(cases$bmi), age = as.numeric(cases$age),
                  male = as.integer(cases$male == 1))
dat$age_centered_sq <- (dat$age - mean(dat$age))^2

set.seed(20260901)
fit <- envelope_fit(dat, outcome = "bmi",
                    covariates = c("age", "male", "age_centered_sq"),
                    p = c(0.70, 0.75, 0.80, 0.85, 0.90, 0.97),
                    lower = 0.05, target = 0.80)
fit
plot(fit)

## ----eval = has_survival------------------------------------------------------
set.seed(20260901)
band <- envelope_bootstrap(fit, R = 199)
band

