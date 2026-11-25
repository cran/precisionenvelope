test_that("the two-stratum recovery example matches the reference implementation", {
  plan <- recovery_plan(b = c(0.1, 1.9), baseline = c(1, 0.5), A = 1, target = 0.8,
                        weights = c(0.5, 0.5))
  expect_identical(plan$status, "planned")
  expect_equal(plan$C_target, 1.5)
  expect_equal(plan$units$final, c(1, 0.6551724137941102), tolerance = 1e-9)
  expect_equal(plan$added_rate, 0.07758620689705509, tolerance = 1e-9)
  expect_equal(plan$added_cost, 0.07758620689705509, tolerance = 1e-9)
  expect_equal(plan$C_achieved, 1.5, tolerance = 1e-9)
  expect_true(plan$C_achieved <= plan$C_target + 1e-9)
  expect_output(print(plan), "Minimum-cost targeted recovery plan")
})

test_that("recovery recognizes a requirement that already holds", {
  plan <- recovery_plan(b = c(1, 1), baseline = c(0.9, 0.9), residual_target = 2)
  expect_identical(plan$status, "already meets the requirement")
  expect_equal(plan$added_rate, 0)
})

test_that("recovery reports a requirement beyond the complete-data limit", {
  plan <- recovery_plan(b = c(1, 2), baseline = c(0.5, 0.5), residual_target = 1)
  expect_identical(plan$status, "infeasible")
})

test_that("costs shift recovery toward cheaper units", {
  cheap <- recovery_plan(b = c(1, 1), baseline = c(0.5, 0.5), residual_target = 1.6,
                         cost = c(1, 4))
  expect_true(cheap$units$final[1] > cheap$units$final[2])
})

test_that("recovery input checks", {
  expect_error(recovery_plan(b = 1, baseline = 0, residual_target = 1), "baseline")
  expect_error(recovery_plan(b = 1, baseline = 0.5), "exactly one")
  expect_error(recovery_plan(b = 1, baseline = 0.5, target = 0.8), "`A`")
})
