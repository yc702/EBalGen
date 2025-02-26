test_that("ebal_wts has correct output", {
  set.seed(1400, kind = "L'Ecuyer-CMRG")
  n = 100
  p = 5
  x = runif(n * p)
  x = matrix(4 * x - 2, n, p)
  trt = rbinom(n,1,0.5)
  H_vars = c(1,2,3)
  target_moments = c(0,0,0)

  ## Exact balancing weights
  expect_no_condition({wts_full_EB <- ebal_wts(x, trt,H_vars, target_moments = target_moments,
                                               H_add_intercept = TRUE,delta=numeric(8))$w})
  expect_no_condition({wts_simple_EB <- ebal_wts_simple(x[,H_vars],target_moments = target_moments,
                                                        H_add_intercept = TRUE,delta=numeric(3))$w})

  expect_equal(round(sum(wts_full_EB)),2*n)
  expect_equal(round(sum(wts_simple_EB)),n)
  expect_failure(expect_equal(wts_full_EB, wts_simple_EB))

  target_moments = c(0,0,0,0)
  expect_error(ebal_wts(x, trt,H_vars, target_moments = target_moments,
                        H_add_intercept = TRUE,delta=numeric(8)),
               "Error: number of target moments must be equal to the number of
         specified covariates in the source")


  ## Approximate balancing weights
  target_moments = c(0,0,0)
  expect_no_condition({wts_full_AB <- ebal_wts(x, trt,H_vars, target_moments = target_moments,
                                               H_add_intercept = TRUE,delta=numeric(8)+0.1)$w})
  expect_no_condition({wts_simple_AB <- ebal_wts_simple(x[,H_vars],target_moments = target_moments,
                                                        H_add_intercept = TRUE,delta=numeric(3)+0.03)$w})

  expect_equal(round(sum(wts_full_AB)),2*n)
  expect_equal(round(sum(wts_simple_AB)),n)
  expect_failure(expect_equal(wts_full_AB, wts_simple_AB))
  expect_equal(wts_simple_EB, wts_simple_AB, tolerance = 0.1)

  H_vars = c(1,2,3,4)
  expect_error(ebal_wts(x, trt,H_vars, target_moments = target_moments,
                        H_add_intercept = TRUE,delta=numeric(8)+0.03),
               "Error: number of target moments must be equal to the number of
         specified covariates in the source")

})
