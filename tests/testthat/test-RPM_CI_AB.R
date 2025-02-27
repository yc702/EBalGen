check_for_mosek <- function() {
    if (!("MOSEK" %in% CVXR::installed_solvers())) {
    testthat::skip("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")
    }
}
test_that("exact CI construction has correct output", {
  set.seed(1400, kind = "L'Ecuyer-CMRG")
  n = 100
  p = 5
  x = runif(n * p)
  x = matrix(4 * x - 2, n, p)
  y = rnorm(n*p)
  trt = rbinom(n,1,0.5)
  H_vars = c(1,2,3)
  target_mean = c(0,0,0)
  target_sd=c(1,1,1)

  expect_no_condition({EB_CI <- RPM_CI(x,y,trt,H_vars, target_mean, target_sd, num_sim=100,
                                       H_add_intercept = TRUE,cluster=1, set_seed=111)})
})

test_that("approx CI construction has correct output", {
  set.seed(1400, kind = "L'Ecuyer-CMRG")
  n = 100
  p = 5
  x = runif(n * p)
  x = matrix(4 * x - 2, n, p)
  y = rnorm(n*p)
  trt = rbinom(n,1,0.5)
  H_vars = c(1,2,3)
  target_mean = c(0,0,0)
  target_sd=c(1,1,1)
  check_for_mosek()
  expect_no_condition({AB_CI <- RPM_AB_delta(x,y,trt,H_vars, target_mean, target_sd, num_sim=100,
                                             H_add_intercept = TRUE,delta=numeric(8)+0.01,cluster=1, set_seed=111)})
  # expect_equal(EB_CI$mean_ATE, AB_CI$mean_ATE, tolerance = 0.2)
  # expect_equal(EB_CI$lb_ATE, AB_CI$lb_ATE, tolerance = 0.2)
  # expect_equal(EB_CI$ub_ATE, AB_CI$ub_ATE, tolerance = 0.2)



})
