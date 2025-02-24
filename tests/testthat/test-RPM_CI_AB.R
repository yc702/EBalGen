test_that("CI construction has correct output", {
  n = 200
  p = 5
  x = runif(n * p)
  x = matrix(4 * x - 2, n, p)
  y = rnorm(n*p)
  trt = rbinom(n,1,0.5)
  H_vars = c(1,2,3)
  target_mean = c(0,0,0)
  target_sd=c(1,1,1)


  expect_no_condition({EB_CI <- RPM_CI(x,y,trt,H_vars, target_mean, target_sd, num_sim=50,
                                       H_add_intercept = TRUE,cluster=1, set_seed=111)})
  expect_no_condition({AB_CI <- RPM_AB_delta(x,y,trt,H_vars, target_mean, target_sd, num_sim=50,
                                             H_add_intercept = TRUE,delta=numeric(8)+0.02,cluster=1, set_seed=111)})
  expect_equal(EB_CI$mean_ATE, AB_CI$mean_ATE, tolerance = 0.1)
  expect_equal(EB_CI$lb_ATE, AB_CI$lb_ATE, tolerance = 0.1)
  expect_equal(EB_CI$ub_ATE, AB_CI$ub_ATE, tolerance = 0.1)



})
