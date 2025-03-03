check_for_mosek <- function() {
  if (!("MOSEK" %in% CVXR::installed_solvers())) {
    testthat::skip("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")
  }
}

test_that("exact ebal_ATE has correct output", {
  set.seed(1400, kind = "L'Ecuyer-CMRG")
  n = 100
  p = 5
  x = runif(n * p)
  x = matrix(4 * x - 2, n, p)
  y = rnorm(n*p)
  trt = rbinom(n,1,0.5)
  H_vars = c(1,2,3)
  target_moments = c(0,0,0)

  ## Exact balancing weights
  wts_full <- ebal_wts(x, trt,H_vars, target_moments = target_moments,
                       H_add_intercept = TRUE,delta=numeric(8))$w
  wts_gen <- wts_full %>% dplyr::as_tibble() %>%
    dplyr::mutate(trts = trt) %>%
    dplyr::group_by(trts) %>%
    dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-trts)

  ATE_wts <- colSums((2 * trt - 1) * y * wts_gen)
  ATE_ebal_EB <- ebal_ATE(x,y,trt,H_vars, target_moments= target_moments,
                          H_add_intercept = TRUE,delta=numeric(8))

  expect_equal(ATE_wts,ATE_ebal_EB$ATE)

})

test_that("approx ebal_ATE has correct output", {
  set.seed(1400, kind = "L'Ecuyer-CMRG")
  n = 100
  p = 5
  x = runif(n * p)
  x = matrix(4 * x - 2, n, p)
  y = rnorm(n*p)
  trt = rbinom(n,1,0.5)
  H_vars = c(1,2,3)
  target_moments = c(0,0,0)
  ## Approximate balancing weights
  check_for_mosek()
  wts_full <- ebal_wts(x, trt,H_vars, target_moments = target_moments,
                       H_add_intercept = TRUE,delta=numeric(8)+0.03)$w
  wts_gen <- wts_full %>% dplyr::as_tibble() %>%
    dplyr::mutate(trts = trt) %>%
    dplyr::group_by(trts) %>%
    dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-trts)

  ATE_wts <- colSums((2 * trt - 1) * y * wts_gen)
  ATE_ebal_AB <- ebal_ATE(x,y,trt,H_vars, target_moments= target_moments,
                          H_add_intercept = TRUE,delta=numeric(8)+0.03)

  expect_equal(ATE_wts,ATE_ebal_AB$ATE)

})
