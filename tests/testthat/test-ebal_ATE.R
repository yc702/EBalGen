
test_that("ebal_ATE has correct output", {
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


  ## Approximate balancing weights
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
  expect_equal(ATE_ebal_EB, ATE_ebal_AB, tolerance = 0.1)


})
