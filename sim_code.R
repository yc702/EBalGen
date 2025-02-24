## Evaluate CI coverage for exact balancing
library(halfmoon)
library(ggplot2)
library(doParallel)
library(doRNG)
library(resample)
library(dplyr)
source("simu_make_data_readme.R")

## ATE
n <- 800
p <- 5
H_vars <- 1:5

setting_x =1
setting_src = 1
setting_trt = 1
setting_main = 1
setting_cate = 1
err_sd =1

set.seed(1400, kind = "L'Ecuyer-CMRG")
cl <- makeCluster(5)
registerDoParallel(cl)
ATET <- foreach(1:2000, .combine = "mean",.packages = c("dplyr","MASS","resample")) %dorng% {
  with(make_data(5e5, p, setting_x = setting_x,
                 setting_src = setting_src, setting_trt=setting_trt,
                 setting_main = setting_main, setting_cate = setting_cate), mean(y1[s != 1] - y0[s != 1]))
}
on.exit(stopCluster(cl))

ATET = mean(ATET)


## CI construction for RPM_CI
combine_output <- foreach(1:100, .combine = rbind,.packages = c("resample"),.errorhandling = 'remove') %dorng%{

  dat <- make_data(n, p, setting_x = setting_x,
                   setting_src = setting_src, setting_trt=setting_trt,
                   setting_main = setting_main, setting_cate = setting_cate)

  xs <- dat$x[dat$s == 1,]
  xt <- dat$x[dat$s != 1,]
  trts <- dat$trt[dat$s == 1]
  ys <- dat$y[dat$s == 1]
  delta = numeric(length(H_vars)+ncol(xs))
  target_moments = colMeans(xt)[H_vars]
  target_sd = colStdevs(xt)[H_vars]

  ATE_CI = RPM_CI(xs, ys,trts, H_vars=H_vars,target_mean=target_moments,
                  target_sd=target_sd,num_sim=100, H_add_intercept=TRUE,
                  cluster=5, set_seed=100)
  return(unlist(ATE_CI))
}




## Evaluate CI coverage for approximate balancing
n <- 400
p <- 5
H_vars <- 1:3
setting_x =2

cl <- makeCluster(5)
registerDoParallel(cl)
ATET <- foreach(1:2000, .combine = "mean",.packages = c("dplyr","MASS","resample")) %dorng% {
  with(make_data(5e5, p, setting_x = setting_x,
                 setting_src = setting_src, setting_trt=setting_trt,
                 setting_main = setting_main, setting_cate = setting_cate), mean(y1[s != 1] - y0[s != 1]))
}
on.exit(stopCluster(cl))

ATET = mean(ATET)


## CI construction for RPM_AB
combine_output <- foreach(1:100, .combine = rbind,.packages = c("resample"),.errorhandling = 'remove') %dorng%{

  dat <- make_data(n, p, setting_x = setting_x,
                   setting_src = setting_src, setting_trt=setting_trt,
                   setting_main = setting_main, setting_cate = setting_cate)

  xs <- dat$x[dat$s == 1,]
  xt <- dat$x[dat$s != 1,]
  trts <- dat$trt[dat$s == 1]
  ys <- dat$y[dat$s == 1]
  delta = numeric(length(H_vars)+ncol(xs))
  target_moments = colMeans(xt)[H_vars]
  target_sd = colStdevs(xt)[H_vars]

  ATE_CI = RPM_AB(xs, ys,trts, H_vars=H_vars,target_mean=target_moments,
                  target_sd=target_sd,num_sim=100, H_add_intercept=TRUE,
                  cluster=5, set_seed=100)
  return(unlist(ATE_CI))
}


## Compare with EB for standard causal inference with package weightit
## ATE
n <- 800
p <- 5
H_vars <- NULL
set.seed(1400, kind = "L'Ecuyer-CMRG")

setting_x =1
setting_src = 1
setting_trt = 1
setting_main = 1
setting_cate = 1
err_sd =1

dat <- make_data(n, p, setting_x = setting_x,
                 setting_src = setting_src, setting_trt=setting_trt,
                 setting_main = setting_main, setting_cate = setting_cate)

xs <- dat$x[dat$s == 1,]
xt <- dat$x[dat$s != 1,]
trts <- dat$trt[dat$s == 1]
ys <- dat$y[dat$s == 1]
delta = numeric(length(H_vars)+ncol(xs))
target_moments = NULL

wts_gen <- ebal_wts(xs, trts,H_vars, target_moments,H_add_intercept = TRUE,delta )$w
wts_gen <- wts_gen %>% dplyr::as_tibble() %>%
  dplyr::mutate(trts = trts) %>%
  dplyr::group_by(trts) %>%
  dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-trts)
head(wts_gen$value,10)

colSums((2 * trts - 1) * ys * wts_gen)


## Entropy balancing using weightit
library(WeightIt)
xs = data.frame(xs)
colnames(xs) = c("X1","X2","X3","X4","X5")
W1 <- weightit(trts ~ X1 + X2 + X3 +
                 X4 + X5, data = xs,
               method = "ebal", estimand = "ATE")$weights

W1 <- W1 %>% dplyr::as_tibble() %>%
  dplyr::mutate(trts = trts) %>%
  dplyr::group_by(trts) %>%
  dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-trts)
head(W1$value,10)

colSums((2 * trts - 1) * ys * W1)

mean((wts_gen$value - W1$value)^2)



## Compare with SBW
n <- 800
p <- 5
H_vars <- 1:5
set.seed(1400, kind = "L'Ecuyer-CMRG")

setting_x =1
setting_src = 1
setting_trt = 1
setting_main = 1
setting_cate = 1
err_sd =1
dat <- make_data(n, p, setting_x = setting_x,
                 setting_src = setting_src, setting_trt=setting_trt,
                 setting_main = setting_main, setting_cate = setting_cate)

xs <- dat$x[dat$s == 1,]
xt <- dat$x[dat$s != 1,]
trts <- dat$trt[dat$s == 1]
ys <- dat$y[dat$s == 1]
delta = numeric(length(H_vars)+ncol(xs))
target_moments = colMeans(xt)[H_vars]
target_sd = colStdevs(xt)[H_vars]
wts_gen <- ebal_wts(xs, trts,H_vars, target_moments,H_add_intercept = TRUE,delta )$w
wts_gen <- wts_gen %>% dplyr::as_tibble() %>%
  dplyr::mutate(trts = trts) %>%
  dplyr::group_by(trts) %>%
  dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-trts)
head(wts_gen$value,10)

colSums((2 * trts - 1) * ys * wts_gen)


library(sbw)
xs <- data.frame(xs)
colnames(xs) = c("X1","X2","X3","X4","X5")
dat.1 <- xs[trts==1,]
dat.0 <- xs[trts==0,]
bal_cov <- c("X1","X2","X3","X4","X5")
bal_alg <- FALSE
bal_tol <- tols <-  rep(0,5)
bal_std <- "manual"
bal <- list(bal_cov = bal_cov, bal_alg = bal_alg, bal_tol = bal_tol, bal_std = bal_std)
wei <- list(wei_sum = TRUE, wei_pos = TRUE)
par_tar <- colMeans(xt)

sbw.results.1 <- sbw(dat = dat.1, bal = bal,
                     par = list(par_est = "aux", par_tar= par_tar),
                     sol = list(sol_nam = "gurobi"),wei = wei)

sbw.results.0 <- sbw(dat = dat.0, bal = bal,
                     par = list(par_est = "aux", par_tar= par_tar),
                     sol = list(sol_nam = "gurobi"),wei = wei)

weighted.df.1 <- sbw.results.1$dat_weights
weighted.df.0 <- sbw.results.0$dat_weights
weighted.df <- rbind(weighted.df.1, weighted.df.0)
weighted.df <- weighted.df[order(as.numeric(rownames(weighted.df))), ]
head(weighted.df$sbw_weights,10)
sum((2 * trts - 1) * ys * weighted.df$sbw_weights)
mean((wts_gen$value - weighted.df$sbw_weights)^2)




## Change another tols/delta
n <- 800
p <- 5
H_vars <- 1:5
set.seed(1400, kind = "L'Ecuyer-CMRG")

setting_x =1
setting_src = 1
setting_trt = 1
setting_main = 1
setting_cate = 1
err_sd =1
dat <- make_data(n, p, setting_x = setting_x,
                 setting_src = setting_src, setting_trt=setting_trt,
                 setting_main = setting_main, setting_cate = setting_cate)

xs <- dat$x[dat$s == 1,]
xt <- dat$x[dat$s != 1,]
trts <- dat$trt[dat$s == 1]
ys <- dat$y[dat$s == 1]
delta = numeric(length(H_vars)+ncol(xs))+0.05
target_moments = colMeans(xt)[H_vars]
target_sd = colStdevs(xt)[H_vars]
wts_gen <- ebal_wts(xs, trts,H_vars, target_moments,H_add_intercept = TRUE,delta )$w
wts_gen <- wts_gen %>% dplyr::as_tibble() %>%
  dplyr::mutate(trts = trts) %>%
  dplyr::group_by(trts) %>%
  dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
  dplyr::ungroup() %>%
  dplyr::select(-trts)
head(wts_gen$value,10)

colSums((2 * trts - 1) * ys * wts_gen)


library(sbw)
xs <- data.frame(xs)
colnames(xs) = c("X1","X2","X3","X4","X5")
dat.1 <- xs[trts==1,]
dat.0 <- xs[trts==0,]
bal_cov <- c("X1","X2","X3","X4","X5")
bal_alg <- FALSE
bal_tol <- tols <-  rep(0,5)+0.05
bal_std <- "manual"
bal <- list(bal_cov = bal_cov, bal_alg = bal_alg, bal_tol = bal_tol, bal_std = bal_std)
wei <- list(wei_sum = TRUE, wei_pos = TRUE)
par_tar <- colMeans(xt)

sbw.results.1 <- sbw(dat = dat.1, bal = bal,
                     par = list(par_est = "aux", par_tar= par_tar),
                     sol = list(sol_nam = "gurobi"),wei = wei)

sbw.results.0 <- sbw(dat = dat.0, bal = bal,
                     par = list(par_est = "aux", par_tar= par_tar),
                     sol = list(sol_nam = "gurobi"),wei = wei)

weighted.df.1 <- sbw.results.1$dat_weights
weighted.df.0 <- sbw.results.0$dat_weights
weighted.df <- rbind(weighted.df.1, weighted.df.0)
weighted.df <- weighted.df[order(as.numeric(rownames(weighted.df))), ]
head(weighted.df$sbw_weights,10)
sum((2 * trts - 1) * ys * weighted.df$sbw_weights)

