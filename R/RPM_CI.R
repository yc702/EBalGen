#' @title Confidence interval for exact balancing ATE
#' @description Compute the exact balancing ATE confidence interval (CI)
#' using resampling-based perturbation method proposed in the paper.
#' @param x A data matrix for the source sample. Each column represents
#' source sample covariate and each row represents an observation.
#' @param y A vector of the source sample response values.
#' @param trt A vector of 0, 1 or FALSE/TRUE of treatment assignment for the source sample.
#' @param H_vars A vector of numbers indexing which covariate in `x`
#' need to be balanced between source and target samples.
#' @param target_mean A vector of mean of the target sample covariates
#' that needs to be balanced between source and target.
#' @param target_sd A vector of standard deviation of the target sample covariates
#' that needs to be balanced between source and target.
#' @param num_sim `numeric` the number of simulations used in bootstrap.
#' @param H_add_intercept `logical` whether to include 1 as intercept in
#' H covariates, default as TRUE.
#' @param cluster Number of parallel running CPU cores, Default: 1
#' @param set_seed Random seed for simulation, Default: 111
#' @return A list containing:
#' \describe{
#'  \item{mean_ATE}{The mean ATE over bootstrap.}
#'  \item{lb_ATE}{The lower bound of 95\% CI.}
#'  \item{ub_ATE}{The upper bound of 95\% CI.}
#'  \item{n_success}{The number of feasible solutions in \code{num_sim} bootstrap.}
#' }
#' @import dplyr
#' @import doRNG
#' @import MASS
#' @import resample
#' @import doParallel
#' @import parallel
#' @import foreach
#' @rdname RPM_CI
#' @export
#'
#' @references
#' Chen, R., Chen, G., & Yu, M. (2023). Entropy balancing for causal
#' generalization with target sample summary information. Biometrics, 79(4)
#' @examples
#' library(EBalGen)
#' set.seed(1)
#' n = 100
#' p = 5
#' x = runif(n * p)
#' x = matrix(4 * x - 2, n, p)
#' y = rnorm(n*p)
#' trt = rbinom(n,1,0.5)
#' H_vars = c(1,2,3)
#' target_mean = c(0,0,0)
#' target_sd=c(1,1,1)
#' ## Get CI for this generalized ATE
#' RPM_CI(x,y,trt,H_vars, target_mean, target_sd, num_sim=50,
#' H_add_intercept = TRUE,cluster=1, set_seed=111)
#'
RPM_CI <- function(x,y,trt,H_vars,target_mean,
                   target_sd,num_sim,
                   H_add_intercept=TRUE,
                   cluster=1, set_seed=111){

  if (length(target_mean) > NCOL(x)){
    stop("Error: number of target moments must be smaller than the number of
         source covariates")
  }

  if (length(target_mean) != length(H_vars)){
    stop("Error: number of target moments must be equal to the number of
         specified covariates in the source")
  }

  nx_s <- nrow(x)
  delta = numeric(length(H_vars)+ncol(x))
  ## Parallel processing
  cl <- makeCluster(cluster)
  registerDoParallel(cl)

  if(!is.null(set_seed)){set.seed(set_seed, kind = "L'Ecuyer-CMRG")}

  perturb_CI <- foreach::foreach(i=1:num_sim, .combine = rbind,
                                 .export=c("ebal_wts",".weighted_ATE"),
                                 .packages = c("dplyr","MASS","resample")) %dorng% {

                                   target_moments <- mvrnorm(1,target_mean,diag(target_sd)%*%cor(x[,H_vars])%*%diag(target_sd))

                                   id_source <- sample(nx_s,nx_s,replace = TRUE)
                                   x_s <- x[id_source,]
                                   trts <- trt[id_source]
                                   y_s <- y[id_source]
                                   wts_gen <- tryCatch(
                                     ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,delta)$w,
                                     error = function(e) {
                                       return(rep(NA,nx_s))
                                     }
                                   )
                                   #wts_gen <- ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,delta)$w
                                   ## R function optim() would sometimes give all wts as 0 but declare convergence
                                   #if (round(sum(wts_gen))!=nx_s*2| all(is.na(wts_gen))){
                                   #  wts_gen=rep(NA,nx_s)
                                   #}
                                   .weighted_ATE(y=y_s,trt=trts,wts=wts_gen)


                                 }
  on.exit(stopCluster(cl))

  ## get the CIs and see how many times the ATET is in the CI
  perturb_CI = data.frame(perturb_CI)
  mean_ATE <- mean(perturb_CI$value,na.rm = TRUE)
  lb_ATE <- quantile(perturb_CI$value,0.025,na.rm = TRUE)
  ub_ATE <- quantile(perturb_CI$value,0.975,na.rm = TRUE)
  n_success <- sum(!is.na(perturb_CI$value))

  ## output weights and dual parameters.
  list(mean_ATE=mean_ATE,
       lb_ATE=lb_ATE,
       ub_ATE=ub_ATE,
       n_success=n_success)
}
