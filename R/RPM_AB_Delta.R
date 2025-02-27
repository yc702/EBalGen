#' @title Confidence interval for approximate balancing
#' @description Compute the approximate balancing confidence interval (CI)
#' using resampling-based perturbation method proposed in the paper.
#' @param x A covariate data matrix for the source sample.
#' @param y A vector of the source sample response values.
#' @param trt A vector of 0, 1 or FALSE/TRUE of treatment assignment for the source sample.
#' @param H_vars A vector of numbers indexing which covariate in \code{x}
#' need to be balanced between source and target samples.
#' @param target_mean A vector of mean of the target sample covariates
#' that needs to be balanced between source and target.
#' @param target_sd A vector of standard deviation of the target sample covariates
#' that needs to be balanced between source and target.
#' @param num_sim `numeric` the number of simulation used for bootstrap.
#' @param H_add_intercept `logical` whether to include 1 as intercept in
#' H covariates, default as TRUE.
#' @param cluster Number of parallel running CPU cores, Default: 1
#' @param set_seed Random seed for simulation, Default: 111
#' @return A list containing:
#' \describe{
#'  \item{mean_ATE}{The mean ATE over bootstrap.}
#'  \item{lb_ATE}{The lower bound of 95\% CI.}
#'  \item{ub_ATE}{The upper bound of 95\% CI.}
#'  \item{use_exact}{The number of times there no need for approximate balancing in the bootstrap}
#'  \item{n_success}{The number of feasible solutions in `num_sim` bootstrap.}
#' }
#' @references
#' Chen, Y., Chen, G., & Yu, M. (2025+). Confidence interval construction for
#' causally generalized estimators with summary-level target data.
#' @examples
#' library(EBalGen)
#' set.seed(1)
#' n = 100
#' p = 5
#' x = runif(n * p)
#' x = matrix(6 * x - 2, n, p)
#' y = rnorm(n*p)
#' trt = rbinom(n,1,0.5)
#' H_vars = c(1,2,3)
#' target_mean = c(0,0,0)
#' target_sd=c(1,1,1)
#' if (requireNamespace(c("dplyr","doRNG","rockchalk","resample","doParallel"),
#'   quietly = TRUE)) {
#'   \dontrun{
#'   RPM_AB(x,y,trt,H_vars, target_mean, target_sd, num_sim=10,
#'   H_add_intercept = TRUE,cluster=1, set_seed=111)
#'}
#' }
#' @rdname RPM_AB
#' @export
#' @importFrom dplyr %>% as_tibble mutate group_by mutate_at ungroup select vars group_cols
#' @import doRNG
#' @importFrom rockchalk mvrnorm
#' @import resample
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom stats cor quantile
#' @import CVXR
RPM_AB <- function(x,y,trt,H_vars,target_mean,
                   target_sd,num_sim,
                   H_add_intercept=TRUE,
                   cluster=1, set_seed=111){

  if (length(target_mean) > NCOL(x)){
    stop("Error: number of target moments must be smaller than the
         number of source covariates")
  }

  if (length(target_mean) != length(H_vars)){
    stop("Error: number of target moments must be equal to the number
         of specified covariates in the source")
  }

  nx_s <- nrow(x)
  sd_s = resample::colStdevs(x)
  delta = numeric(ncol(x)+length(H_vars))

  exact_bal <- tryCatch(
    ebal_wts(x, trt, H_vars, target_moments=target_mean,H_add_intercept = TRUE, delta),
    error = function(e) {
      return(NA)
    }
  )

  # wts_gen <- approx_bal(xs, trts, colMeans(xt)[H_vars],H_add_intercept = TRUE,delta=delta*0)$w
  ## R function optim() would sometimes give all wts as 0 but declare convergence
  #if (round(sum(wts_gen))!=nx_s*2| all(is.na(wts_gen))){
  #  wts_gen=rep(NA,nx_s)
  #}

  if (length(exact_bal)<2){
    delta = c(sd_s[H_vars],sd_s[H_vars],sd_s[-H_vars])
  } else{

    inv_lamb <- abs(1/(exact_bal$theta))[-c(1,length(H_vars)+2)]
    delta <- inv_lamb
  }

  ## Parallel processing
  cl <- makeCluster(cluster)
  registerDoParallel(cl)

  if(!is.null(set_seed)){set.seed(set_seed, kind = "L'Ecuyer-CMRG")}
  perturb_CI <- foreach(1:num_sim, .combine = rbind,
                        .packages = c("dplyr","rockchalk","resample","CVXR"),
                        .export=c("ebal_wts","weighted_ATE"),
                        .errorhandling = 'remove') %dorng% {

                          target_moments <- rockchalk::mvrnorm(1,target_mean,diag(target_sd)%*%stats::cor(x[,H_vars])%*%diag(target_sd))

                          id_source <- sample(nx_s,nx_s,replace = TRUE)
                          x_s <- x[id_source,]
                          trts <- trt[id_source]
                          y_s <- y[id_source]

                          constant = 0

                          #wts_gen <- ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,
                          #                    delta=numeric(length(H_vars)+ncol(x)))$w
                          wts_gen <- tryCatch(
                            ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,
                                     delta=numeric(length(H_vars)+ncol(x)))$w,
                            error = function(e) {
                              return(rep(NA,nx_s))
                            }
                          )
                          ## R function optim() would sometimes give all wts as 0 but declare convergence
                          #if (round(sum(wts_gen))!=nx_s*2| all(is.na(wts_gen))){
                          #  wts_gen=rep(NA,nx_s)
                          #}

                          while(all(is.na(wts_gen))){
                            constant = constant + 0.1
                            #wts_gen <- ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,
                            #                    delta=delta*constant)$w

                            wts_gen <- tryCatch(
                              ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,
                                       delta=delta*constant)$w,
                              error = function(e) {
                                return(rep(NA,nx_s))
                              }
                            )
                          }

                          approx_ATE <- weighted_ATE(y_s,trts,wts_gen)

                          c(approx_ATE,constant)
                        }
  on.exit(stopCluster(cl))
  ## get the CIs and see how many times the ATET is in the CI
  perturb_CI <- data.frame(perturb_CI)
  n_success <- sum(!is.na(perturb_CI[,1]))
  mean_ATE <- mean(perturb_CI[,1],na.rm=TRUE)
  lb_ATE <- stats::quantile(perturb_CI[,1],0.025,na.rm=TRUE)
  ub_ATE <- stats::quantile(perturb_CI[,1],0.975,na.rm=TRUE)
  delta_0 = sum(perturb_CI[,2]==0)

  list(mean_ATE=mean_ATE,lb_ATE=lb_ATE, ub_ATE=ub_ATE,
       n_success=n_success,use_exact=delta_0)

}



#' @title Confidence interval for approximate balancing for a specified margin
#' @description Compute the approximate balancing confidence interval (CI) for a specified margin
#' using resampling-based perturbation method proposed in the paper.
#' @param x A covariate data matrix for the source sample.
#' @param y Logical vector of treatment assignment for the source sample.
#' @param trt A vector of 0, 1 or FALSE/TRUE of treatment assignment for the source sample.
#' @param H_vars A vector of numbers indexing which covariate in \code{x}
#' need to be balanced between source and target samples.
#' @param target_mean A vector of mean of the target sample covariates
#' that needs to be balanced between source and target.
#' @param target_sd A vector of standard deviation of the target sample covariates
#' that needs to be balanced between source and target.
#' @param num_sim `numeric` the number of simulation used for bootstrap.
#' @param H_add_intercept `logical` whether to include 1 as intercept in
#' H covariates, default as TRUE.
#' @param delta A vector specifying the approximate balancing tolerance margin.
#' The vector has a total length of 2H+G, where H represents the number of
#' covariates balanced between the source (treatment and control) and the
#' target moments, and G represents the covariates balanced solely between
#' the source treatment and control groups. If exact balancing, delta are all zeros.
#' @param cluster Number of parallel running CPU cores, Default: 1
#' @param set_seed Random seed for simulation, Default: 111
#' @return A list containing:
#' \describe{
#'  \item{mean_ATE}{The mean ATE over bootstrap.}
#'  \item{lb_ATE}{The lower bound of 95\% CI.}
#'  \item{ub_ATE}{The upper bound of 95\% CI.}
#'  \item{n_success}{The number of feasible solutions in `num_sim` bootstrap.}
#' }
#'
#' @references
#' Chen, Y., Chen, G., & Yu, M. (2025+). Confidence interval construction for
#' causally generalized estimators with summary-level target data.
#' @examples
#' library(EBalGen)
#' set.seed(1)
#' n = 100
#' p = 5
#' x = runif(n * p)
#' x = matrix(6 * x - 2, n, p)
#' y = rnorm(n*p)
#' trt = rbinom(n,1,0.5)
#' H_vars = c(1,2,3)
#' target_mean = c(0,0,0)
#' target_sd=c(1,1,1)
#' delta = numeric(8)+0.1
#' if (requireNamespace(c("dplyr","doRNG","rockchalk","resample","doParallel"),
#'   quietly = TRUE)) {
#'   \dontrun{
#'   RPM_AB_delta(x,y,trt,H_vars, target_mean, target_sd, num_sim=10,
#'   H_add_intercept = TRUE,delta,cluster=1, set_seed=111)
#'}
#' }
#' @rdname RPM_AB_delta
#' @export
#' @importFrom dplyr %>% as_tibble mutate group_by mutate_at ungroup select vars group_cols
#' @import doRNG
#' @importFrom rockchalk mvrnorm
#' @import resample
#' @import doParallel
#' @import parallel
#' @import foreach
#' @importFrom stats cor quantile
#' @import CVXR
RPM_AB_delta <- function(x,y,trt,H_vars,target_mean,
                         target_sd,num_sim,
                         H_add_intercept=TRUE,delta,
                         cluster=1, set_seed=111){

  if (length(target_mean) > NCOL(x)){
    stop("Error: number of target moments must be smaller than the
         number of source covariates")
  }

  if (length(target_mean) != length(H_vars)){
    stop("Error: number of target moments must be equal to the number of
         specified covariates in the source")
  }

  nx_s <- nrow(x)
  ## Parallel processing
  cl <- makeCluster(cluster)
  registerDoParallel(cl)
  if(!is.null(set_seed)){set.seed(set_seed, kind = "L'Ecuyer-CMRG")}

  perturb_CI <- foreach(1:num_sim, .combine = rbind,
                        .packages = c("dplyr","rockchalk","CVXR"),
                        .export=c("ebal_wts","weighted_ATE"),
                        .errorhandling = 'remove') %dorng% {

                          target_moments <- rockchalk::mvrnorm(1,target_mean,diag(target_sd)%*%stats::cor(x[,H_vars])%*%diag(target_sd))

                          id_source <- sample(nx_s,nx_s,replace = TRUE)
                          x_s <- x[id_source,]
                          trts <- trt[id_source]
                          y_s <- y[id_source]

                          #wts_gen <- ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,delta)$w
                          wts_gen <- tryCatch(
                            ebal_wts(x=x_s, trt=trts,H_vars, target_moments,H_add_intercept,delta)$w,
                            error = function(e) {
                              return(rep(NA,nx_s))
                            }
                          )
                          ## R function optim() would sometimes give all wts as 0 but declare convergence
                          #if (round(sum(wts_gen))!=nx_s*2| all(is.na(wts_gen))){
                          #  wts_gen=rep(NA,nx_s)
                          #}


                          # while(all(is.na(wts_gen))){
                          #   constant = constant + 0.1
                          #   wts_gen <- approx_bal(x_s, trts, target_moments,H_add_intercept = TRUE,delta=delta*constant)$w
                          # }

                          approx_ATE <- weighted_ATE(y_s,trts,wts_gen)


                        }
  on.exit(stopCluster(cl))
  perturb_CI <- data.frame(perturb_CI)
  ## get the CIs and see how many times the ATET is in the CI
  if (nrow(perturb_CI) ==0 | sum(!is.na(perturb_CI[,1]))==0){
    stop("Error: there is no solution for this specified tolerence delta")
  } else{

    n_success <- sum(!is.na(perturb_CI[,1]))
    mean_ATE <- mean(perturb_CI[,1],na.rm=TRUE)
    lb_ATE <- stats::quantile(perturb_CI[,1],0.025,na.rm=TRUE)
    ub_ATE <- stats::quantile(perturb_CI[,1],0.975,na.rm=TRUE)

    list(mean_ATE=mean_ATE,lb_ATE=lb_ATE, ub_ATE=ub_ATE,n_success=n_success)
  }


}

