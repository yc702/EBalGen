#' @title Stabilize weights and get ATE
#' @param y A vector of source response values.
#' @param trt A vector of 0, 1 or FALSE/TRUE of treatment assignment for the source sample.
#' @param wts A vector of source weights.
#'
#' @return ATE numeric value generalized from source to the target sample.
#' @examples
#' library(EBalGen)
#' set.seed(1)
#' n = 100
#' p = 5
#' x = runif(n * p)
#' x = matrix(4 * x - 2, n, p)
#' y = rnorm(n*p)
#' trt = rbinom(n,1,0.5)
#' wts = plogis(rnorm(n))
#' .weighted_ATE(y, trt,wts)
#' @noRd
#' @import dplyr
#' @importFrom dplyr %>%
.weighted_ATE <- function(y,trt,wts){
  wts_gen <- wts %>% dplyr::as_tibble() %>%
    dplyr::group_by(trt) %>%
    dplyr::mutate_at(vars(-group_cols()), function(z) z / sum(z)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-trt)

  colSums((2 * trt - 1) * y * wts_gen)

}



#' @title Entropy balancing ATE for causal generalization
#' @description Compute the exact and approximate Average Treatment Effect (ATE)
#' using entropy balancing weights proposed in the paper.
#'
#' If the tolerance margin \code{delta} is all 0, it computes the exact balancing ATE.
#' Otherwise, it computes the approximate balancing ATE. If exact balancing
#' does not yield a feasible solution, the standard deviation of \code{x} is
#' used as \code{delta} , which convert exact into approximate balancing.
#'
#' If the specified \code{delta} does not yield a feasible solution, for approximate balancing,
#' the constant is increased (starting from 1) by 1 times delta until a solution is found.
#' For exact balancing that later uses the standard deviation as \code{delta},
#' the constant is increased (starting from 0) by 0.1 times delta until a solution is achieved.
#' @param x A data matrix for the source sample. Each column represents
#' source sample covariate and each row represents an observation.
#' @param y A vector of the source sample response values.
#' @param trt A vector of 0, 1 or TRUE/FALSE of treatment assignment for the source sample.
#' @param H_vars A vector of numbers indexing which covariate in \code{x}
#' need to be balanced between source and target samples.
#' @param target_moments A vector of first moments of the target sample covariates
#' that needs to be balanced between source and target.
#' @param H_add_intercept `logical` whether to include 1 as intercept in
#' H covariates, default as TRUE.
#' @param delta A vector specifying the approximate balancing tolerance margin.
#' The vector has a total length of H+H+G, where H represents the number of
#' covariates balanced between the source (treatment and control) and the
#' target moments, and G represents the covariates balanced solely between
#' the source treatment and control groups. If exact balancing, delta are all zeros.
#'
#' @return A list containing:
#' \describe{
#'   \item{ate_est}{ATE for causal generalization.}
#'  \item{constant}{The final constant used for the approximate balancing
#'  tolerance margin if no feasible solution is achieved with the specified `delta`.
#'  If the specified `delta` results in a feasible solution, the constant remains 0.
#'  Otherwise, the constant is incrementally increased, multiplying `delta`
#'  until a feasible solution is found.}
#'  }
#' @details
#' `ebal_ATE` computes causal generalized ATE from source to the target sample by
#' estimating entropy balancing weights. See Chen, Chen, Yu (2023) method details
#' for the exact balancing set up.
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
#' target_moments = c(0,0,0)
#'
#' ## Exact balancing ATE
#' delta = numeric(8)
#' library(dplyr)
#' ebal_ATE(x,y,trt,H_vars, target_moments = target_moments, H_add_intercept = TRUE,delta)
#'
#' ## Approximate balancing ATE requires installation of MOSEK for optimization.
#' delta = numeric(8)+0.1
#' if (requireNamespace("CVXR", quietly = TRUE)) {
#'   library(rlang)
#'   library(dplyr)
#'   if (!("MOSEK" %in% CVXR::installed_solvers())) {
#'       rlang::abort("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")}
#'   ebal_ATE(x,y,trt,H_vars, target_moments = target_moments, H_add_intercept = TRUE,delta)
#'
#' }
#' @rdname ebal_ATE
#' @export
#' @import dplyr
#' @import resample
ebal_ATE <- function(x,y,trt,H_vars,
                     target_moments,
                     H_add_intercept=TRUE,delta){

  if (length(target_moments) > NCOL(x)){
    stop("Error: number of target moments must be smaller than the
         number of source covariates")
  }

  if (length(target_moments) != length(H_vars)){
    stop("Error: number of target moments must be equal to the number of
         specified covariates in the source")
  }

  nx_s <- nrow(x)
  sd_s = resample::colStdevs(x)

  ## Get entropy balancing weights with the specified delta
  wts_gen <- tryCatch(
    ebal_wts(x, trt, H_vars,target_moments,H_add_intercept = TRUE,delta)$w,
    error = function(e) {
      return(rep(NA,nx_s))
    }
  )
  #wts_gen <- ebal_wts(x, trt, H_vars,target_moments,H_add_intercept = TRUE,delta)$w

  ## If there is no feasible solution in exact balancing, use SD of x as delta
  if (all(delta==0)){
    delta = c(sd_s[H_vars],sd_s[H_vars],sd_s[-H_vars])
    # wts_gen <- approx_bal(xs, trts, colMeans(xt)[H_vars],H_add_intercept = TRUE,delta=delta*0)$w
    ## R function optim() would sometimes give all wts as 0 but declare convergence
    #if (round(sum(wts_gen))!=nx_s*2| all(is.na(wts_gen))){
    #  wts_gen=rep(NA,nx_s)
    #}
    constant = 0
    while(all(is.na(wts_gen))){
      ## increment by 0.1 of the constant \times SD of x as delta
      constant = constant + 0.1
      #wts_gen <- ebal_wts(x, trt, H_vars, target_moments,H_add_intercept = TRUE,delta=delta*constant)$w
      wts_gen <- tryCatch(
        ebal_wts(x, trt, H_vars, target_moments,H_add_intercept = TRUE,delta=delta*constant)$w,
        error = function(e) {
          return(rep(NA,nx_s))
        }
      )
    }

  } else {
    # wts_gen <- approx_bal(xs, trts, colMeans(xt)[H_vars],H_add_intercept = TRUE,delta=delta*0)$w
    ## R function optim() would sometimes give all wts as 0 but declare convergence
    #if (round(sum(wts_gen))!=nx_s*2| all(is.na(wts_gen))){
    #  wts_gen=rep(NA,nx_s)
    #}

    constant = 1
    while(all(is.na(wts_gen))){
      ## increment by 1 of the constant \times delta
      constant = constant + 1
      #wts_gen <- ebal_wts(x, trt, H_vars, target_moments,H_add_intercept = TRUE,delta=delta*constant)$w
      wts_gen <- tryCatch(
        ebal_wts(x, trt, H_vars, target_moments,H_add_intercept = TRUE,delta=delta*constant)$w,
        error = function(e) {
          return(rep(NA,nx_s))
        }
      )
    }
    constant=ifelse(constant==1,0,constant)
  }
  ate_est <- .weighted_ATE(y=y,trt=trt,wts=wts_gen)

  list(ATE = ate_est, constant = constant)
}
