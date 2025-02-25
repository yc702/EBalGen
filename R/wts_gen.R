#' @title Entropy balancing weights for causal generalization
#' @description Compute the exact and approximate entropy balancing weights
#' proposed in the paper. If the tolerance margin \code{delta} is all 0,
#' it computes exact balancing weights. Otherwise, it computes approximate
#' balancing weights.
#' @param x A data matrix for the source sample. Each column represents
#' source sample covariate and each row represents an observation.
#' @param trt A vector of 0, 1 or FALSE/TRUE of treatment assignment for the source sample.
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
#' @return A list containing:
#' \item{w}{A vector of entropy balancing weights for causal generalization.}
#' \item{theta}{The dual parameters estimated from the optimization process.}
#' @examples
#' library(EBalGen)
#' set.seed(1)
#' n = 100
#' p = 5
#' x = runif(n * p)
#' x = matrix(4 * x - 2, n, p)
#' trt = rbinom(n,1,0.5)
#' H_vars = c(1,2,3)
#' target_moments = c(0,0,0)
#' delta = numeric(8)
#'
#' ## Exact balancing weights
#' ebal_wts(x, trt,H_vars, target_moments = target_moments, H_add_intercept = TRUE,delta)
#'
#' ## Approximate balancing weights requires installation of MOSEK for optimization.
#' delta = numeric(8)+0.1
#' if (requireNamespace("CVXR", quietly = TRUE)) {
#'   library(rlang)
#'   if (!("MOSEK" %in% CVXR::installed_solvers())) {
#'       rlang::abort("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")}
#'   ebal_wts(x, trt,H_vars, target_moments = target_moments, H_add_intercept = TRUE,delta)
#'
#' }
#' @rdname ebal_wts
#' @export
ebal_wts <- function(x, trt,H_vars,
                     target_moments = NULL,
                     H_add_intercept = TRUE,
                     delta) {
  if (length(target_moments) > NCOL(x)){
    stop("Error: number of target moments must be smaller than the
         number of source covariates")
  }

  if (length(target_moments) != length(H_vars)){
    stop("Error: number of target moments must be equal to the number of
         specified covariates in the source")
  }


  # If H_add_intercept is TRUE, include 1 in H, so that the sums of the weights
  # in the treated and control group will be kept as 1 in the computation
  if (H_add_intercept)
  {

    delta <- c(0,delta[1:length(target_moments)],0,
               delta[(length(target_moments)+1):length(delta)])
    x <- cbind(1, x)
    target_moments <- c(1, target_moments)
    H_vars <- c(1,H_vars+1)

  }

  if (length(H_vars)==0){
    H1 <- x[trt == 1, numeric(0), drop = FALSE]
    H0 <- x[trt != 1, numeric(0), drop = FALSE]
    G1 <- x[trt == 1, drop = FALSE]
    G0 <- x[trt != 1, drop = FALSE]
  } else {
    H1 <- x[trt == 1, H_vars, drop = FALSE]
    H0 <- x[trt != 1, H_vars, drop = FALSE]
    G1 <- x[trt == 1, -H_vars, drop = FALSE]
    G0 <- x[trt != 1, -H_vars, drop = FALSE]
  }


  n <- NROW(x)
  K_h <- length(target_moments)
  K_g <- NCOL(x) - K_h

  ## Primal optimization problem
  primal_dual_args <- function(theta) {
    lambda1 <- theta[(0:K_h)[-1]]
    lambda0 <- theta[K_h + (0:K_h)[-1]]
    gamma <- theta[2 * K_h + (0:K_g)[-1]]
    if (K_g==0){
      list(lambda1 = lambda1, lambda0 = lambda0, gamma = gamma,
           w1 = drop(exp(H1 %*% lambda1  - 1)),
           w0 = drop(exp(H0 %*% lambda0  - 1)))
    } else if (K_h==0){
      list(lambda1 = lambda1, lambda0 = lambda0, gamma = gamma,
           w1 = drop(exp( G1 %*% gamma - 1)),
           w0 = drop(exp(G0 %*% gamma - 1)))
    } else {
      list(lambda1 = lambda1, lambda0 = lambda0, gamma = gamma,
           w1 = drop(exp(H1 %*% lambda1 + G1 %*% gamma - 1)),
           w0 = drop(exp(H0 %*% lambda0 - G0 %*% gamma - 1)))
    }

  }

  ## Dual optimization problem
  dual_fn <- function(theta) {
    args <- primal_dual_args(theta)
    (sum(args$w1) + sum(args$w0)) / n -
      sum((args$lambda1 + args$lambda0) * target_moments)
  }

  ## When delta are all zero, perform exact balancing
  if (all(delta==0)){
    ## Gradient of the dual problem
    dual_grad <- function(theta) {
      args <- primal_dual_args(theta)
      c(colSums(args$w1 * H1) / n - target_moments,
        colSums(args$w0 * H0) / n - target_moments,
        colSums(args$w1 * G1) / n - colSums(args$w0 * G0) / n)
    }

    # Optimization is done with the built-in optim function
    opt <- optim(rep(0, 2 * K_h + K_g), fn = dual_fn, gr = dual_grad, method = "BFGS",
                 control = list(trace = 0,
                                reltol = 1e-10,
                                maxit = 1e4))

    # Output weights and dual parameters
    args <- primal_dual_args(opt$par)
    w <- numeric(n)
    w[trt == 1] <- args$w1
    w[trt != 1] <- args$w0
    if (round(sum(w,na.rm = TRUE))!=n*2 | opt$convergence!=0){
      stop("Error: there is no solution for this exact balancing")
    } else{

      list(w = w,
           theta = opt$par)
    }



  } else{
    if (!requireNamespace("CVXR", quietly = TRUE)) {
      # rlang::abort(
      #   "Package 'CVXR' must be installed when using approximate balancing.")
      stop("Package 'CVXR' must be installed when using approximate balancing.")
    }

    if (!("MOSEK" %in% CVXR::installed_solvers())) {
      # rlang::abort(
      #   "MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")
      stop("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")

    }
    ## When delta are not all zero, perform approximate balancing
    p = 2 * K_h + K_g
    theta <- CVXR::Variable(p)
    loss <- dual_fn(theta)

    ## The dual problem is the obj function, which is exact loss plus L1 penalty.
    obj <- loss + sum(CVXR::multiply(delta,abs(theta)))
    prob <- CVXR::Problem(CVXR::Minimize(obj))

    ## Use MOSEK to solve the optimization as the problem contains exp.

    # result <- CVXR::solve(prob,solver="MOSEK")
    args_theta <- tryCatch(
      (CVXR::solve(prob,solver="MOSEK"))$getValue(theta),
      error = function(e) {
        return(NA)
      }
    )
    args <- primal_dual_args(args_theta)
    wts_gen <- numeric(n)
    wts_gen[trt == 1] <- args$w1
    wts_gen[trt != 1] <- args$w0
    if (round(sum(wts_gen,na.rm = TRUE))!=n*2 | all(is.na(wts_gen))){
      stop("Error: there is no solution for this approximate balancing")
    } else{
      # Output weights and dual parameters
      list(w=wts_gen,
           theta = args_theta)
    }

  }

}


#' @title Simple entropy balancing weights for causal generalization
#' @description Compute the simple exact and approximate entropy balancing weights.
#' The resulting weights calibrate the whole x sample to the target moments (not
#' distinguishing treated and control), which is equivalent to an exponential
#' tilting calibration. If the tolerance margin \code{delta} is all 0,
#' it computes exact balancing weights. Otherwise, it computes approximate balancing weights.
#' @param x A data matrix for the source sample. Each column represents
#' source sample covariate and each row represents an observation.
#' @param target_moments A vector of first moments of the target sample covariates
#' that needs to be balanced between source and target.
#' @param H_add_intercept `logical` whether to include 1 as intercept in
#' H covariates, default as TRUE.
#' @param delta A vector specifying the approximate balancing tolerance margin.
#' The vector has a total length of \code{2*ncol(x)}. If exact balancing, delta are all zeros.
#' @return A list containing:
#' \item{w}{A vector of entropy balancing weights for causal generalization.}
#' \item{theta}{The dual parameters estimated from the optimization process.}
#' @examples
#' library(EBalGen)
#' set.seed(1)
#' n = 100
#' p = 5
#' x = runif(n * p)
#' x = matrix(4 * x - 2, n, p)
#' target_moments = c(0,0,0,0,0)
#' delta = numeric(5)
#' ebal_wts_simple(x, target_moments = target_moments, H_add_intercept = TRUE,delta)
#'
#' ## Approximate balancing weights requires installation of MOSEK for optimization.
#' delta = numeric(5)+0.1
#' if (requireNamespace("CVXR", quietly = TRUE)) {
#'   library(rlang)
#'   if (!("MOSEK" %in% CVXR::installed_solvers())) {
#'       rlang::abort("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")}
#'   ebal_wts_simple(x, target_moments = target_moments, H_add_intercept = TRUE,delta)
#'
#' }
#' @rdname ebal_wts_simple
#' @export
ebal_wts_simple <- function(x,target_moments = NULL,
                            H_add_intercept = TRUE,delta) {
  if (length(target_moments) > NCOL(x)){
    stop("Error: number of target moments must be smaller than the number of source covariates")
  }


  n <- NROW(x)

  # If H_add_intercept is TRUE, include 1 in H, so that the sums of the weights
  # in the treated and control group will be kept as 1 in the computation
  if (H_add_intercept)
  {

    delta <- c(0,delta)
    x <- cbind(1, x)
    target_moments <- c(1, target_moments)

  }
  # Function for computing the primal arguments from the dual arguments
  primal_w <- function(lambda) drop(exp(x %*% lambda - 1))

  # Set up the dual problem: functions for computing the dual objective
  dual_fn <- function(lambda) mean(primal_w(lambda)) - sum(target_moments * lambda)

  if (all(delta==0)){
    ## Gradient of the dual problem
    dual_grad <- function(lambda) colMeans(primal_w(lambda) * x) - target_moments

    # Optimization is done with the built-in optim function
    opt <- optim(rep(0, NCOL(x)), fn = dual_fn, gr = dual_grad, method = "BFGS",
                 control = list(trace = 0,
                                reltol = 1e-10,
                                maxit = 1e4))

    if (round(sum(primal_w(opt$par)))!=n | opt$convergence!=0){
      stop("Error: there is no solution for this exact balancing")
    } else{

      # Output
      list(w = primal_w(opt$par),
           theta = opt$par)
    }


  } else {
      if (!requireNamespace("CVXR", quietly = TRUE)) {
        # rlang::abort(
        #   "Package 'CVXR' must be installed when using approximate balancing.")
        stop("Package 'CVXR' must be installed when using approximate balancing.")
      }

      if (!("MOSEK" %in% CVXR::installed_solvers())) {
        # rlang::abort(
        #   "MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")
        stop("MOSEK solver is not installed. Please install MOSEK to use it with CVXR.")

      }
    p = length(target_moments)
    theta <- CVXR::Variable(p)
    loss <- dual_fn(theta)

    ## The dual problem is the obj function, which is exact loss plus L1 penalty.
    obj <- loss + sum(CVXR::multiply(delta,abs(theta)))
    prob <- CVXR::Problem(CVXR::Minimize(obj))

    ## Use MOSEK to solve the optimization as the problem contains exp.
    #result <- CVXR::solve(prob,solver="MOSEK")
    args_theta <- tryCatch(
      (CVXR::solve(prob,solver="MOSEK"))$getValue(theta),
      error = function(e) {
        return(NA)
      }
    )
    wts_gen <- primal_w(args_theta)
    if (round(sum(wts_gen,na.rm = TRUE))!=n | all(is.na(wts_gen))){
      stop("Error: there is no solution for this approximate balancing")
    } else{
      # Output weights and dual parameters
      list(w=wts_gen,
           theta = args_theta)
    }
  }

}
