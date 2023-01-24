#' Sure independence screening followed by lasso
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}.
#' @param family Either a character string representing one of the built-in exponential families, including "gaussian", "poisson", and "binomial". Default is "gaussian".
#' @param num_keep Number of variables to keep in the screening phase
#' @param ... other arguments to be passed to the \code{glmnet} calls, such as \code{alpha} or \code{penalty.factor}
#'
#' @return An object of S3 class "\code{sis_lasso}".
#'  \describe{
#'   \item{\code{n}}{The sample size.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{a0}}{Intercept value.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'  }
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 200
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- sis_lasso(x = x, y = y)
#'
#' @import glmnet
#' @export
sis_lasso <- function(x, y, family = "gaussian", num_keep = NULL, ...){

  n <- nrow(x)
  p <- ncol(x)
  stopifnot(n == length(y))

  if(is.null(num_keep))
    num_keep <- ceiling(n / log(n))

  # step 1: screen interactions
  offset <- matrix(0, nrow = n, ncol = 1)
  if(inherits(x, "sparseMatrix")){
    idx <- screen_sparse_cpp(x = x, y = y, nlev = 1, offset = offset,  num_keep = num_keep, square = FALSE, main_effect = FALSE, family = family)
  }else{
    idx <- screen_cpp(x = x, y = y, nlev = 1, offset = offset, num_keep = num_keep, square = FALSE, main_effect = FALSE,  family = family)
  }

  # step 2: fit cv.glment with all main effects and selected interactions
  x <- myscale(x)
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")

  # construct design matrix of pairwise interactions
  xx <- myscale(x[, idx[, 1]] * x[, idx[, 2]])
  col_mean <- c(col_mean, attr(x = xx, which = "scaled:center"))
  col_sd <- c(col_sd, attr(x = xx, which = "scaled:scale"))

  design <- cbind(x, xx)

  fit <- glmnet::cv.glmnet(x = design, y = y, family = family,
                           intercept = FALSE,
                           standardize = FALSE, ...)

  coef <- fit$glmnet.fit$beta[, which.min(fit$cvm)]
  # scale estimates back to the original scale of design matrix
  coef <- coef / col_sd
  a0 <- as.numeric(fit$glmnet.fit$a0[which.min(fit$cvm)] - crossprod(col_mean, coef))

  idx_all <- rbind(cbind(rep(0, p), seq(p)), idx[, 1:2])
  compact <- cbind(idx_all[which(coef != 0), , drop = FALSE], coef[coef != 0])
  rownames(compact) <- NULL
  colnames(compact) <- c("index_1", "index_2", "coefficient")

  # finally return the best lambda
  out <- list(n = n,
              p = p,
              a0 = a0,
              compact = compact,
              call = match.call())
  class(out) <- "other"
  return(out)
}
