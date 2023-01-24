#' Calculate prediction from a \code{cv.sprinter} object.
#'
#' @param object a fitted \code{cv.sprinter} object.
#' @param newdata a design matrix of all the \code{p} main effects of some new observations of which predictions are to be made.
#' @param ... additional argument (not used here, only for S3 generic/method consistency)
#' @return The prediction of \code{newdata} by the cv.sprinter fit \code{object}. The prediction is on the scale of the linear predictors, i.e., \code{type = "link"} as in \code{predict.glm} or \code{predict.glmnet}. For example, for "gaussian" family it gives the fitted response values, and for "binomial" family the predictions are of log-odds (logit).
#'
#'  The output is a vector of size \code{n} in the "gaussian", "poisson", and "binomial" families, or a matrix of dimensions \code{n}-by-\code{[object$nlev-1]} in "ordinal" family.
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] + 2 * x[, 2] - 3 * x[, 1] * x[, 2] + rnorm(n)
#' mod <- cv.sprinter(x = x, y = y)
#' fitted <- predict(mod, newdata = x)
#'
#' @export
predict.cv.sprinter <- function(object, newdata, ...) {
  
  # input check
  stopifnot(ncol(newdata) == object$p)
  
  n <- nrow(newdata)
  nlev <- object$nlev
  
  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)
  
  idx <- object$compact[, 1:2, drop = FALSE]
  # selected indices for main effects
  idxm <- idx[idx[, 1] == 0, 2]
  # selected index pairs for interactions
  idxi <- idx[idx[, 1] != 0, , drop = FALSE]
  
  if(nrow(idxi) == 1){
    xint <- matrix(xm[, idxi[, 1]] * xm[, idxi[, 2]], ncol = 1)
  }else{
    xint <- xm[, idxi[, 1]] * xm[, idxi[, 2]]    
  }
  
  design <- cbind(newdata[, idxm], xint)
  
  # output is a vector of size n in the "gaussian", "poisson", and "binomial" families
  # or a matrix of dimensions n-by-[nlev-1] in "ordinal" family
  mu <- if(object$family == "ordinal") matrix(object$a0, nrow = n, ncol = nlev-1, byrow = TRUE) + (design %*% object$compact[, 3])[ , rep(1, nlev-1), drop=FALSE] else as.numeric(object$a0 + design %*% object$compact[, 3])
  
  return(mu)
}
