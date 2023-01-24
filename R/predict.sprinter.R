#' Calculate prediction from a \code{sprinter} object.
#'
#' @param object a fitted \code{sprinter} object.
#' @param newdata a design matrix of all the \code{p} main effects of some new observations of which predictions are to be made.
#' @param ... additional argument (not used here, only for S3 generic/method consistency)
#' @return The prediction of \code{newdata} by the sprinter fit \code{object}. The prediction is on the scale of the linear predictors, i.e., \code{type = "link"} as in \code{predict.glm} or \code{predict.glmnet}. For example, for "gaussian" family it gives the fitted response values, and for "binomial" family the predictions are of log-odds (logit).
#'
#'  The output is a list consisting of \code{length(object$lambda1)} elements. For the "gaussian", "poisson", and "binomial" families, each element is a matrix of dimensions \code{n}-by-\code{nrow(object$lambda3)}. For "ordinal" family, each element is a list with \code{length(object$lambda3)} components, each of which is an \code{n}-by-\code{[object$nlev-1]} matrix.
#'
#' @examples
#' n <- 100
#' p <- 200
#' x <- matrix(rnorm(n * p), n, p)
#' y <- x[, 1] + 2 * x[, 2] - 3 * x[, 1] * x[, 2] + rnorm(n)
#' mod <- sprinter(x = x, y = y)
#' fitted <- predict(mod, newdata = x)
#'
#' @export
predict.sprinter <- function(object, newdata, ...) {
  
  # input check
  stopifnot(ncol(newdata) == object$p)
  
  n <- nrow(newdata)
  nlev <- object$nlev
  nlam1 <- length(object$lambda1)
  nlam3 <- length(object$lambda3[,1])
  
  # need to standardize the main effects to construct interactions
  xm <- myscale(newdata)
  
  if(object$square){
    x_step1 <- cbind(newdata, xm^2)
  }else{
    x_step1 <- newdata
  }
  
  # output is a list consisting of nlam1 elements.
  out <- vector("list", length = nlam1)
  
  # For the "gaussian", "poisson", and "binomial" families, each element is a matrix of dimensions n-by-nlam3
  # For "ordinal" family, each element is a list with lambda3 components, each of which is an n-by-[nlev-1] matrix.
  # we add the prediction in Step1 (for each lambda1 value) and in Step3 (for a path of lambda3 values)
  if(object$family == "ordinal"){
    for(k in seq(nlam1)){
      mu_step1 <- matrix(object$step1$a0[, k], nrow = n, ncol = nlev-1, byrow = TRUE) + x_step1 %*% (matrix(object$step1$beta[, k], ncol = 1)[ , rep(1, nlev-1), drop=FALSE])
      
      idx <- object$step2[[k]]
      xint <- xm[, idx[, 1]] * xm[, idx[, 2]]
      design <- cbind(x_step1, xint)
      
      out[[k]] <- vector("list", length = nlam3)
      for(j in seq(nlam3)){
        mu_step3 <- matrix(object$step3[[k]]$a0[, j], nrow = n, ncol = nlev-1, byrow = TRUE) + design %*% (matrix(object$step3[[k]]$coef[, j], ncol = 1)[ , rep(1, nlev-1), drop=FALSE])
        mu <- mu_step1 + mu_step3
        out[[k]][[j]]  <- mu
      }
    }
  }else{
    for(k in seq(nlam1)){
      mu_step1 <- as.numeric(object$step1$a0[k] + x_step1 %*% object$step1$beta[, k])
      
      idx <- object$step2[[k]]
      xint <- xm[, idx[, 1]] * xm[, idx[, 2]]
      design <- cbind(x_step1, xint)
      
      mu_step3 <- t(object$step3[[k]]$a0 + t(as.matrix(design %*% object$step3[[k]]$coef)))
      out[[k]] <- mu_step1 + mu_step3
      colnames(out[[k]]) <- NULL
    }
  }
  
  return(out)
  
}
