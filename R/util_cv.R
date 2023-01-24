#' Compute the error used in each fold of cross-validation
#'
#' @param y A response vector of size \code{n}. Usually used as the response in the held-out set.
#' @param mu A predicted linear predictor values. Usually used as the linear predictions for the held-out set over a path of tuning parameter values.
#' @param family Either a character string representing one of the built-in exponential families. Default is "gaussian".
#' @param type.measure type of validation error used for cross-validation. Currently only supporting 2 options. The default is \code{type.measure = "deviance"}, which uses squared-error for gaussian models, deviance for logistic, poisson and ordinal logistic regression. \code{type.measure = "auc"} is for two-class logistic regression only, and gives area under the ROC curve.
#'
compute_error <- function(y, mu, family, type.measure){
  if(family == "gaussian")
    return(as.numeric(colMeans((y - mu)^2)))
  else if(family == "binomial"){
    if(type.measure == "deviance"){
      return(compute_deviance_binomial(y, mu))
    }
    else if(type.measure == "auc"){
      return(compute_auc(y, mu))
    }
  }
  else if(family == "poisson"){
    if(type.measure == "deviance"){
      return(compute_deviance_poisson(y, mu))
    }
  }
  else if(family == "ordinal"){
    if(type.measure == "deviance"){
      return(compute_deviance_ordinal(y, mu))
    }
  }
}

#' Compute deviance of the logistic regression model
#'
#' @param y A response vector of size \code{n}. Usually used as the response in the held-out set.
#' @param mu A predicted linear predictor values. Usually used as the linear predictions for the held-out set over a path of tuning parameter values, i.e., \code{mu} is a matrix of \code{n} rows and the number of tuning parameters columns.
#'
#' @return The output is deviance of the logistic regression model for each value of tuning parameters used to compute each column of mu, i.e., \code{mu} is a matrix of \code{n} rows and the number of tuning parameters columns.
#'
#' For logistic regression the deviance formula is minus twice the log-likelihood, which is consistent with glmnet
#'
compute_deviance_binomial <- function(y, mu){
  # another definition -2 (X beta)^T y - 2 sum_{i = 1}^n log (1 + exp(X_i^T beta))
  #  dev <- -2 * as.numeric(crossprod(mu, y) + colSums(log(1 + exp(mu))))

  pp <- 1 / (1 + exp(-mu))
  pp[pp == 0] <- .Machine$double.eps * 5
  pp[pp == 1] <- 1 - .Machine$double.eps * 5
  dev <- -2 * as.numeric(colSums(((y == 1) * log(pp) + (y == 0) * log(1 - pp))))
  return(dev)
}

#' Compute the NEGATIVE value of the area under the ROC curve (AUC) for the binary logistic model
#'
#' @param y A response vector of size \code{n}. Usually used as the response in the held-out set.
#' @param mu A predicted linear predictor values. Usually used as the linear predictions for the held-out set over a path of tuning parameter values, i.e., \code{mu} is a matrix of \code{n} rows and the number of tuning parameters columns.
#'
#' @return The output is the negative value of AUC for each value of tuning parameters used to compute each column of mu. We return the negative value to be consistent with the cross-validation procedure which selects the tuning parameter that minimizes CV-estimate of test error.
#'
#' Estimate of AUC could be inaccurate if \code{n} is small.
#'
compute_auc <- function(y, mu){

  nlam <- ncol(mu)
  # sort fitted mu for each lambda
  # note that the fitted probability is monotonically non-decreasing with the value of mu
  # so sorting the probability non-decreasingly is equivalent to sorting mu non-decreasingly
  musort <- apply(mu, 2, sort, decreasing = TRUE, index.return = TRUE)
  idx_musort <- matrix(unlist(lapply(musort, `[[`, 2)), ncol = nlam)

  # sort true y by using the order of fitted mu
  roc_y <- matrix(y[idx_musort], ncol = nlam)

  # apply Trapezoidal rule to calculate auc (area under roc curve)
  stack_x <- apply(roc_y == 0, 2, cumsum)/apply(roc_y == 0, 2, sum)
  stack_y <- apply(roc_y == 1, 2, cumsum)/apply(roc_y == 1, 2, sum)
  diff_x <- stack_x[-1, ] - stack_x[-nrow(stack_x), ]
  mean_y <- (stack_y[-1, ] + stack_y[-nrow(stack_y), ])/2
  auc <- apply(diff_x * mean_y, 2, sum)

  return(-auc)
}

#' Compute deviance of the Poisson regression model
#'
#' @param y A response vector of size \code{n}. Usually used as the response in the held-out set.
#' @param mu A predicted linear predictor values. Usually used as the linear predictions for the held-out set over a path of tuning parameter values, i.e., \code{mu} is a matrix of \code{n} rows and the number of tuning parameters columns.
#'
#' @return The output is deviance of the Poisson regression model for each value of tuning parameters used to compute each column of mu.
compute_deviance_poisson <- function(y, mu){

  pp <-  exp(mu)
  pp[pp == 0] <- .Machine$double.eps * 5

  # poisson deviance: 2*colSums(y * log(y / pp) + pp - y) where y > 0
  # deviance can be omitted some terms if use for comparison
  dev <- 2 * as.numeric(colSums(- y * log(pp) + pp))
  return(dev)
}

#' Compute deviance of the ordinal logistic regression model
#'
#' @param y A response vector of size \code{n}. Usually used as the response in the held-out set.
#' @param mu A predicted linear predictor values in the form of an \code{n}-by-\code{[object$nlev-1]} matrix. Usually used as the linear predictions for the held-out set over a path of tuning parameter values.
#'
#' @return The output is deviance of the ordinal logistic regression model for each value of tuning parameters.
compute_deviance_ordinal <- function(y, mu){

  # transform y to yMat with dimension n * nlev
  nlev <- nlevels(y)
  nobs <- length(y)
  yMat <- matrix(0, nrow=nobs, ncol=nlev)
  yInt <- as.integer(y)
  yMat[cbind(1:nobs, yInt)] <- 1

  # ordinal logistic regression deviance
  delta = t(apply(mu, 1, invLogit))
  pMat = delta - cbind(0, delta[,-ncol(delta)])
  pMat <- cbind(pMat, 1-rowSums(pMat))
  pMat[pMat == 0] <- .Machine$double.eps * 5
  pMat[pMat == 1] <- 1 - .Machine$double.eps * 5
  pMat <- pMat / rowSums(pMat)
  dev <- -2 * sum(yMat * log(pMat))

  return(dev)
}


