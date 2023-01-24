#' Running sprinter with cross-validation
#'
#' The main cross-validation function to select the best sprinter fit for a path of tuning parameters.
#'
#' @param x An \code{n} by \code{p} design matrix of main effects. Each row is an observation of \code{p} main effects.
#' @param y A response vector of size \code{n}. \code{y} is a vector of continuous data in "gaussian" family; a vector of positive count data in "poisson" family; a vector of 0-1 binary data in "binomial" family. Notice that \code{y} should be a factor when using the "ordinal" family for ordinal response.
#' @param family Either a character string representing one of the built-in exponential families, including "gaussian", "poisson", "binomial", "ordinal". Default is "gaussian".
#' @param nlev The number of levels for response variable. Set to 1 for "gaussian" and "poisson" family.
#' @param type.measure type of validation error used for cross-validation. Currently only supporting 2 options. The default is \code{type.measure = "deviance"}, which uses squared-error for gaussian models, deviance for logistic, poisson and ordinal logistic regression. \code{type.measure = "auc"} is for two-class logistic regression only, and gives area under the ROC curve.
#' @param square Indicator of whether squared effects should be fitted in Step 1. Default to be FALSE.
#' @param num_keep A user specified number of candidate interactions to keep in Step 2. If \code{num_keep} is not specified (as default), it will be set to \code{round[n / log n]}.
#' @param lambda1 Tuning parameter values for Step 1. \code{lambda1} is a vector. Default to be NULL, and the program will compute its own \code{lambda1} based on \code{nlam1} and \code{lam_min_ratio}.
#' @param lambda3 Tuning parameter values for Step 3. \code{lambda3} is a matrix, where the k-th column is the list of tuning parameter in Step 3 corresponding to Step 1 using \code{lambda1[k]}. Default to be NULL, and the program will compute its own \code{lambda3} based on \code{nlam3} and \code{lam_min_ratio}.
#' @param cv_step1 Indicator of whether cross-validation of \code{lambda1} should be carried out in Step 1 before subsequent steps. Default is \code{FALSE}.
#' @param nlam1 the number of values in \code{lambda1}. If not specified, they will be all set to \code{10}.
#' @param nlam3 the number of values in each column of \code{lambda3}. If not specified, they will be all set to \code{100}.
#' @param lam_min_ratio The ratio of the smallest and the largest values in \code{lambda1} and each column of \code{lambda2}. The largest value is usually the smallest value for which all coefficients are set to zero. Default to be \code{1e-2} in the \code{n} < \code{p} setting.
#' @param nfold Number of folds in cross-validation. Default value is 5. If each fold gets too view observation, a warning is thrown and the minimal \code{nfold = 3} is used.
#' @param foldid A vector of length \code{n} representing which fold each observation belongs to. Default to be \code{NULL}, and the program will generate its own randomly.
#' @param verbose If \code{TRUE}, a progress bar shows the progress of the fitting.
#' @param ... other arguments to be passed to the \code{glmnet} calls, such as \code{alpha} or \code{penalty.factor}
#'
#' @return An object of S3 class "\code{sprinter}".
#'  \describe{
#'   \item{\code{n}}{The number of observations in the dataset.}
#'   \item{\code{p}}{The number of main effects.}
#'   \item{\code{square}}{The \code{square} parameter passed into sprinter.}
#'   \item{\code{family}}{The exponential \code{family} used in sprinter.}
#'   \item{\code{nlev}}{The number of levels of the response variable.}
#'   \item{\code{cv_step1}}{The \code{cv_step1} parameter passed into sprinter, indicating whether to perform cross-validation of \code{lambda1} in Step 1.}
#'   \item{\code{a0}}{Estimate of intercept corresponding to the CV-selected model.}
#'   \item{\code{compact}}{A compact representation of the selected variables. \code{compact} has three columns, with the first two columns representing the indices of a selected variable (main effects with first index = 0), and the last column representing the estimate of coefficients.}
#'   \item{\code{fit}}{The whole \code{glmnet} fit object.}
#'   \item{\code{num_keep}}{The value of \code{num_keep}.}
#'   \item{\code{type.measure}}{The type of validation error used for cross-validation.}
#'   \item{\code{cvm}}{The averaged estimated error on the test sets over K folds.}
#'   \item{\code{foldid}}{Fold assignment. A vector of length \code{n}.}
#'   \item{\code{i_lambda1_best}}{The index in \code{lambda1} that is chosen by CV by minimizing cvm.}
#'   \item{\code{i_lambda3_best}}{The index in \code{lambda3} that is chosen by CV by minimizing cvm.}
#'   \item{\code{lambda1_best}}{The value of \code{lambda1} that is chosen by CV by minimizing cvm.}
#'   \item{\code{lambda3_best}}{The value of \code{lambda3} that is chosen by CV by minimizing cvm.}
#'   \item{\code{call}}{Function call.}
#'  }
#' @seealso
#'   \code{\link{predict.cv.sprinter}}
#' @examples
#' set.seed(123)
#' n <- 100
#' p <- 100
#' # dense input
#' x <- matrix(rnorm(n * p), n, p)
#' mu <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5]
#'
#' # Gaussian
#' y <- mu + rnorm(n)
#' mod_gau <- cv.sprinter(x = x, y = y)
#'
#' # Binomial
#' prob = exp(mu) / (1 + exp(mu))
#' set.seed(123)
#' y = rbinom(n = n, prob = prob, size = 1)
#' mod_bin <- cv.sprinter(x = x, y = y, family = "binomial", type.measure = "deviance")
#'
#' # sparse input
#' library(Matrix)
#' x <- Matrix::Matrix(0, n, p)
#' idx <- cbind(sample(seq(n), size = 10, replace = TRUE), sample(seq(p), size = 10, replace = TRUE))
#' x[idx] <- 1
#' y <- x[, 1] - 2 * x[, 2] + 3 * x[, 1] * x[, 3] - 4 * x[, 4] * x[, 5] + rnorm(n)
#' mod <- sprinter(x = x, y = y)
#'
#' @import glmnet
#' @export
cv.sprinter <- function(x, y, family = "gaussian", nlev = 1, type.measure = "deviance", square = FALSE, num_keep = NULL, lambda1 = NULL, lambda3 = NULL, cv_step1 = FALSE, nlam1 = 10, nlam3 = 100, lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), nfold = 5, foldid = NULL, verbose = FALSE, ...){
  
  n <- nrow(x)
  p <- ncol(x)
  q <- ifelse(square, p * (p - 1) / 2, p * (p - 1) / 2 + p)
  
  if(family == "binomial"){
    nlev <- 2
  } else if(family == "ordinal"){
    if (!is.factor(y)) stop("y should be a factor in ordinal family.")
    nlev <- nlevels(y)
  }
  
  if (is.null(num_keep)) {
    num_keep <- ceiling(n / log(n))    
  } else {
    stopifnot(num_keep > 0 & num_keep <= q)
  }
  
  # first fit the sprinter using all data with all lambdas
  fit <- sprinter(x = x, y = y, family = family, nlev = nlev, 
                  square = square, num_keep = num_keep,
                  lambda1 = lambda1, lambda3 = lambda3,
                  cv_step1 = cv_step1,
                  nlam1 = nlam1, nlam3 = nlam3, 
                  lam_min_ratio = lam_min_ratio, ...)
  
  nlam1 <- length(fit$lambda1)
  nlam3 <- nrow(fit$lambda3)
  
  # use cross-validation to select the best lambda
  if (is.null(foldid)){
    # foldid is a vector of values between 1 and nfold
    # identifying what fold each observation is in.
    # If supplied, nfold can be missing.
    foldid <- sample(seq(nfold), size = n, replace = TRUE)
  }
  
  # estimated error of lasso estimate of coef
  err <- vector("list", nfold)
  err_mat <- matrix(0, nrow = nlam1, ncol = nlam3)
  for (i in seq(nfold)){
    if(verbose)
      cat(paste("cv fold ", i, ":"), fill = TRUE)
    # train on all but i-th fold
    id_tr <- (foldid != i)
    id_te <- (foldid == i)
    # training/testing data set
    # standardize the training data
    x_tr <- myscale(x[id_tr, ])
    # use the scaling for training data
    x_te <- myscale(x[id_te, ],
                    center = attr(x = x_tr, which = "scaled:center"),
                    scale = attr(x = x_tr, which = "scaled:scale"))
    y_tr <- y[id_tr]
    y_te <- y[id_te]
    
    # get the fit using lasso on training data
    # and fit/obj on the test data
    fit_tr <- sprinter(x = x_tr, y = y_tr,
                       family = fit$family,
                       nlev = fit$nlev,      
                       square = fit$square,
                       num_keep = fit$num_keep,
                       lambda1 = fit$lambda1,
                       lambda3 = fit$lambda3, ...)
    
    err[[i]] <- matrix(NA, nrow = nlam1, ncol = nlam3)
    
    pred_te <- predict.sprinter(object = fit_tr, newdata = x_te)
    
    for(k in seq(nlam1)){
      if(family == "ordinal"){
        for(j in seq(nlam3)){
          err[[i]][k, j] <- compute_error(y = y_te, mu = pred_te[[k]][[j]], family = family, type.measure = type.measure) 
        }
      }else{
        if(ncol(pred_te[[k]]) < nlam3){
          # this only happens when glmnet returns a numerical warning as follows:
          # Fortran error code; Convergence for i-th lambda value not reached after maxit=100000 iterations; solutions for larger lambdas returned
          # when this happens, not all nlam3 lambda values in Step3 is returned
          # we just copy the prediction from the last viable lambda values to fill in the empty slots
          last_col <- ncol(pred_te[[k]])
          ncol_to_add <- nlam3 - last_col
          pred_te[[k]] <- cbind(pred_te[[k]], matrix(rep(pred_te[[k]][, last_col], ncol_to_add), ncol = ncol_to_add))
        }
        err[[i]][k, ] <- compute_error(y = y_te, mu = pred_te[[k]],
                                       family = family,
                                       type.measure = type.measure)
      }
    }
    err_mat <- err_mat + err[[i]]
  }
  
  # extract information from CV
  # the mean cross-validated error, a nlam1-by-nlam3 matrix
  err_mat <- err_mat / nfold
  
  # the index of best lambda
  # if there is ties
  # we select the one with smaller lambda1 value and larger lambda 3 value
  # i.e., we select the model with more main effects and fewer interactions
  ibest <- tail(which(err_mat == min(err_mat), arr.ind = TRUE), 1)
  colnames(ibest) <- NULL
  lambda1_best <- ibest[1]
  lambda3_best <- ibest[2]
  
  # get a compact output
  idx_i <- fit$step2[[lambda1_best]]
  
  if(square){
    idx <- rbind(cbind(rep(0, p), seq(p), rep(0, p)),
                 cbind(seq(p), seq(p), rep(0, p)),
                 idx_i)
  }else{
    idx <- rbind(cbind(rep(0, p), seq(p), rep(0, p)), idx_i)
  }
  
  # beta is the estimate from Step 1
  beta <- fit$step1$beta[, lambda1_best]
  # coef is the estimate from Step 3
  coef <- fit$step3[[lambda1_best]]$coef[, lambda3_best]
  # we now combine coef from Step 3 and beta from Step 1
  coef[1:length(beta)] <- coef[1:length(beta)] + beta
  
  compact <- cbind(idx[which(coef != 0), 1:2, drop = FALSE], coef[coef != 0])
  colnames(compact) <- c("index_1", "index_2", "coefficient")
  
  # combine intercept from Step1 and Step3
  if(family == "ordinal"){
    a0 <- matrix(fit$step1$a0[, lambda1_best] +  fit$step3[[lambda1_best]]$a0[, lambda3_best], ncol = 1)
    rownames(a0) <- paste0("Intercept:", 1:(nlev-1))
  }else{
    a0 <- fit$step1$a0[lambda1_best] +  fit$step3[[lambda1_best]]$a0[lambda3_best]
  }
  
  # finally return the best lambda
  out <- list(n = n,
              p = p,
              square = square,
              family = family,
              nlev = nlev,
              cv_step1 = cv_step1,
              a0 = a0,
              compact = compact,
              fit = fit,
              num_keep = num_keep,
              type.measure = type.measure,
              cvm = err_mat,
              foldid = foldid,
              i_lambda1_best = lambda1_best,
              i_lambda3_best = lambda3_best,
              lambda1_best = fit$lambda1[lambda1_best],
              lambda3_best = fit$lambda3[lambda3_best, lambda1_best],
              call = match.call())
  class(out) <- "cv.sprinter"
  return(out)
}
