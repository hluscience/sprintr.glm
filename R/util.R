#' Self-defined scale function
#'
#' @export
myscale <- function(x, center = TRUE, scale = TRUE){
  # Treat specifically the case where columns have zero sds and sparse matrices
  if(inherits(x, "sparseMatrix")){
    # if x is either sparse or binary
    mean_vec <- rep(0, ncol(x))
    sd_vec <- rep(1, ncol(x))
    attr(x = x, which = "scaled:center") <- mean_vec
    attr(x = x, which = "scaled:scale") <- sd_vec
  }
  else{
    n <- nrow(x)
    x <- scale(x, center = center, scale = scale)
    # for those columns with 0 sds
    # add small Gaussian noise to it
    mean_vec <- attr(x = x, which = "scaled:center")
    sd_vec <- attr(x = x, which = "scaled:scale")
    ind <- which(sd_vec == 0)
    if(length(ind) > 0){
      warning("Found columns with zero sd! Replaced this column with a column of all zeros!")
      submat <- matrix(0, nrow = n, ncol = length(ind))
      x[, ind] <- submat
      attr(x = x, which = "scaled:center")[ind] <- 0
      attr(x = x, which = "scaled:scale")[ind] <- 1
    }
  }
  return(x)
}

#' A general function for calculating the lambda sequence for lasso problem
#' @export
get_lambda <- function(x, y, offset = NULL, intercept = FALSE,
                       family = "gaussian",
                       nlam = 100,
                       lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04)){
  n <- nrow(x)

  # get a path of tuning parameters
  if(family == "ordinal"){
    lam_max <- lam_max_ordinet(x = x, y = y, offset = offset, intercept = intercept)
  }else{
    if (is.null(offset)) offset <- rep(0, n)
    if(inherits(x, "sparseMatrix")){
      # if x is either sparse or binary
      if(family == "gaussian")
        lam_max <- max(abs(Matrix::crossprod(x, y - offset))) / n
      else if(family == "binomial")
        lam_max <- max(abs(Matrix::crossprod(x, y - 1 / (1 + exp(-offset))))) / n
      else if(family == "poisson")
        lam_max <- max(abs(Matrix::crossprod(x, y - exp(offset)))) / n
    }else{
      if(family == "gaussian")
        lam_max <- max(abs(crossprod(x, y - offset))) / n
      else if(family == "binomial")
        lam_max <- max(abs(crossprod(x, y - 1 / (1 + exp(-offset))))) / n
      else if(family == "poisson"){
        lam_max <- max(abs(crossprod(x, y - exp(offset)))) / n
      }
    }
  }

  if(lam_max == 0)
    stop('the calculated maximum of tuning parameter is 0! Please check your response and data matrix.')
  return(lam_max * exp(seq(0, log(lam_min_ratio), length = nlam)))
}

