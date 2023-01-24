#' Running all pairs lasso
#'
#' @import glmnet
#' @export
apl <- function(x, y, family = "gaussian", type.measure = "deviance", ...){
  
  # x is the unstandardized design matrix
  p <- ncol(x)
  n <- nrow(x)
  q <- (p^2 + 3 * p) / 2
  
  # fit cv.glment with all main effects, square effects and all interactions
  x <- myscale(x)
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")
  
  # pairwise index of interactions
  idx <- t(combn(p, 2))
  # construct design matrix of pairwise interactions
  xx <- myscale(cbind(x^2, x[, idx[, 1]] * x[, idx[, 2]]))
  # now idx contains all index pairs
  idx <- rbind(cbind(rep(0, p), seq(p)), cbind(seq(p), seq(p)), idx)
  
  col_mean <- c(col_mean, attr(x = xx, which = "scaled:center"))
  col_sd <- c(col_sd, attr(x = xx, which = "scaled:scale"))
  
  design <- cbind(x, xx)
  
  # run cv.glmnet
  fit <- glmnet::cv.glmnet(x = design, y = y, family = family,
                           type.measure = type.measure,
                           intercept = TRUE,
                           standardize = FALSE, ...)
  
  ibest <- which.min(fit$cvm)
  beta <- as.numeric(fit$glmnet.fit$beta[, ibest])
  # scale estimates back to the original scale of x
  beta <- beta / col_sd
  a0 <- fit$glmnet.fit$a0[ibest] + as.numeric(- crossprod(col_mean, beta))
  
  compact <- cbind(idx[beta != 0, , drop = FALSE], beta[beta != 0])
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

#' Running lasso on main effects only
#' @import glmnet
#' @export
mel <- function(x, y, family = "gaussian", type.measure = "deviance", ...){
  
  n <- nrow(x)
  p <- ncol(x)
  
  x <- myscale(x)
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")
  
  fit <- glmnet::cv.glmnet(x = x, y = y, family = family, 
                           type.measure = type.measure,
                           intercept = TRUE,
                           standardize = FALSE, ...)
  
  ibest <- which.min(fit$cvm)
  beta <- as.numeric(fit$glmnet.fit$beta[, ibest])
  
  # scale estimates back to the original scale of x
  beta <- beta / col_sd
  a0 <- fit$glmnet.fit$a0[ibest] + as.numeric(- crossprod(col_mean, beta))
  
  idx <- cbind(rep(0, p), seq(p))
  compact <- cbind(idx[beta != 0, , drop = FALSE], beta[beta != 0])
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