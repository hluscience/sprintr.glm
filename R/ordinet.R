# We use ordinet, a modification of the ordinalNet R Package created by Wurm, Michael J., Rathouz, Paul J., and Hanlon, Bret M., for ordinal logistic regression.


# Penalty term is in elastic net format
# returns lambda * (alpha*|betaHat|_1 + (1-alpha)/2*|betaHat|_2^2)
getPenalty <- function(betaHat, lambdaMod, alpha)
{
  lambdaMod[betaHat == 0] <- 0  
  l1 <- sum(abs(betaHat) * lambdaMod)
  l22 <- sum(betaHat^2 * lambdaMod)
  pen1 <- alpha * l1
  pen2 <- .5 * (1-alpha) * l22
  pen <- pen1 + pen2
  pen
}

getScoreInfo <- function(x, xMat, yMat, pMat, pMin, intercept)
{
  nobs <- nrow(x)
  nlev <- ncol(yMat)
  
  # Fitted probabilities less than pMin are set to pMin; everything is rescaled to sum to 1
  pMatFull <- cbind(pMat, 1-rowSums(pMat))
  pMatFull[pMatFull < pMin] <- pMin
  pMatFull <- pMatFull / rowSums(pMatFull) # scale to have sum == 1
  pkplusone <- pMatFull[, nlev]  # vector of fitted probabilities for class K+1
  pMat <- pMatFull[, -nlev, drop=FALSE]
  delta <- t(apply(pMat, 1, cumsum))
  
  # Apply block-diagonal sparse matrix multiplication to calculate score = X'[Q'P]
  # where Q is a sparse block-diagonal matrix, with each block diagonal element being the Jacobian of the inverse link for one observation,
  # and P is a vector that concatenates the log-likelihood gradients of n observations
  # Construct P
  d <- yMat[, nlev] / pkplusone
  d[yMat[, nlev] == 0] <- 0  # 0/0 = 0
  uMat <- yMat[, -nlev, drop=FALSE] / pMat
  uMat[yMat[, -nlev] == 0] <- 0  # 0/0 = 0
  uminusdMat <- uMat - d
  P <- c(t(uminusdMat))
  
  # Construct Q
  ttipDiag <- diag(1, nrow=nlev-1)
  ttipOffDiag <- cbind(-ttipDiag[, -1], 0)
  ttip <- ttipDiag + ttipOffDiag
  Q <- Matrix::t(Matrix::bdiag(lapply(split(matrix(rep(c(t(delta * (1-delta))), each=nlev-1) * rep(ttip, times = nobs), ncol = nlev - 1, byrow = TRUE), rep(c(1:nobs), each=nlev-1)), matrix, nrow = nlev-1)))
  
  # Calculate score = X'[Q'P]
  prod <- Matrix::crossprod(Q, P) 
  scoreIntercept <- rowSums(matrix(prod, nrow = nlev - 1, byrow = FALSE))
  scoreNoneIntercept <- as.vector(Matrix::t(x) %*% colSums(matrix(prod, nrow = nlev - 1, byrow = FALSE)))
  score <- if(intercept) c(scoreIntercept, scoreNoneIntercept) else scoreNoneIntercept
  
  # Apply block-diagonal sparse matrix multiplication to calculate W = Q'sigInvQ
  wts <- if (is.null(attr(yMat, "wts"))) rowSums(yMat) else attr(yMat, "wts")
  wpMat <- wts / pMat
  wpkplusone <- wts / pkplusone
  sigInv <- Matrix::Diagonal(x=c(t(wpMat))) +
    Matrix::bdiag(lapply(split(matrix(rep(wpkplusone, each=(nlev-1)*(nlev-1)), ncol = nlev - 1, byrow=TRUE), rep(c(1:nobs),each=nlev-1)), matrix, nrow = nlev-1))
  W <- Matrix::crossprod(Q, sigInv) %*% Q
  
  # W is a block-diagonal sparse matrix that can be utilized to get fisher information by using the matrix sum operation 
  infoIntercept <- matrix(0, nrow = nlev-1, ncol = nlev-1)
  widx <- 1:(nobs*(nlev-1))
  for (iLev in 1:nlev-1){
    for (jLev in 1:nlev-1){
      if (iLev == nlev - 1) i <- 0 else i <- iLev
      if (jLev == nlev - 1) j <- 0 else j <- jLev
      ix <- widx %% (nlev-1)
      infoIntercept[iLev, jLev] <- sum(W[widx[which(ix == i)], widx[which(ix == j)]])
    }
  }
  wRsumMat <- matrix(Matrix::rowSums(W), nrow = nlev - 1, byrow = FALSE)
  wTSumSeq <- colSums(wRsumMat)
  InfoCross <- wRsumMat%*%x
  InfoNoneIntercept <- Matrix::t(wTSumSeq*x)%*%x
  info <- if(intercept) rbind(cbind(infoIntercept, InfoCross), cbind(Matrix::t(InfoCross),InfoNoneIntercept)) else InfoNoneIntercept
  
  list(score=score, info=as.matrix(info))
  
}

getLoglik <- function(pMat, yMat)
{
  pkplusone <- 1 - rowSums(pMat)
  pMatFull <- cbind(pMat, pkplusone)
  if (any(pMatFull < 0)) return(-Inf)
  llMat <- yMat * log(pMatFull)
  llMat[yMat==0] <- 0  # -Inf*0 = 0
  llik <- sum(llMat)
  llik
}

# Returns approximate log-likelihood (as a function of beta, up to a constant)
getLoglikApprox <- function(betaHatActive, scoreActive, infoActive)
{
  -sum(betaHatActive * (infoActive %*% betaHatActive)) + sum(scoreActive * betaHatActive)
}

getDeltaNames <- function(linkFamily, reverse, nlev)
{
  index <- if (reverse) nlev:2 else 1:(nlev-1)
  deltaNames <- sapply(index, function(i)
  {
    if (linkFamily=="cumulative") {
      if (reverse) {
        paste0("P[Y>=", i, "]")  # P(Y>=i)
      } else {
        paste0("P[Y<=", i, "]")  # P(Y<=i)
      }
    }  
  })
  
  deltaNames
}

softThresh <- function(z, g) sign(z)*max(0, abs(z)-g)

invLogit <- function(x) 1 / (1+exp(-x))

yFactorToMatrix <- function(y)
{ 
  nobs <- length(y)
  nlev <- length(levels(y))
  yMat <- matrix(0, nrow=nobs, ncol=nlev, dimnames=list(NULL, levels(y)))
  yInt <- as.integer(y)
  yMat[cbind(1:nobs, yInt)] <- 1
  yMat
}

# Creates a list of functions consisting of
# g (link), h (inverse link), and getQ (Jacobian of inverse link, dh/ddeta^T).
# delta is a transformation of p to which the element-wise link function is applied (e.g. logit)
# tt(p) = delta, ttinv(delta) = p, and ttinvprime(delta) = dttinv/ddelta^T
makeLinkfun <- function(linkFamily, link)
{
  if (linkFamily == "cumulative") {
    linkfun <- makeLinkCumulative(link)
  }
  linkfun
}

# Cumulative probability family
makeLinkCumulative <- function(link)
{
  lf <- stats::make.link(link)
  
  # prob -> cumulative prob
  tt <- function(p) cumsum(p)
  
  # cumulative prob -> prob
  ttinv <- function(delta) delta - c(0, delta[-length(delta)])
  
  ttinvprime <- function(delta)
  {
    k <- length(delta)
    ttipDiag <- diag(1, nrow=k)
    ttipOffDiag <- cbind(-ttipDiag[, -1], 0)
    ttip <- ttipDiag + ttipOffDiag
    ttip
  }
  
  g <- function(p) lf$linkfun(tt(p))
  h <- function(eta) ttinv(lf$linkinv(eta))
  getQ <- function(eta) rep(lf$mu.eta(eta), each=length(eta)) * ttinvprime(lf$linkinv(eta))
  
  list(g=g, h=h, getQ=getQ)
}

# Coordinate descent inner loop function
# Note: should not need to check for improvement because each coordinate step necessarily improves the approximate objective
cdIn <- function(wtsum, betaHat, score, info, alpha, lambdaMod, positiveID, threshIn, maxiterIn)
{
  # Update Active Set
  activeSet <- which((betaHat!=0 | lambdaMod==0) & diag(info)!=0)
  # check lambdaMod==0 because if data are balanced, an intercept could be initialized exactly to zero
  if (length(activeSet) == 0) betaHat <- rep(0, length(betaHat)) else betaHat[-activeSet] <- 0
  betaHatActive <- betaHat[activeSet]
  scoreActive <- score[activeSet]
  infoActive <- info[activeSet, activeSet, drop=FALSE]
  infoInactive <- info[-activeSet, activeSet, drop=FALSE]
  lambdaModActive <- lambdaMod[activeSet]
  lambdaModInactive <- if (length(activeSet) == 0) lambdaMod else lambdaMod[-activeSet]
  positiveIDActive <- positiveID[activeSet]
  positiveIDInactive <- if (length(activeSet) == 0) positiveID else positiveID[-activeSet]

  # softThreshTerms vector does not change during the inner loop, even if active set changes
  softThreshTerms <- c(info[, activeSet, drop=FALSE] %*% betaHatActive) + score
  # softThreshTerms = I(beta^(r)) %*% beta^(r) + U(beta^(r))
  softThreshTermsActive <- softThreshTerms[activeSet]
  softThreshTermsInactive <- if (length(activeSet) == 0) softThreshTerms else softThreshTerms[-activeSet]

  # Initialize quadratic approximation to the log-likelihood and objective
  loglikApprox <- getLoglikApprox(betaHatActive, scoreActive, infoActive)
  penalty <- getPenalty(betaHat, lambdaMod, alpha)
  objApprox <- -loglikApprox/wtsum + penalty

  iterIn <- 0
  kktAll <- FALSE
  while (!kktAll && iterIn<maxiterIn)
  {
    conv <- FALSE
    while (!conv && iterIn<maxiterIn)
    {
      iterIn <- iterIn + 1
      for (i in seq_along(activeSet))
      {
        numTerm <- softThreshTermsActive[i] - sum(infoActive[i, -i, drop=FALSE] * betaHatActive[-i])
        denTerm <- infoActive[i, i]
        penTerm <- wtsum * lambdaModActive[i]
        betaHatActive[i] <- softThresh(numTerm, penTerm*alpha) / (denTerm + penTerm*(1-alpha))
        if (positiveIDActive[i]) betaHatActive[i] <- max(0, betaHatActive[i])
      }

      # get_lambda step 1: betaHatActive[i] <- (infoActive[i, i] * betaHatActive[i] + score)/infoActive[i, i]
      betaHatOld <- betaHat
      betaHat[activeSet] <- betaHatActive
      loglikApprox <- getLoglikApprox(betaHatActive, scoreActive, infoActive)
      penalty <- getPenalty(betaHat, lambdaMod, alpha)
      objApproxOld <- objApprox
      objApprox <- -loglikApprox / wtsum + penalty
      dif <- abs((objApproxOld - objApprox) / (abs(objApproxOld) + 1e-100))
      conv <- dif < threshIn

    }  # end while (!conv && iterIn<maxiterIn)

    kkt <- rep(TRUE, length(betaHat))
    kktInactiveTerms <- if (length(activeSet) == 0) softThreshTermsInactive else softThreshTermsInactive - c(infoInactive %*% betaHatActive)
    kktInactiveTerms[!positiveIDInactive] <- abs(kktInactiveTerms[!positiveIDInactive])
    kktInactive <- kktInactiveTerms <= wtsum * lambdaModInactive * alpha
    if (length(activeSet) == 0) kkt <- kktInactive else kkt[-activeSet] <- kktInactive
    kktAll <- all(kkt)
    if (!kktAll)
    {
      iterIn <- iterIn - 1  # repeat the iteration if kkt conditions are not satisfied
      activeSet <- union(activeSet, which(!kkt))
      betaHatActive <- betaHat[activeSet]
      scoreActive <- score[activeSet]
      infoActive <- info[activeSet, activeSet, drop=FALSE]
      infoInactive <- info[-activeSet, activeSet, drop=FALSE]
      softThreshTermsActive <- softThreshTerms[activeSet]
      softThreshTermsInactive <- softThreshTerms[-activeSet]
      lambdaModActive <- lambdaMod[activeSet]
      lambdaModInactive <- lambdaMod[-activeSet]
      positiveIDActive <- positiveID[activeSet]
      positiveIDInactive <- positiveID[-activeSet]
    }

  }  # end while (!kktAll && iterIn<maxiterIn)

  list(betaHat=betaHat, iterIn=iterIn)
}

# Coordinate descent outer loop function
cdOut <- function(betaHat, lambdaMod, positiveID,
                  x, xMat, yMat, offset, intercept, alpha, linkfun,
                  pMin, threshOut, threshIn, maxiterOut, maxiterIn)
{
  nobs <- nrow(yMat)
  wts <- if (is.null(attr(yMat, "wts"))) rowSums(yMat) else attr(yMat, "wts")
  wtsum <- if (is.null(attr(yMat, "wtsum"))) sum(wts) else attr(yMat, "wtsum")

  # Could carry over epsMat, pMat, and loglik from previous cdOut iteration
  # betaHat = betaStart from mirlsNet
  betaNonzeroIndex <- which(betaHat != 0)
  etaMat <- matrix(xMat[, betaNonzeroIndex, drop = FALSE] %*% betaHat[betaNonzeroIndex], nrow = nobs, byrow = TRUE) + offset
  pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
  loglik <- getLoglik(pMat, yMat)
  # when calculate lambda_max, penalty = 0
  penalty <- getPenalty(betaHat, lambdaMod, alpha)
  obj <- -loglik/wtsum + penalty

  conv <- FALSE
  iterOut <- 0
  while (!conv && iterOut<maxiterOut)
  {
    iterOut <- iterOut + 1

    # Update score and info
    si <- getScoreInfo(x, xMat, yMat, pMat, pMin, intercept)
    score <- si$score
    info <- si$info

    # Run coordinate descent inner loop
    betaHatOld <- betaHat
    cdInResult <- cdIn(wtsum, betaHat, score, info, alpha, lambdaMod, positiveID, threshIn, maxiterIn)
    betaHat <- cdInResult$betaHat
    iterIn <- cdInResult$iterIn
    betaNonzeroIndex <- which(betaHat != 0)
    etaMat <- matrix(xMat[, betaNonzeroIndex, drop = FALSE] %*% betaHat[betaNonzeroIndex], nrow = nobs, byrow = TRUE) + offset
    pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))

    # Update log-likelihood and objective
    loglikOld <- loglik
    objOld <- obj
    loglik <- getLoglik(pMat, yMat)
    penalty <- getPenalty(betaHat, lambdaMod, alpha)
    obj <- -loglik / wtsum + penalty

    # Take half steps if obj does not improve. Loglik is set to -Inf
    # if any fitted probabilities are negative, which can happen for
    # the nonparallel or semiparallel cumulative probability model.
    nhalf <- 0
    while (obj > objOld && nhalf < 10) {
      nhalf <- nhalf + 1
      betaHat <- (betaHat + betaHatOld) / 2
      betaNonzeroIndex <- which(betaHat != 0)
      etaMat <- matrix(xMat[, betaNonzeroIndex, drop=FALSE] %*% betaHat[betaNonzeroIndex], nrow=nobs, byrow=TRUE) + offset
      pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
      loglik <- getLoglik(pMat, yMat)
      penalty <- getPenalty(betaHat, lambdaMod, alpha)
      obj <- -loglik / wtsum + penalty
    }
    dif <- (objOld - obj) / (abs(objOld) + 1e-100)
    conv <- dif < threshOut

    # Convergence is declared if objective worsens. In this case, set betaHat
    # to previous value. (Typically means model is saturated.)
    objImproved <- obj <= objOld
    if (!objImproved)
    {
      betaHat <- betaHatOld
      loglik <- loglikOld
    }
  }  # end while (!conv && iterOut<maxiterOut)

  # Opting not to return penalty or obj because they depend on covariate scaling.
  list(betaHat=betaHat, loglik=loglik, iterOut=iterOut, iterIn=iterIn, dif=dif)
}

# General optimization algorithm for multinomial regression models via coordinate descent
mirlsNet <- function(x, xMat, yMat, offset, intercept, alpha,
                     penaltyFactors, positiveID, linkfun, betaStart,
                     lambdaVals, nlam, lam_min_ratio, includeLambda0, alphaMin, pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn)
{
  nobs <- nrow(yMat)
  wtsum <- if (is.null(attr(yMat, "wtsum"))) sum(yMat) else attr(yMat, "wtsum")
  if (!is.null(lambdaVals)) lambdaVals <- sort(lambdaVals, decreasing=TRUE)
  lambdaNum <- if (is.null(lambdaVals)) nlam + includeLambda0 else length(lambdaVals)
  fits <- vector("list", length=lambdaNum)

  # If lambdaVals=NULL, need to determine the minimum lambda value that sets all penalized coefficients to zero
  if (is.null(lambdaVals))
  {
    lambdaMod <- ifelse(penaltyFactors==0, 0, Inf)  # to find solution with only unpenalized terms
    fits[[1]] <- cdOut(betaHat=betaStart, lambdaMod, positiveID,
                       x, xMat, yMat, offset, intercept, max(alpha, alphaMin), linkfun,
                       pMin, threshOut, threshIn, maxiterOut, maxiterIn)

    betaStart <- fits[[1]]$betaHat

    # Calculate starting lambda value
    etaMat <- matrix(xMat %*% betaStart, nrow=nobs, byrow=TRUE) + offset
    pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
    si <- getScoreInfo(x, xMat, yMat, pMat, pMin, intercept)
    # betaHat is zero for all penalized terms, so the soft threshold argument is just the score function
    penID <- penaltyFactors != 0
    lambdaMaxVals <- si$score[penID] / (wtsum * max(alpha, alphaMin) * penaltyFactors[penID])
    lambdaMaxVals[positiveID[penID]] <- pmax(0, lambdaMaxVals[penID & positiveID])
    lambdaMaxVals <- abs(lambdaMaxVals)
    lambdaMax <- max(lambdaMaxVals)
    lambdaMin <- lambdaMax * lam_min_ratio
    lambdaVals <- exp(seq(log(lambdaMax), log(lambdaMin), length.out=nlam))
    if (includeLambda0) lambdaVals <- c(lambdaVals, 0)
  }

  # If alpha < alphaMin, the model needs to be re-fit the first lambda value
  # using alpha instead of alphaMin
  if (alpha < alphaMin)
    fits[1] <- list(NULL)

  # fits[[1]] is NULL if alpha < alphaMin or if lambdaVals is specified by user.
  llik <- if (is.null(fits[[1]])) -Inf else fits[[1]]$loglik
  for (i in (1+!is.null(fits[[1]])):lambdaNum)
  {
    # If relative change in loglik is < stopThresh, then use the current fit
    # for all remaining lambda. Do not stop if loglik stays the same, because
    # this can happen if the first several lambda values produce null models,
    # e.g. in cross validation.
    if ((i > 2) && llikOld != llik && (abs((llikOld - llik) / llikOld) < stopThresh))
    {
      fits[[i]] <- fits[[i-1]]
    } else
    {
      lambdaMod <- lambdaVals[i] * penaltyFactors
      lambdaMod <- ifelse(penaltyFactors==0, 0, lambdaVals[i] * penaltyFactors)
      fits[[i]] <- cdOut(betaHat=betaStart, lambdaMod, positiveID,
                         x, xMat, yMat, offset, intercept, alpha, linkfun,
                         pMin, threshOut, threshIn, maxiterOut, maxiterIn)

      betaStart <- fits[[i]]$betaHat
      llikOld <- llik
      llik <- fits[[i]]$loglik
    }
  }

  iterOut <- sapply(fits, function(x) x$iterOut)
  iterIn <- sapply(fits, function(x) x$iterIn)
  dif <- sapply(fits, function(x) x$dif)
  betaHat <- t(sapply(fits, function(f) f$betaHat))
  loglik <- sapply(fits, function(f) f$loglik)
  list(lambdaVals=lambdaVals, betaHat=betaHat, loglik=loglik, iterOut=iterOut, iterIn=iterIn)
}

# Ordinal regression models with elastic net penalty
ordinet <- function(x, y, offset = NULL, alpha=1, intercept = TRUE, standardize=TRUE, penaltyFactors=NULL, positiveID=NULL, linkFamily = "cumulative", reverse=FALSE, link = "logit", parallelTerms=TRUE, nonparallelTerms=FALSE, parallelPenaltyFactor=1, lambdaVals=NULL, nlam=20, lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), includeLambda0=FALSE, alphaMin=0.01, pMin=1e-8, stopThresh=1e-8, threshOut=1e-8, threshIn=1e-8, maxiterOut=100, maxiterIn=100) 
{
  
  args <- as.list(environment()) 
  # Initial argument checks
  if (!is.factor(y) && !is.matrix(y))
    stop("y should be a factor or matrix in ordinal family")
  
  # Variable definitions
  yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)
  wts <- attr(yMat, "wts") <- rowSums(yMat)
  wtsum <- attr(yMat, "wtsum") <- sum(wts)
  nvar <- ncol(x) # nvar = 20
  nobs <- nrow(x)
  nlev <- ncol(yMat) # nlev = 5
  if (reverse) yMat <- yMat[, nlev:1] # backward/forward
  if (is.null(offset)) offset <- matrix(0, nrow = nobs, ncol = nlev-1)
  
  # Create linkfun
  linkfun <- makeLinkfun(linkFamily, link)
  
  if (inherits(x, "sparseMatrix") | standardize == FALSE) {
    # if x is either sparse or binary | standardize == FALSE
    xMeans <- rep(0, ncol(x))
    xSD <- rep(1, ncol(x))
    xStd <- x
  } else {
    xStd <- scale(x, center = TRUE, scale = TRUE)
    # for those columns with 0 sds
    # add small Gaussian noise to it
    xMeans <- attr(x = xStd, which = "scaled:center")
    xSD <- attr(x = xStd, which = "scaled:scale")
  }
  
  x1 <- diag(nlev-1)  # intercept columns
  xList <- lapply(1:nrow(x), function(i)
  {
    xi <- x[i, ]
    x2 <- if (!parallelTerms) NULL else rbind(xi)[rep(1, nlev-1), , drop=FALSE]
    x3 <- if (!nonparallelTerms) NULL else makeNonparallelBlock(xi, nlev)
    xListi <- if(intercept) cbind(x1, x2, x3) else cbind(x2, x3)
    rownames(xListi) <- NULL
    xListi
  })
  xMat <- if(inherits(x, "sparseMatrix")) Matrix::Matrix(do.call(rbind, xList), sparse = TRUE) else do.call(rbind, xList)
  
  # Augment penaltyFactors to include all model coefficients
  if (is.null(penaltyFactors)) penaltyFactors <- rep(1, nvar)
  penaltyFactorsParallel <- if (parallelTerms) penaltyFactors * parallelPenaltyFactor else NULL
  penaltyFactorsNonparallel <- if(nonparallelTerms) rep(penaltyFactors, nlev-1) else NULL
  penaltyFactors <- if(intercept) c(rep(0, nlev-1), penaltyFactorsParallel, penaltyFactorsNonparallel) else c(penaltyFactorsParallel, penaltyFactorsNonparallel)
  
  # Augment positiveID to include all model coefficients
  if (is.null(positiveID)) positiveID <- rep(FALSE, nvar)
  positiveID <- if(intercept) c(rep(FALSE, nlev-1), rep(positiveID, parallelTerms + nonparallelTerms*(nlev-1))) else c(rep(positiveID, parallelTerms + nonparallelTerms*(nlev-1)))
  
  if(intercept)
  {
    # Initialize coefficient values to intercept-only model
    yFreq <- colSums(yMat) / wtsum
    interceptStart <- linkfun$g(yFreq[-nlev]) 
    interceptStart <- pmin(100, pmax(-100, interceptStart))
    noninterceptStart <- rep(0, nvar*(parallelTerms + nonparallelTerms*(nlev-1)))
    betaStart <- c(interceptStart, noninterceptStart)
  }else{
    betaStart <- rep(0, nvar*(parallelTerms + nonparallelTerms*(nlev-1)))
  }
  
  
  # Fit solution path
  mirlsNetFit <- mirlsNet(x, xMat, yMat, offset, intercept, alpha, 
                          penaltyFactors, positiveID, linkfun, betaStart,
                          lambdaVals, nlam, lam_min_ratio, includeLambda0, alphaMin,
                          pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn) 
  
  betaHat <- mirlsNetFit$betaHat
  lambdaVals <- mirlsNetFit$lambdaVals
  loglik <- mirlsNetFit$loglik
  iterOut <- mirlsNetFit$iterOut
  iterIn <- mirlsNetFit$iterIn
  
  # Change coefficient estimates back to original scale if standardize=TRUE
  intercepts0 <- if (intercept) betaHat[ , 1:(nlev-1), drop=FALSE] else NULL 
  nonintercepts0 <- if (intercept) betaHat[ , -(1:(nlev-1)), drop=FALSE] else betaHat
  unscaleFact <- xMeans / xSD 
  intAdjust <- matrix(0, nrow=nrow(betaHat), ncol=nlev-1)
  if (parallelTerms) intAdjust <- intAdjust +
    (nonintercepts0[ , 1:nvar, drop=FALSE] %*% unscaleFact)[ , rep(1, nlev-1), drop=FALSE]
  if (nonparallelTerms) intAdjust <- intAdjust + sapply(1:(nlev-1), function(i) {
    nonintercepts0[ , (nvar*(i-1+parallelTerms)+1):(nvar*(i+parallelTerms)), drop=FALSE] %*% unscaleFact
  })
  intercepts <- intercepts0 - intAdjust 
  nonintercepts <- if (standardize) t(t(nonintercepts0) / xSD) else nonintercepts0
  coefs <- cbind(intercepts, nonintercepts)
  
  # Create coefficient column names
  catOrder <- if (reverse) nlev:2 else 1:(nlev-1)
  interceptNames <- if (intercept) paste0("(Intercept):", catOrder) else NULL
  xNames <- if (is.null(colnames(x))) paste0("X", 1:nvar) else colnames(x)
  parallelNames <- nonparallelNames <- NULL
  if (parallelTerms) parallelNames <- xNames
  if (nonparallelTerms) nonparallelNames <- paste0(rep(xNames, nlev-1), ":", rep(catOrder, each=nvar))
  colnames(coefs) <- c(interceptNames, parallelNames, nonparallelNames)
  
  fit <- list(coefs=coefs, lambdaVals=lambdaVals, loglik=loglik,
              nlev=nlev, nvar=nvar, xNames=xNames, args = args, iterOut=iterOut, iterIn=iterIn) 
  return(fit)
}

# Uses K-fold cross validation to obtain out-of-sample log-likelihood
# Lambda is tuned within each cross validation fold.
cv.ordinet <-  function(x, y, offset = NULL, alpha=1, intercept = TRUE, standardize=FALSE, penaltyFactors=NULL, positiveID=NULL, linkFamily="cumulative", reverse=FALSE, link="logit", parallelTerms=TRUE, nonparallelTerms=FALSE, parallelPenaltyFactor=1, lambdaVals=NULL, nlam=20, lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04), includeLambda0=FALSE, alphaMin=0.01, Min=1e-8, stopThresh=1e-8, threshOut=1e-8, threshIn=1e-8, maxiterOut=100, maxiterIn=100, folds=NULL, nFolds=5) 
{
  
  yMat <- if (is.matrix(y)) y else yFactorToMatrix(y)  
  fit <- ordinet(x, y, offset, alpha, intercept, standardize, 
                 penaltyFactors, positiveID, linkFamily, reverse,
                 link,
                 parallelTerms, nonparallelTerms, parallelPenaltyFactor,
                 lambdaVals, nlam, lam_min_ratio, includeLambda0, alphaMin,
                 pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn)
  
  if (is.null(lambdaVals)) lambdaVals <- fit$lambdaVals
  
  if (is.null(folds))
  {
    n <- nrow(x)
    randIndex <- sample(n)
    folds <- split(randIndex, rep(1:nFolds, length.out=n))
  } else
  {
    nFolds <- length(folds)
  }
  
  nlam <- length(lambdaVals)
  loglik <- matrix(nrow=nlam, ncol=nFolds)
  colnames(loglik) <- paste0("fold", 1:nFolds)
  rownames(loglik) <- paste0("lambda", 1:nlam)
  for (i in 1:nFolds)
  {
    testFold <- folds[[i]]
    xTrain <- x[-testFold, , drop=FALSE]
    xTest <- x[testFold, , drop=FALSE]
    yTrain <- if (is.matrix(y)) y[-testFold, , drop=FALSE] else y[-testFold]
    yMatTest <- yMat[testFold, , drop=FALSE]
    fitTrain <- ordinet(xTrain, yTrain, offset, alpha, intercept, standardize, 
                        penaltyFactors, positiveID, linkFamily, reverse,
                        link, 
                        parallelTerms, nonparallelTerms, parallelPenaltyFactor,
                        lambdaVals, nlam, lam_min_ratio, includeLambda0, alphaMin,
                        pMin, stopThresh, threshOut, threshIn, maxiterOut, maxiterIn)
    
    for (j in 1:nlam){
      pHatFull <- predict.ordinet(object=fitTrain, newx=xTest, whichLambda=j)
      pHat <- pHatFull[, -ncol(pHatFull), drop=FALSE]
      loglik[j, i] <- getLoglik(pHat, yMatTest)
    } 
  }    
  
  bestLambdaIndex <- as.numeric(which.max(rowMeans(loglik*is.finite(loglik),na.rm=TRUE)))
  
  out <- list(loglik=loglik, bestLambdaIndex=bestLambdaIndex, lambdaVals=lambdaVals, folds=folds, fit=fit)
  return(out)
}

# Get a sequence of lambda
lam_max_ordinet <- function(x, y, offset = NULL, alpha=1, intercept = TRUE, penaltyFactors=NULL, positiveID=NULL, linkFamily="cumulative", link="logit", parallelTerms=TRUE, nonparallelTerms=FALSE, parallelPenaltyFactor=1, alphaMin=0.01, pMin=1e-8, stopThresh=1e-8, threshOut=1e-8, threshIn=1e-8, maxiterOut=100, maxiterIn=100) 
{
  
  # Variable definitions
  yMat <- yFactorToMatrix(y)
  wts <- attr(yMat, "wts") <- rowSums(yMat)
  wtsum <- attr(yMat, "wtsum") <- sum(wts)
  nobs <- nrow(x)
  nvar <- ncol(x) 
  nlev <- ncol(yMat) 
  
  # Create linkfun
  linkfun <- makeLinkfun(linkFamily, link)
  
  x1 <- diag(nlev-1)  # intercept columns
  xList <- lapply(1:nrow(x), function(i)
  {
    xi <- x[i, ]
    x2 <- if (!parallelTerms) NULL else rbind(xi)[rep(1, nlev-1), , drop=FALSE]
    x3 <- if (!nonparallelTerms) NULL else makeNonparallelBlock(xi, nlev)
    xListi <- if(intercept) cbind(x1, x2, x3) else cbind(x2, x3)
    rownames(xListi) <- NULL
    xListi
  })
  xMat <- if(inherits(x, "sparseMatrix")) Matrix::Matrix(do.call(rbind, xList), sparse = TRUE) else do.call(rbind, xList)
  
  # Augment penaltyFactors to include all model coefficients
  if (is.null(penaltyFactors)) penaltyFactors <- rep(1, nvar)
  penaltyFactorsParallel <- if (parallelTerms) penaltyFactors * parallelPenaltyFactor else NULL
  penaltyFactorsNonparallel <- if(nonparallelTerms) rep(penaltyFactors, nlev-1) else NULL
  penaltyFactors <- if(intercept) c(rep(0, nlev-1), penaltyFactorsParallel, penaltyFactorsNonparallel) else c(penaltyFactorsParallel, penaltyFactorsNonparallel)
  
  # Augment positiveID to include all model coefficients
  if (is.null(positiveID)) positiveID <- rep(FALSE, nvar)
  positiveID <- if(intercept) c(rep(FALSE, nlev-1), rep(positiveID, parallelTerms + nonparallelTerms*(nlev-1))) else c(rep(positiveID, parallelTerms + nonparallelTerms*(nlev-1)))
  
  if(intercept)
  {
    yFreq <- colSums(yMat) / wtsum
    interceptStart <- linkfun$g(yFreq[-nlev]) 
    interceptStart <- pmin(100, pmax(-100, interceptStart))
    noninterceptStart <- rep(0, nvar*(parallelTerms + nonparallelTerms*(nlev-1)))
    betaStart <- c(interceptStart, noninterceptStart)
  }else{
    betaStart <- rep(0, nvar*(parallelTerms + nonparallelTerms*(nlev-1)))
  }
  
  # Determine the minimum lambda value that sets all penalized coefficients to zero
  if (is.null(offset))
  {
    offset <- matrix(0, nrow = nobs, ncol = nlev-1)
    lambdaMod <- ifelse(penaltyFactors==0, 0, Inf)  # to find solution with only unpenalized terms
    fit <- cdOut(betaHat = betaStart, lambdaMod, positiveID,
                 x, xMat, yMat, offset, intercept, max(alpha, alphaMin), linkfun,
                 pMin, threshOut, threshIn, maxiterOut, maxiterIn) 
    betaStart <- fit$betaHat
    
    # Calculate starting lambda value
    etaMat <- matrix(xMat %*% betaStart, nrow=nobs, byrow=TRUE)
  }else{
    etaMat <- matrix(xMat %*% betaStart, nrow=nobs, byrow=TRUE) + offset
  }
  
  pMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
  si <- getScoreInfo(x, xMat, yMat, pMat, pMin, intercept)
  # betaHat is zero for all penalized terms, so the soft threshold argument is just the score function
  penID <- penaltyFactors != 0
  lambdaMaxVals <- si$score[penID] / (wtsum * max(alpha, alphaMin) * penaltyFactors[penID])
  lambdaMaxVals[positiveID[penID]] <- pmax(0, lambdaMaxVals[penID & positiveID])
  lambdaMaxVals <- abs(lambdaMaxVals)
  lambdaMax <- max(lambdaMaxVals)
  
  return(lambdaMax)
}

# Make predictions for "ordinet" object
predict.ordinet <- function(object, newx=NULL, whichLambda=NULL)
{
  # Extract variables from ordinet object
  nlev <- object$nlev
  nvar <- object$nvar
  xNames <-  object$xNames
  parallelTerms <- object$args$parallelTerms
  nonparallelTerms <- object$args$nonparallelTerms
  reverse <- object$args$reverse
  linkFamily <- object$args$linkFamily
  link <- object$args$link
  linkfun <- makeLinkfun(linkFamily, link)
  
  # Create coefficient matrix
  betaHat <- object$coefs[whichLambda, ]
  intercepts <- betaHat[1:(nlev-1)]
  nonintercepts <- matrix(0, nrow=nvar, ncol=nlev-1)
  if (parallelTerms) nonintercepts <- nonintercepts + betaHat[nlev:(nlev-1+nvar)]
  if (nonparallelTerms) nonintercepts <- nonintercepts + betaHat[-(1:(nlev-1+nvar*parallelTerms))]
  betaMat <- rbind(intercepts, nonintercepts)
  rownames(betaMat) <- c("(Intercept)", xNames)
  deltaNames <- getDeltaNames(linkFamily, reverse, nlev)
  colnames(betaMat) <- paste0(link, "(", deltaNames, ")")
  
  # Compute prediction values
  etaMat <- cbind(1, newx) %*% betaMat
  deltaNames <- getDeltaNames(linkFamily, reverse, nlev)
  colnames(etaMat) <- colnames(betaMat)
  probMat <- do.call(rbind, lapply(1:nrow(etaMat), function(i) linkfun$h(etaMat[i, ])))
  probMat <- cbind(probMat, 1-rowSums(probMat))
  if (reverse) probMat <- probMat[, nlev:1]
  colnames(probMat) <- paste0("P[Y=", 1:nlev, "]")
  return(probMat)
}

sprinter.ordinet <- function(x, y, nlev = 1, square = FALSE, num_keep = NULL, lambda1 = NULL, lambda3 = NULL, cv_step1 = FALSE, nlam1 = 10, nlam3 = 100, lam_min_ratio = ifelse(nrow(x) < ncol(x), 0.01, 1e-04)){
  
  n <- nrow(x)
  p <- ncol(x)
  q <- ifelse(square, p * (p - 1) / 2, p * (p - 1) / 2 + p)
  if (!is.factor(y)) stop("y should be a factor in ordinal family")
  nlev <- nlevels(y)
  
  if(is.null(num_keep)){
    num_keep <- ceiling(n / log(n))
  }else{
    stopifnot(num_keep > 0 & num_keep <= q)}
  
  # we always standardize the design matrix to get main effects
  # interactions are built upon standardized main effects
  x <- myscale(x)
  # xm is the (standardized) design matrix of main effects
  col_mean <- attr(x = x, which = "scaled:center")
  col_sd <- attr(x = x, which = "scaled:scale")
  xm <- x
  
  # The First Step
  # run lasso on
  # (1) main effects M (square == FALSE)
  # (2) main effects M + squared main effects M^2 (square == TRUE)
  if(square){
    x_sq <- myscale(x^2)
    col_mean <- c(col_mean, attr(x_sq, which = "scaled:center"))
    col_sd <- c(col_sd, attr(x_sq, which = "scaled:scale"))
    x <- cbind(x, x_sq)
  }
  
  # construct xMat with intercept columns
  x1 <- diag(nlev-1) # intercept columns
  xList <- lapply(1:nrow(x), function(i)
  {
    xi <- x[i, ]
    xListi <- cbind(x1, rbind(xi)[rep(1, nlev-1), , drop=FALSE])
    rownames(xListi) <- NULL
    xListi
  })
  xMat <- if(inherits(x, "sparseMatrix")) Matrix::Matrix(do.call(rbind, xList), sparse = TRUE) else do.call(rbind, xList)
  # transform y yo yMat with dimension n * nlev
  yMat <- Matrix::Matrix(yFactorToMatrix(y), sparse = TRUE)
  
  # initiate lambda 1
  if(is.null(lambda1)){
    lambda1 <- get_lambda(x = x, y = y, family = "ordinal", intercept = TRUE, nlam = nlam1, lam_min_ratio = lam_min_ratio)
  }
  
  
  # Step 1
  if(cv_step1){
    fit <- cv.ordinet(x, y, lambdaVals = lambda1, intercept = TRUE, standardize=FALSE)
    
    # grab coefficient estimate
    theta <- matrix(fit$fit$coefs[fit$bestLambdaIndex, ], ncol = 1)
    colnames(theta) <- NULL
    rownames(theta) <- NULL
    
    # fitted value of mu, i.e., the linear part
    mu <- matrix(xMat %*% theta, ncol = 1)
    
    # Change coefficient estimates back to original scale if standardize=TRUE
    intercepts0 <- matrix(theta[1:(nlev-1)], ncol = 1)
    nonintercepts0 <- matrix(theta[-(1:(nlev-1))], ncol = 1)
    unscaleFact <- col_mean / col_sd
    intAdjust <- matrix((unscaleFact %*% nonintercepts0)[rep(1, nlev-1)], ncol = 1)
    intercepts <- intercepts0 - intAdjust
    rownames(intercepts) <- paste0("Intercept:", 1:(nlev-1))
    nonintercepts <- nonintercepts0 / col_sd
    
    # update lambda1
    lambda1 <- fit$lambdaVals[fit$bestLambdaIndex]
    
  }else{
    fit <- ordinet(x = x, y = y, lambdaVals = lambda1, intercept = TRUE, standardize = FALSE)
    
    # grab coefficient estimate
    theta <- t(fit$coefs)
    colnames(theta) <- NULL
    rownames(theta) <- NULL
    
    # fitted mu 
    mu <- as.matrix(xMat %*% theta)
    
    # Change coefficient estimates back to original scale if standardize=TRUE
    intercepts0 <- theta[1:(nlev-1), , drop=FALSE] # nrow = nlev - 1, ncol = nlam
    nonintercepts0 <- theta[-(1:(nlev-1)), , drop=FALSE] # nrow = p, ncol = nlam
    unscaleFact <- col_mean / col_sd
    intAdjust <- (unscaleFact %*% nonintercepts0)[rep(1, nlev-1), , drop=FALSE]
    intercepts <- intercepts0 - intAdjust 
    rownames(intercepts) <- paste0("Intercept:", 1:(nlev-1))
    nonintercepts <- nonintercepts0 / col_sd
  }
  
  # update lambda1 and initiate lambda3
  nlam1 <- length(lambda1)
  if(is.null(lambda3))
    lambda3 <- matrix(NA, nlam3, nlam1)
  
  stopifnot(is.matrix(lambda3))
  nlam3 <- nrow(lambda3)
  
  # output from step 1
  step1 <- list()
  step1$mu <- mu
  step1$beta <- nonintercepts
  step1$a0 <- intercepts
  
  # pre-specify the returns
  step2 <- vector("list", nlam1)
  step3 <- vector("list", nlam1)
  
  for (k in seq(nlam1)){
    # The Second Step:
    # find num_keep higher order terms from
    # (1) squared main effects M^2 + Interaction effects I
    #     (square == FALSE)
    # (2) Interaction effects I (square == TRUE)
    offset <- matrix(step1$mu[,k], nrow=n, byrow=TRUE)
    if(inherits(x, "sparseMatrix")){
      idx <- screen_sparse_cpp(x = xm, y = y, yMat = yMat, nlev = nlev, offset = offset, num_keep = num_keep, square = square, family = "ordinal")
    }else{
      idx <- screen_cpp(x = xm, y = y, nlev = nlev, offset = offset, num_keep = num_keep, square = square, family = "ordinal")
    }
    
    # remove NA if there is any
    idx <- idx[!is.na(idx[, 3]), ]
    
    # preparing for Step 3
    # idx has two columns
    # which are the j,k indices of nonzero elements
    # main effect index is of form (0, k)
    # interaction effect index is of form (j, k) for j < k
    idx <- idx[order(idx[, 3], decreasing = TRUE), , drop = FALSE]
    colnames(idx) <- c("index_1", "index_2", "score")
    step2[[k]] <- idx
    
    # construct design matrix of selected interactions
    if(nrow(idx) == 1){
      design <- myscale(matrix(xm[, idx[, 1]] * xm[, idx[, 2]], ncol = 1))
    }else{
      design <- myscale(xm[, idx[, 1]] * xm[, idx[, 2]])
    }
    
    # grab scales 
    col_mean3 <- c(col_mean,
                   attr(design, which = "scaled:center"))
    col_sd3 <- c(col_sd,
                 attr(design, which = "scaled:scale"))
    
    # the total design matrix
    design <- cbind(x, design)
    
    # construct xMat without intercept columns
    xList <- lapply(1:nrow(design), function(i)
    {
      xi <- design[i, ]
      xListi <- rbind(xi)[rep(1, nlev-1), , drop=FALSE]
      rownames(xListi) <- NULL
      xListi
    })
    xMat <- if(inherits(x, "sparseMatrix")) Matrix::Matrix(do.call(rbind, xList), sparse = TRUE) else do.call(rbind, xList)
    
    # The Third Step:
    # in Step3, we fit the response with offset being mu from Step 1
    if(any(is.na(lambda3[, k]))){
      lambda3[, k] <- get_lambda(x = design, y = y, offset = offset, family = "ordinal", intercept = FALSE, nlam = nlam3, lam_min_ratio = lam_min_ratio)
    }
    
    # in Step 3 we no longer fit intercept, which has been captured in Step 1
    fit <-  ordinet(x = design, y = y, offset = offset, lambdaVals = lambda3[, k], intercept = FALSE, standardize = FALSE)
    
    # grab coefficient estimate
    coef <- t(fit$coefs)
    colnames(coef) <- NULL
    rownames(coef) <- NULL
    
    # scale estimates back to the original scale of design
    nonintercepts0 <- coef
    unscaleFact <- col_mean3 / col_sd3
    intAdjust <- (unscaleFact %*% nonintercepts0)[rep(1, nlev-1), , drop=FALSE]
    intercepts <- - intAdjust 
    rownames(intercepts) <- paste0("Intercept:", 1:(nlev-1))
    nonintercepts <- nonintercepts0 / col_sd3
    step3[[k]]$coef <- nonintercepts
    step3[[k]]$a0 <- intercepts 
    
    # number of non-zero main effects & interactions
    step3[[k]]$nzm <- colSums(as.matrix(coef[1:p, ] != 0))
    step3[[k]]$nzi <- colSums(as.matrix(coef[(p + 1): ncol(coef), ] != 0))
  }
  
  
  result <- list(n = n, p = p, square = square,
                 family = family, nlev = nlev,
                 cv_step1 = cv_step1,
                 step1 = step1, lambda1 = lambda1,
                 step2 = step2, num_keep = num_keep,
                 step3 = step3, lambda3 = lambda3,
                 call = match.call())
  return(result)
}
