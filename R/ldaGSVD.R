#' Linear Discriminant Analysis using the Generalized Singular Value
#' Decomposition
#'
#' @description `r lifecycle::badge('experimental')` Fit an LDA/GSVD model.
#'
#' @details
#'
#' Traditional Fisher's Linear Discriminant Analysis (LDA) ceases to work when
#' the within-class scatter matrix is singular. The Generalized Singular Value
#' Decomposition (GSVD) is used to address this issue. GSVD simultaneously
#' diagonalizes both the within-class and between-class scatter matrices without
#' the need to invert a singular matrix. This method is believed to be more
#' accurate than PCA-LDA (as in `MASS::lda`) because it also considers the
#' information in the between-class scatter matrix.
#'
#' @param data a data frame that contains both predictors and the response.
#'   Missing values are NOT allowed.
#' @param formula 123
#' @param method default to be all
#' @param strict 123
#'
#' @returns An object of class `ldaGSVD` containing the following components:
#' * `scaling`: a matrix which transforms the training data to LD scores, normalized so that the within-group scatter matrix is proportional to the identity matrix.
#' * `formula`: the formula passed to the [ldaGSVD()]
#' * `terms`: a object of class `terms` derived using the input `formula` and the training data
#' * `prior`: a `table` of the estimated prior probabilities.
#' * `groupMeans`: a matrix that records the group means of the training data on the transformed LD scores.
#' * `xlevels`: a list records the levels of the factor predictors, derived using the input `formula` and the training data
#'
#' @export
#'
#' @references Ye, J., Janardan, R., Park, C. H., & Park, H. (2004). \emph{An
#'   optimization criterion for generalized discriminant analysis on
#'   undersampled problems}. IEEE Transactions on Pattern Analysis and Machine
#'   Intelligence
#'
#'   Howland, P., Jeon, M., & Park, H. (2003). \emph{Structure preserving dimension
#'   reduction for clustered text data based on the generalized singular value
#'   decomposition}. SIAM Journal on Matrix Analysis and Applications
#'
#' @examples
#' fit <- ldaGSVD(Species~., data = iris)
#' # prediction
#' predict(fit,iris)
ldaGSVD <- function(formula, data, method = "all", strict = 10, stepTimeCapInMins = 20, ...){
  print(strict)
  method = match.arg(method, c("step", "all"))
  modelFrame <- model.frame(formula, data, na.action = "na.fail")
  Terms <- terms(modelFrame)
  response <- droplevels(as.factor(modelFrame[,1])) # some levels are branched out
  prior <- table(response, dnn = NULL) / length(response) # estimated prior

  # Design Matrix
  m <- scale(model.matrix(formula, data)) # constant cols would be changed to NaN in this step
  cnames <- colnames(m)
  currentVarList <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x))))) # remove constant columns and intercept

  if(method == "step"){
    #> Output: currentVarList, which contains indices of the selected variables
    #> RESPECTIVELY in the design matrix, some columns of m might be removed
    stepRes <- stepVarSelByF(m = m, response = response, currentCandidates = currentVarList,
                             strict = strict, stepTimeCapInMins = stepTimeCapInMins)
    currentVarList <- stepRes$currentVarList

    if(length(currentVarList) != 0){
      #> modify the design matrix and formula to make it more compact
      #> so that only the selected variables are included in the design matrix,
      #> and eventually make the prediction faster
      selectedVarRawIdx <- unique(sort(attributes(m)$assign[currentVarList])) # MUST be from the modelFrame where the factors are not dummied
      formula <- as.formula(paste(colnames(modelFrame)[1],"~", paste(colnames(modelFrame)[1+selectedVarRawIdx], collapse="+")))
      modelFrame <- model.frame(formula, data, na.action = "na.fail")
      Terms <- terms(modelFrame)
      m <- scale(model.matrix(formula, data)) # This double scaling is not optimal,
      # but prevent losing all the attributes due to subseting
      selectedVarName <- setdiff(stepRes$stepInfo$var, stepRes$stepInfo$var[stepRes$stepInfo$FtoRemove != 0]) # all vars that are not being removed
      currentVarList <- which(colnames(m) %in% selectedVarName)
    }else{
      #> When no variable is selected, use only the best single variable
      #> instead of using the full model to save some time
      currentVarList <- stepRes$bestVar
    }
  }

  varSD <- attr(m,"scaled:scale")[currentVarList]
  varCenter <- attr(m,"scaled:center")[currentVarList]
  m <- m[,currentVarList, drop = FALSE]

  # Step 1: SVD on the combined matrix H
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  Hb <- sqrt(tabulate(response)) * groupMeans # grandMean = 0 if scaled
  fitSVD <- svd(rbind(Hb, m - groupMeans[response,]))
  rankT <- sum(fitSVD$d >= max(dim(fitSVD$u),dim(fitSVD$v)) * .Machine$double.eps * fitSVD$d[1])

  # Step 2: SVD on the P matrix
  fitSVDp <- svd(fitSVD$u[seq_len(nlevels(response)), seq_len(rankT), drop = FALSE], nu = 0L)
  rankAll <- min(nlevels(response)-1, rankT) # This is not optimal, but rank(Hb) takes time
  # Fix the variance part
  unitSD <- pmin(diag(sqrt((length(response) - nlevels(response)) / abs(1 - fitSVDp$d^2 + 1e-15)), nrow = rankAll),1e15) # Scale to unit var
  scalingFinal <- (fitSVD$v[,seq_len(rankT), drop = FALSE] %*% diag(1 / fitSVD$d[seq_len(rankT)], nrow = rankT) %*% fitSVDp$v)[,seq_len(rankAll), drop = FALSE] %*% unitSD
  rownames(scalingFinal) <- cnames[currentVarList]

  ### TESTING ###
  mNew <- m %*% scalingFinal
  groupMeansNew <- tapply(c(mNew), list(rep(response, dim(mNew)[2]), col(mNew)), function(x) mean(x, na.rm = TRUE))
  HbNew <- sqrt(tabulate(response)) * groupMeansNew # grandMean = 0 if scaled
  Sb <- t(HbNew) %*% HbNew
  HwNew <- mNew - groupMeansNew[response,]
  Sw <- t(HwNew) %*% HwNew
  print(sum(diag(solve(Sb + Sw, Sb))))
  ###############

  groupMeans <- groupMeans %*% scalingFinal
  rownames(groupMeans) <- levels(response)
  colnames(groupMeans) <- colnames(scalingFinal) <- paste("LD", seq_len(ncol(groupMeans)), sep = "")

  res <- list(scaling = scalingFinal, formula = formula, terms = Terms, prior = prior,
              groupMeans = groupMeans, xlevels = .getXlevels(Terms, modelFrame),
              varIdx = currentVarList, varSD = varSD, varCenter = varCenter)
  if(method == "step"){
    res$stepInfo = stepRes$stepInfo
  }
  class(res) <- "ldaGSVD"
  return(res)
}

#' Predictions from a fitted ldaGSVD object
#'
#' Prediction of test data using a fitted ldaGSVD object
#'
#' Unlike the original paper, which uses the k-nearest neighbor (k-NN) as the
#' classifier, we use a faster and more straightforward likelihood-based method.
#' One limitation of the traditional likelihood-based method for LDA is that it
#' ceases to work when there are Linear Discriminant (LD) directions with zero
#' variance in the within-class scatter matrix. However, when using LDA/GSVD,
#' all chosen LD directions possess non-zero variance in the between-class
#' scatter matrix. This implies that LD directions with zero variance in the
#' within-class scatter matrix will yield the highest Fisher's ratio. Therefore,
#' to get these directions higher weights, we manually adjust the zero variance
#' to `1e-15` for computational reasons.
#'
#' @param object a fitted model object of class `ldaGSVD`, which is assumed to
#'   be the result of the [ldaGSVD()] function.
#' @param newdata data frame containing the values at which predictions are
#'   required. Missing values are NOT allowed.
#' @param type character string denoting the type of predicted value returned.
#'   The default is to return the predicted class (`type` = 'response'). The
#'   predicted posterior probabilities for each class will be returned if `type`
#'   = 'prob'.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function returns different values based on the `type`, if
#' * `type = 'response'`: vector of predicted responses.
#' * `type = 'prob'`: a data frame of the posterior probabilities. Each class takes a column.
#' @export
#'
#' @references Ye, J., Janardan, R., Park, C. H., & Park, H. (2004). \emph{An
#'   optimization criterion for generalized discriminant analysis on
#'   undersampled problems}. IEEE Transactions on Pattern Analysis and Machine
#'   Intelligence
#'
#'   Howland, P., Jeon, M., & Park, H. (2003). \emph{Structure preserving dimension
#'   reduction for clustered text data based on the generalized singular value
#'   decomposition}. SIAM Journal on Matrix Analysis and Applications
#'
#' @examples
#' fit <- ldaGSVD(Species~., data = iris)
#' predict(fit,iris)
#' # output prosterior probabilities
#' predict(fit,iris,type = "prob")
predict.ldaGSVD <- function(object, newdata, type = c("response", "prob"), ...){
  type <- match.arg(type, c("response", "prob"))
  # add one extra check for levels of the predictors
  LDscores <- getLDscores(modelLDA = object, data = newdata)
  loglikelihood <- LDscores %*% t(object$groupMeans) + matrix(log(object$prior) - 0.5 * rowSums(object$groupMeans^2), nrow(LDscores), length(object$prior), byrow = TRUE)
  # Computation Optimization 2: Prevent a very large likelihood due to exponential
  likelihood <- exp(loglikelihood - apply(loglikelihood, 1, max))
  posterior <- likelihood / apply(likelihood, 1, sum)
  if(type == "prob") return(posterior)
  return(rownames(object$groupMeans)[max.col(posterior, ties.method = "first")])
}


#' @export
print.ldaGSVD <- function(x, ...){
  cat("\nObserved proportions of groups:\n")
  print(x$prior)
  cat("\n\nGroup means of LD scores:\n")
  print(x$groupMeans)
  cat("\n\nScaling (coefficients) of LD scores:\n")
  print(x$scaling)
  invisible(x)
}


# helper functions --------------------------------------------------------

getPillai <- function(Sw, St){
  #> return -1 if the St is not full rank
  #> sometimes diag has negative value due to instability of R
  #> We ignore those for now
  tryCatch({
    sum(diag(solve(St, St - Sw)))
  }, error = function(e) {-1})
}

selectVar <- function(currentVar, newVar, Sw, St, direction = "forward"){
  #> return the column index
  #> return 0 if all var makes St = 0
  if(direction == "forward"){
    lambdaAll <- sapply(newVar, function(i) getPillai(Sw[c(i,currentVar),c(i,currentVar), drop = FALSE], St[c(i,currentVar),c(i,currentVar), drop = FALSE]))
  }else{
    lambdaAll <- sapply(currentVar, function(i) getPillai(Sw[setdiff(currentVar, i),setdiff(currentVar, i), drop = FALSE], St[setdiff(currentVar, i),setdiff(currentVar, i), drop = FALSE]))
  }
  currentVarIdx <- which.max(lambdaAll)

  return(list(stopflag = (lambdaAll[currentVarIdx] == -1),
              varIdx = newVar[currentVarIdx],
              statistics = lambdaAll[currentVarIdx]))
}


stepVarSelByF <- function(m, response, currentCandidates, strict = 10, stepTimeCapInMins = 20){
  idxOriginal <- currentCandidates
  m <- m[,currentCandidates, drop = FALSE] # all volumns should be useful

  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response,]

  # Initialize
  n = nrow(m); g = nlevels(response); p = 0; currentVarList = c()
  currentPillai <- 0; kRes <- 1; currentCandidates <- seq_len(ncol(m))
  Sw <- St <- matrix(NA, nrow = ncol(m), ncol = ncol(m))
  diag(Sw) <- apply(mW^2,2,sum); diag(St) <- apply(m^2,2,sum)

  # Empirical: Calculate the threshold for pillaiToEnter
  pillaiThreshold <- 1 / (1 + (n-g) / (abs(J-2)+1) / qf(1 - 0.1 / (ncol(m)+1), abs(J-2)+1, n-g)) / currentCandidates^(0.3)
  # print(pillaiThreshold)

  stepInfo <- data.frame(var = character(2*ncol(m)),
                         pillaiToEnter = 0,
                         pillaiToRemove = 0)

  timeOld <- Sys.time()

  # Stepwise selection starts!
  while(length(currentCandidates) != 0){
    nCandidates <- length(currentCandidates)
    p = p + 1
    if(p>strict) break

    timeNew <- Sys.time()
    if(difftime(timeNew, timeOld, units = "mins") > stepTimeCapInMins) break # when runtime is above some threshold

    selectVarInfo <- selectVar(currentVar = currentVarList,
                               newVar = currentCandidates,
                               Sw = Sw,
                               St = St)
    bestVar <- selectVarInfo$varIdx
    if(selectVarInfo$stopflag) break # If St = 0, stop

    # update the Pillai's trace
    pillaiDiff <- selectVarInfo$statistics - currentPillai
    # print(pillaiDiff)
    currentPillai <- selectVarInfo$statistics

    # Check the stopping rule
    # if(pillaiDiff < pillaiThreshold[p]) break # If no significant variable selected, stop
    if(pillaiDiff <= 0) break # If no significant variable selected, stop

    # Add the variable into the model
    currentVarList <- c(currentVarList, bestVar)
    currentCandidates <- setdiff(currentCandidates, bestVar)
    stepInfo$var[kRes] <- colnames(m)[bestVar]
    stepInfo$pillaiToEnter[kRes] <- pillaiDiff
    kRes <- kRes + 1

    # Update the Sw and St on the new added column
    Sw[currentCandidates, bestVar] <- Sw[bestVar, currentCandidates] <- as.vector(t(mW[, currentCandidates, drop = FALSE]) %*% mW[,bestVar, drop = FALSE])
    St[currentCandidates, bestVar] <- St[bestVar, currentCandidates] <- as.vector(t(m[, currentCandidates, drop = FALSE]) %*% m[,bestVar, drop = FALSE])

    # Removing process
    # if(length(currentVarList)>1){
    #   selectVarInfoOut <- selectVar(currentVar = currentVarList,
    #                                 newVar = currentVarList,
    #                                 Sw = Sw,
    #                                 St = St,
    #                                 direction = "backward")
    #
    #   pillaiDiffOut <- currentPillai - selectVarInfoOut$statistics
    #
    #   # Remove the variable from the model
    #   if(pillaiDiffOut < pillaiThreshold[p]){
    #     currentVarList <- setdiff(currentVarList, selectVarInfoOut$varIdx)
    #     # cat("Variable is removed:", selectVarInfoOut$varIdx, "\n")
    #     stepInfo$var[kRes] <- colnames(m)[selectVarInfoOut$varIdx]
    #     stepInfo$pillaiToRemove[kRes] <- pillaiDiffOut
    #     kRes <- kRes + 1
    #   }
    # }
  }

  # Remove the empty rows in the stepInfo if stepLDA does not select all variables
  if(any(stepInfo$var == "")) stepInfo <- stepInfo[stepInfo$var != "",]

  # why return bestVar: in case no variable is significant, use this
  return(list(currentVarList = idxOriginal[currentVarList], stepInfo = stepInfo, bestVar = idxOriginal[bestVar]))
}
