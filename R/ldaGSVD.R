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
ldaGSVD <- function(formula, data, method = "all", strict = TRUE){
  method = match.arg(method, c("step", "all"))
  modelFrame <- model.frame(formula, data, na.action = "na.fail")
  Terms <- terms(modelFrame)
  response <- droplevels(as.factor(modelFrame[,1])) # some levels are branched out
  prior <- table(response, dnn = NULL) / length(response) # estimated prior

  # Design Matrix
  m <- scale(model.matrix(formula, data))
  cnames <- colnames(m)
  currentVarList <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x))))) # remove constant columns and intercept

  if(method == "step"){
    stepRes <- stepVarSelByF(m = m, response = response, currentCandidates = currentVarList, strict = strict)
    currentVarList <- stepRes$currentVarList

    if(length(currentVarList) != 0){
      #> modify the design matrix to make it more compact
      #> so that only the selected variables are included in the design matrix
      selectedVarRawIdx <- unique(sort(attributes(m)$assign[currentVarList]))
      formula <- as.formula(paste(colnames(modelFrame)[1],"~", paste(colnames(modelFrame)[1+selectedVarRawIdx], collapse="+"), "-1"))
      modelFrame <- model.frame(formula, data, na.action = "na.fail")
      Terms <- terms(modelFrame)
      m <- scale(model.matrix(formula, data))
      cnames <- colnames(m)
      currentVarList <- which(cnames %in% stepRes$stepInfo$var)
    }else{ # When no variable is selected
      # warning("None of the variables is significant. The full model is fitted instead.")
      currentVarList <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x)))))
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

getLambda <- function(Sw, St){
  # return 1 if the St is linear correlated
  detMw <- det(Sw)
  detMt <- det(St)
  if(detMt == 0) return(1)
  rawAlpha <- detMw / detMt
  if(rawAlpha < 0) return(1)
  return(rawAlpha)
}

selectVar <- function(currentVar, newVar, Sw, St, direction = "forward"){
  #> return the column index
  #> return 0 if all var makes St = 0
  if(direction == "forward"){
    lambdaAll <- sapply(newVar, function(i) getLambda(Sw[c(i,currentVar),c(i,currentVar), drop = FALSE], St[c(i,currentVar),c(i,currentVar), drop = FALSE]))
  }else{
    lambdaAll <- sapply(currentVar, function(i) getLambda(Sw[setdiff(currentVar, i),setdiff(currentVar, i), drop = FALSE], St[setdiff(currentVar, i),setdiff(currentVar, i), drop = FALSE]))
  }
  currentVarIdx <- which.min(lambdaAll)
  return(list(stopflag = (lambdaAll[currentVarIdx] == 1),
              varIdx = newVar[currentVarIdx],
              statistics = max(lambdaAll[currentVarIdx], 1e-10)))
}

getCheckIdx <- function(currentVarList, currentCandidates){
  # When dynamically update the Sw and St, this function is used
  # to return all used elements in those matrices
  if(length(currentVarList) == 0) return(data.frame(Var1 = currentCandidates, Var2 = currentCandidates))
  else return(expand.grid(currentVarList, currentCandidates))
}

fillNAinS <- function(idxCheck, m, mW, envir = parent.frame()){
  #> Warning: To avoid having a copy of the p*p scatter matrices,
  #> this function use `assign` to modify the Sw and St without explicitly returning
  Sw <- get("Sw", envir = envir)
  St <- get("St", envir = envir)
  for(i in seq_len(nrow(idxCheck))){
    if(is.na(Sw[idxCheck$Var1[i], idxCheck$Var2[i]])){
      Sw[idxCheck$Var1[i], idxCheck$Var2[i]] <- Sw[idxCheck$Var2[i], idxCheck$Var1[i]] <- sum(mW[,idxCheck$Var1[i]] * mW[,idxCheck$Var2[i]]) / (nrow(m))
      St[idxCheck$Var1[i], idxCheck$Var2[i]] <- St[idxCheck$Var2[i], idxCheck$Var1[i]] <- sum(m[,idxCheck$Var1[i]] * m[,idxCheck$Var2[i]]) / (nrow(m))
    }
  }
  assign("Sw", Sw, envir = envir)
  assign("St", St, envir = envir)
  return(NULL)
}

stepVarSelByF <- function(m, response, currentCandidates, strict = TRUE){
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response,]

  # Initialize
  n = nrow(m); g = nlevels(response); p = 0; currentVarList = c();
  currentLambda <- 1; kRes <- 1; Sw <- St <- matrix(NA, nrow = ncol(m), ncol = ncol(m))
  stepInfo <- data.frame(var = character(2*ncol(m)),
                         FtoEnter = 0,
                         FtoRemove = 0,
                         pValue = 0)

  # Stepwise selection starts!
  while(length(currentCandidates) != 0){
    nCandidates <- length(currentCandidates)
    p = p + 1
    if(p >= n - g + 1) break # F-statistic can not be calculated

    # Update the Sw and St if needed
    idxCheck <- getCheckIdx(currentVarList = currentVarList, currentCandidates = currentCandidates)
    fillNAinS(idxCheck = idxCheck, m = m, mW = mW)
    selectVarInfo <- selectVar(currentVar = currentVarList,
                               newVar = currentCandidates,
                               Sw = Sw,
                               St = St)
    if(selectVarInfo$stopflag) break # If St = 0, stop

    # Get the F-to-enter
    partialLambda <- selectVarInfo$statistics / currentLambda
    currentLambda <- selectVarInfo$statistics
    FtoEnter <- (n - g - p + 1) / (g - 1) * (1 - partialLambda) / partialLambda
    pValue <- pf(FtoEnter,df1 = g - 1, df2 = n - g - p + 1, lower.tail = FALSE) * nCandidates
    if(FtoEnter <= 4) break # If no significant variable selected, stop
    if(strict & pValue > 1) break # Bonferroni's corrected p, stop

    # Add the variable into the model
    currentVarList <- c(currentVarList, selectVarInfo$varIdx)
    currentCandidates <- setdiff(currentCandidates, selectVarInfo$varIdx)
    stepInfo$var[kRes] <- colnames(m)[selectVarInfo$varIdx]
    stepInfo$FtoEnter[kRes] <- FtoEnter
    stepInfo$pValue[kRes] <- pValue
    kRes <- kRes + 1

    # Get the F-to-remove
    if(length(currentVarList)>1){
      selectVarInfoOut <- selectVar(currentVar = currentVarList,
                                    newVar = currentVarList,
                                    Sw = Sw,
                                    St = St,
                                    direction = "backward")
      partialLambdaOut <- currentLambda / selectVarInfoOut$statistics
      FtoRemove <- (n - g - p + 1) / (g - 1) * (1 - partialLambdaOut) / partialLambdaOut

      # Remove the variable from the model
      if(FtoRemove <= 3.996){
        currentVarList <- setdiff(currentVarList, selectVarInfoOut$varIdx)
        # cat("Variable is removed:", selectVarInfoOut$varIdx, "\n")
        stepInfo$var[kRes] <- colnames(m)[selectVarInfoOut$varIdx]
        stepInfo$FtoRemove[kRes] <- FtoRemove
        kRes <- kRes + 1
      }
    }
  }

  # Remove the empty rows in the stepInfo if stepLDA does not select all variables
  if(any(stepInfo$var == "")) stepInfo <- stepInfo[stepInfo$var != "",]

  return(list(currentVarList = currentVarList, stepInfo = stepInfo))
}
