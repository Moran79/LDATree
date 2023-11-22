
# Main function -----------------------------------------------------------

getSplitFun <- function(x, response, method, modelLDA){
  if(method == "FACT") return(getSplitFunFACT(x = x,
                                              response = response,
                                              modelLDA = modelLDA))
  else if(method == "groupMean") return(getSplitFunGroupMean(x = x,
                                                       response = response,
                                                       modelLDA = modelLDA))
  else if(method == "mixed") return(getSplitFunMixed(x = x,
                                                     response = response,
                                                     modelLDA = modelLDA))
}


# FACT --------------------------------------------------------------------

getSplitFunFACT <- function(x, response, modelLDA){
  #> This function is called only when building the tree
  predictedOutcome <- predict(modelLDA, x)
  if(length(unique(predictedOutcome)) == 1) return(NULL)

  #> If there are some classes not being predicted
  #> we will assign them to the class with the second largest posterior prob
  idxPred <- which(names(modelLDA$prior) %in% predictedOutcome)

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
    predictedProb <- predict(modelLDA, fixedData,type = "prob")[,idxPred, drop = FALSE]
    predictedOutcome <- max.col(predictedProb, ties.method = "first")
    finalList <- lapply(seq_along(idxPred), function(i) which(i == predictedOutcome))
    # if(any(sapply(finalList, length) == 0)) browser()
    return(lapply(seq_along(idxPred), function(i) which(i == predictedOutcome)))
  }
}


# mixed -------------------------------------------------------------------

getSplitFunMixed <- function(x, response, modelLDA){
  if(modelLDA$pValue<5e-4) return(getSplitFunFACT(x = x,
                                                  response = response,
                                                  modelLDA = modelLDA))
  else return(getSplitFunGroupMean(x = x,
                                   response = response,
                                   modelLDA = modelLDA))
}

# linear regression line of the group means -------------------------------

getSplitFunGroupMean <- function(x, response, modelLDA){
  #> This function is called only when building the tree
  #> Fixed version

  modelFrame <- model.frame(formula = ~.-1, data = x) # get all columns without intercept
  Terms <- terms(modelFrame)
  m <- model.matrix(Terms, modelFrame)
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))

  currentCandidates <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x)))))
  rankVar <- stepVarSelByFsmall(m = scale(m),
                                response = response,
                                currentCandidates = currentCandidates,
                                k = min(ncol(m), nlevels(response)))
  numOfPredictors <-  length(rankVar) - 1

  lmCoef <- NULL
  while(is.null(lmCoef) & numOfPredictors > 0){
    Y <- groupMeans[, rankVar[1], drop = FALSE] # use the most important variable as the y
    X <- cbind(numeric(nrow(groupMeans)) + 1, groupMeans[, rankVar[seq_len(numOfPredictors)+1], drop = FALSE])
    lmCoef <- tryCatch({
      solve(t(X) %*% X) %*% t(X) %*% Y
    }, error = function(e) {NULL})
    if(is.null(lmCoef)) numOfPredictors <- numOfPredictors - 1 # remove one variable if not invertible
  }
  if(is.null(lmCoef)) return(NULL) # no split


  # generate the final formula, with intercept
  splitCoef <- numeric(ncol(m) + 1)
  splitCoef[c(1, rankVar[seq_len(numOfPredictors)+1]+1)] <- lmCoef
  splitCoef[rankVar[1]+1] <- -1

  attr(splitCoef, "dmInfo") <- list(terms = Terms, xlevels = .getXlevels(Terms, modelFrame))

  # Stop the split if all points belong to one side
  projectionOnSplit <- unname(as.vector(cbind(1,m) %*% splitCoef))
  currentList <- list(which(projectionOnSplit < 0), which(projectionOnSplit >= 0))
  if(any(sapply(currentList, length) == 0)) return(NULL) # all points go to one side

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
    m <- cbind(1, getDesignMatrix(modelLDA = attr(splitCoef, "dmInfo"), data = fixedData, scale = FALSE))
    projectionOnSplit <- unname(as.vector(m %*% splitCoef))
    # return the relative index
    return(list(which(projectionOnSplit < 0), which(projectionOnSplit >= 0)))
  }
  attr(res, "splitCoef") <- splitCoef # record the split function's form
  return(res)
}

stepVarSelByFsmall <- function(m, response, currentCandidates, k){
  idxOriginal <- currentCandidates
  m <- m[,currentCandidates, drop = FALSE] # all columns should be useful

  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response, , drop = FALSE]

  # Initialize
  n = nrow(m); g = nlevels(response); p = 0; currentVarList = c()
  previousPillai <- previousDiff <- previousDiffDiff <- numeric(ncol(m)+1);
  previousDiff[1] <- Inf; diffChecker <- 0
  kRes <- 1; currentCandidates <- seq_len(ncol(m))
  Sw <- St <- matrix(NA, nrow = ncol(m), ncol = ncol(m))
  diag(Sw) <- apply(mW^2,2,sum); diag(St) <- apply(m^2,2,sum)

  # # Empirical: Calculate the threshold for pillaiToEnter
  # pillaiThreshold <- 1 / (1 + (n-g) / (abs(g-2)+1) / qf(1 - 0.1 / (ncol(m)+1), abs(g-2)+1, n-g)) / currentCandidates^(0.25)

  stepInfo <- data.frame(var = character(2*ncol(m)),
                         pillaiToEnter = 0,
                         pillaiToRemove = 0,
                         pillai = 0)

  stopFlag <- 0
  # Stepwise selection starts!
  while(length(currentVarList) < k & length(currentCandidates) > 0){

    nCandidates <- length(currentCandidates)
    p = p + 1

    # selectVarInfo <- tryCatch({
    #   selectVar(currentVar = currentVarList,
    #             newVar = currentCandidates,
    #             Sw = Sw,
    #             St = St)
    # }, error = function(e) {browser()})

    selectVarInfo <- selectVar(currentVar = currentVarList,
                               newVar = currentCandidates,
                               Sw = Sw,
                               St = St)
    bestVar <- selectVarInfo$varIdx
    if(selectVarInfo$stopflag){ # If St = 0, stop. [Might never happens, since there are other variables to choose]
      stopFlag <- 1
      break
    }

    # get the difference in Pillai's trace
    previousDiff[p+1] <- selectVarInfo$statistics - previousPillai[p]
    previousDiffDiff[p+1] <- previousDiff[p+1] - previousDiff[p]
    diffChecker <- ifelse(abs(previousDiffDiff[p+1]) < 0.001, diffChecker + 1, 0)

    if(previousDiff[p+1] > 10 * previousDiff[p]){ # Correlated variable(s) is included
      currentCandidates <- setdiff(currentCandidates, selectVarInfo$maxVarIdx)
      p <- p - 1; next
    }

    # Add the variable into the model
    previousPillai[p+1] <- selectVarInfo$statistics
    currentVarList <- c(currentVarList, bestVar)
    currentCandidates <- setdiff(currentCandidates, bestVar)
    stepInfo$var[kRes] <- colnames(m)[bestVar]
    stepInfo$pillaiToEnter[kRes] <- previousDiff[p+1]
    stepInfo$pillai[kRes] <- previousPillai[p+1]
    kRes <- kRes + 1

    # Update the Sw and St on the new added column
    Sw[currentCandidates, bestVar] <- Sw[bestVar, currentCandidates] <- as.vector(t(mW[, currentCandidates, drop = FALSE]) %*% mW[,bestVar, drop = FALSE])
    St[currentCandidates, bestVar] <- St[bestVar, currentCandidates] <- as.vector(t(m[, currentCandidates, drop = FALSE]) %*% m[,bestVar, drop = FALSE])
  }

  # Remove the empty rows in the stepInfo if stepLDA does not select all variables
  stepInfo <- stepInfo[seq_along(currentVarList),]

  # why return bestVar: in case no variable is significant, use this
  return(idxOriginal[currentVarList])
}
