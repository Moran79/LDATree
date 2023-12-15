
# Main function -----------------------------------------------------------


getSplitFun <- function(datX, response, method, modelLDA){
  if(method == "FACT") return(getSplitFunFACT(datX = datX,
                                              response = response,
                                              modelLDA = modelLDA))
  else if(method == "groupMean") return(getSplitFunGroupMean(datX = datX,
                                                       response = response,
                                                       modelLDA = modelLDA))
  else if(method == "mixed") return(getSplitFunMixed(datX = datX,
                                                     response = response,
                                                     modelLDA = modelLDA))
}


# FACT --------------------------------------------------------------------


getSplitFunFACT <- function(datX, response, modelLDA){
  #> This function is called only when building the tree
  predictedOutcome <- predict(modelLDA, datX)
  # browser()
  if(length(unique(predictedOutcome)) == 1) return(NULL)

  #> If there are some classes not being predicted
  #> we will assign them to the class with the second largest posterior prob
  idxPred <- which(names(modelLDA$prior) %in% predictedOutcome)
  splitResInTraining <- lapply(seq_along(idxPred), function(i) which(names(modelLDA$prior)[i] == predictedOutcome))

  res <- function(datX, missingReference){
    fixedData <- getDataInShape(data = datX, missingReference = missingReference)
    predictedProb <- predict(modelLDA, fixedData,type = "prob")[,idxPred, drop = FALSE]
    predictedOutcome <- max.col(predictedProb, ties.method = "first")
    return(lapply(seq_along(idxPred), function(i) which(i == predictedOutcome)))
  }

  attr(res, "splitResInTraining") <- splitResInTraining # record the split function's form
  return(res)
}


# mixed -------------------------------------------------------------------


getSplitFunMixed <- function(datX, response, modelLDA){
  if(modelLDA$pValue< 5e-4) return(getSplitFunFACT(datX = datX,
                                                  response = response,
                                                  modelLDA = modelLDA)) #
  else return(getSplitFunGroupMean(datX = datX,
                                   response = response,
                                   modelLDA = modelLDA))
}


# linear regression line of the group means -------------------------------


getSplitFunGroupMean <- function(datX, response, modelLDA){
  #> This function is called only when building the tree
  #> in order to make it easier to write and debug, we force stepLDA
  #> to select at least J variables.

  modelFrame <- model.frame(formula = ~.-1, data = datX) # get all columns without intercept
  Terms <- terms(modelFrame)
  m <- model.matrix(Terms, modelFrame)
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))

  currentCandidates <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x)))))
  rankVar <- stepVarSelByFsmall(m = scale(m),
                                response = response,
                                currentCandidates = currentCandidates,
                                k = min(ncol(m), nlevels(response)))
  numOfPredictors <-  length(rankVar) - 1

  # if(nrow(datX) == 528) browser()

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

  res <- function(datX, missingReference){
    fixedData <- getDataInShape(data = datX, missingReference = missingReference)
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


# New split checking ------------------------------------------------------


generateSplitNchildren <- function(datX,
                                   response,
                                   idxCol,
                                   idxRow,
                                   treeType,
                                   ldaType,
                                   fastTree,
                                   missingMethod,
                                   splitMethod,
                                   minNodeSize,
                                   bootstrap = TRUE){
  #> This function can be used to either:
  #> 1. evaluate the split on the bootstrap data, or
  #> 2. distribute the split using all data [Not yet added] [refresh current level and parentIndex]

  treeList = structure(list(), class = "SingleTreee") # save the tree
  #> Bootstrap the data
  idxRowTrain <- sample(idxRow, replace = TRUE)
  idxRowOOB <- setdiff(idxRow, idxRowTrain)
  #> Exit1: stop the split because there is not enough data
  if(length(idxRowOOB) < 5) return(1)

  #> Generate the first node
  treeList[[1]] <- new_TreeeNode(datX = datX,
                                 response = response,
                                 idxCol = idxCol,
                                 idxRow = idxRowTrain,
                                 treeType = treeType,
                                 ldaType = ldaType,
                                 fastTree = fastTree,
                                 missingMethod = missingMethod,
                                 splitMethod = splitMethod,
                                 maxTreeLevel = 1,
                                 minNodeSize = minNodeSize,
                                 currentLevel = 0,
                                 parentIndex = 0)

  #> Exit2: stop the split because the split is not found
  if(treeList[[1]]$stopFlag != 0) return(1)

  testResBefore <- predict(treeList, datX[idxRowOOB,, drop = FALSE]) == as.character(response[idxRowOOB])
  # cat("Before: ", mean(testResBefore), "\n")

  #> Distribute the training set
  trainIndex <- treeList[[1]]$splitFun(datX = datX[idxRowTrain,,drop = FALSE], missingReference = treeList[[1]]$misReference)

  #> Get child nodes
  childNodes <- lapply(seq_along(trainIndex), function(i) new_TreeeNode(datX = datX,
                                                                        response = response,
                                                                        idxCol = idxCol,
                                                                        idxRow = idxRowTrain[trainIndex[[i]]],
                                                                        treeType = treeType,
                                                                        ldaType = ldaType,
                                                                        fastTree = fastTree,
                                                                        missingMethod = missingMethod,
                                                                        splitMethod = splitMethod,
                                                                        maxTreeLevel = 0,
                                                                        minNodeSize = minNodeSize,
                                                                        currentLevel = 1,
                                                                        parentIndex = 1))
  #> Put nodes into the tree
  for(i in seq_along(childNodes)) treeList[[i + 1]] <- childNodes[[i]]
  treeList[[1]]$children <- 1 + seq_along(childNodes)
  testResAfter <- predict(treeList, datX[idxRowOOB,, drop = FALSE]) == as.character(response[idxRowOOB])
  # cat("After: ", mean(testResAfter), "\n")

  # tTestRes <- t.test(testResAfter, testResBefore, paired = TRUE, alternative = "greater")
  tTestRes <- tryCatch({t.test(testResAfter, testResBefore, paired = TRUE, alternative = "greater")},
                       error = function(e){list(p.value = pt(sqrt(length(testResAfter)), df = length(testResAfter) - 1, lower.tail = FALSE))})
  return(tTestRes$p.value)
}

