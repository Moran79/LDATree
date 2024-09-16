# constant In Group fix ---------------------------------------------------------

fixConstantGroup <- function(data, response, tol = 1e-8){
  #> Remove the -1 from the formula to make it more general,
  #> or the first level of the first categorical variable will be kept
  #> I also find sometimes it will throw constant error, so I think I just
  #> do a while loop for now, and improve the error in the future

  while(TRUE){

    m <- model.matrix.lm(~., data = data, na.action = "na.pass")

    response <- as.factor(response) # make sure it is a factor, since we are using indexing in tapply
    #> It is OK if the response has empty levels
    #> because groupMeans[response,] will not choose any information from the empty level
    groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
    sdInGroup <- apply(m - groupMeans[response,], 2, function(x) sd(x, na.rm = TRUE))
    idxConstM <- setdiff(which(sdInGroup < tol), 1) # remove the intercept, get the indices from design matrix
    idxConstRaw <- attr(m,"assign")[idxConstM] # from design matrix to raw data
    idxConstRawUniq <- unique(idxConstRaw) # get the unique columns. there will be multiple entries for factors

    if(length(idxConstRawUniq) ==0) {break}

    cat("Constant Group Program Starts!")
    cat("Affected Variable(s):", colnames(data)[idxConstRawUniq], "\n")
    #> add one observation is harder to code, and the response will be changed
    #> instead, we only modified observations

    fixList <- vector(mode = "list", length = length(idxConstRawUniq))
    colTypes <- getNumFlag(data[, idxConstRawUniq])
    for(i in seq_along(idxConstRawUniq)){
      colIdx <- idxConstRawUniq[i]
      fixList[[i]] <- list(colIdx = colIdx,
                           resolved = FALSE,
                           numOrNot = colTypes[i])
      if(!colTypes[i]){
        # Set up candidates for levels
        idxM <- idxConstM[which(idxConstRaw == colIdx)]
        fixList[[i]]$levels <- sub(colnames(data)[colIdx], "", colnames(m)[idxM])
      }
    }

    idxPluralityClass <- which(response == getMode(response))
    idxModified <- sample(idxPluralityClass, max(table(idxConstRaw)))
    for(o_0 in seq_along(idxModified)){
      currResult <- fixConstantGroupHelper(data = data, fixList = fixList, idxModified = idxModified[o_0])
      data <- currResult$data
      fixList <- currResult$fixList
    }
  }
  return(droplevels(data)) # levels with few obs are likely to be deleted
}

fixConstantGroupHelper <- function(data, fixList, idxModified){
  for(i in seq_along(fixList)){
    if(!fixList[[i]]$resolved){ # Only those still need fix
      currIdx <- fixList[[i]]$colIdx
      if(fixList[[i]]$numOrNot){ # for num
        valueCandidates <- unique(data[, currIdx])
        stopifnot(length(valueCandidates) > 1)
        newValue <- sample(setdiff(valueCandidates, data[idxModified, currIdx]), 1)
        data[idxModified, currIdx] <- newValue
        fixList[[i]]$resolved <- TRUE
      }else{ # for cat
        currLevel <- fixList[[i]]$levels[1]
        fixList[[i]]$levels <- fixList[[i]]$levels[-1] # delete the first levels
        data[, currIdx] <- as.factor(data[, currIdx])
        data[idxModified, currIdx] <- ifelse(data[idxModified, currIdx] == currLevel,
                                             levels(data[, currIdx])[1], currLevel)
        if(length(fixList[[i]]$levels) == 0) {fixList[[i]]$resolved <- TRUE}
      }
    }
  }
  return(list(data = data, fixList = fixList))
}

fixReferenceWithData <- function(data, ref){
  #> ref contains a wilder range of levels
  #> some of the levels might not exist in the data
  #> we delete those extra levels in ref
  #> But also, some columns in ref has already been deleted from data (constant)
  ref <- ref[,match(colnames(data), colnames(ref)), drop = FALSE]
  for(i in seq_along(ref)){
    if(!is.null(levels(ref[,i]))){ # categorical vars
      if(colnames(ref)[i] %in% colnames(data)){ # if the variable is not constant, so that it remains in data
        levelsData <- levels(data[, colnames(ref)[i]])
        stopifnot(!is.null(levelsData)) # check if it is a factor
        if(!setequal(levels(ref[,i]), levelsData)){
          ref[,i] <- factor(as.character(ref[,i]), levels = levelsData)
          if(is.na(ref[1,i])) {ref[1,i] <- getMode(data[, colnames(ref)[i]])}
        }
      }
    }
  }
  return(ref)
}


# Plot Indexing Helper: from BFT to LFT -----------------------------------


transIdForBiPlot <- function(treeeList){
  # Keep the tree structure, but change the indexing system
  res <- (seq_along(treeeList) == 1) + 0
  for(i in seq_along(treeeList)){
    treeeNode <- treeeList[[i]]
    terminalFlag <- is.null(treeeNode$children)
    # Index for binary tree
    if(!terminalFlag){
      currentIdx <- res[i]
      res[treeeNode$children] <- 2 * currentIdx + c(0,1)
    }
    # future: index for multi-split tree
  }
  return(res)
}

# Pruning -----------------------------------------------------------------

#> In the future, I might have to change the name of alpha to\
#> something like <drop in training error>

pruneByTrainErrDrop <- function(treeeList, pThreshold = 0.01, verbose = TRUE){
  stopifnot(class(treeeList) == "SingleTreee")
  lenBefore <- length(treeeList)
  treeeList <- makeAlphaMono(treeeList = treeeList)
  treeeList <- cutNodes(treeeList = treeeList, pThreshold = pThreshold)
  treeeList <- dropNodes(treeeList = treeeList)
  lenAfter <- length(treeeList)
  if(verbose) cat(paste(lenBefore - lenAfter, 'nodes are being pruned.\n'))
  return(treeeList)
}

makeAlphaMono <- function(treeeList){
  #> Purpose: make the alpha monotonic

  for(i in rev(seq_along(treeeList))){
    if(is.null(treeeList[[i]]$pruned)){ # Only loop over the unpruned nodes
      if(is.null(treeeList[[i]]$children)){
        treeeList[[i]]$alpha <- 1 # for terminal nodes
      } else{ # for intermediate nodes
        childrenAlpha <- sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$alpha)
        treeeList[[i]]$alpha <- min(c(1, treeeList[[i]]$alpha, unlist(childrenAlpha)), na.rm = TRUE)
      }
    }
  }
  return(treeeList)
}

getTerminalNodes <- function(currentIdx, treeeList, keepNonTerminal = FALSE){
  #> Get all terminal nodes that are offsprings from currentIdx
  treeeNode <- treeeList[[currentIdx]]
  if(is.null(treeeNode$children)){
    return(currentIdx)
  }else{
    nonTerminalNodes <- if(keepNonTerminal) currentIdx
    terminalNodes <- do.call(c, sapply(treeeNode$children,function(x) getTerminalNodes(currentIdx = x, treeeList = treeeList, keepNonTerminal = keepNonTerminal), simplify = FALSE))
    return(c(nonTerminalNodes, terminalNodes))
  }
}

cutNodes <- function(treeeList, pThreshold){
  for(i in rev(seq_along(treeeList))){
    treeeNode <- treeeList[[i]]
    # not yet pruned + non-terminal node + p-value above threshold
    cutFlag <- is.null(treeeNode$pruned) & !is.null(treeeNode$children) & treeeNode$alpha > pThreshold # R rounding error, 1e-10 needed
    if(cutFlag){
      allChildren <- getTerminalNodes(currentIdx = i, treeeList = treeeList, keepNonTerminal = TRUE)
      for(j in setdiff(allChildren, i)) {treeeList[[j]]$pruned <- TRUE}
      treeeList[[i]]["children"] <- list(NULL) # cut the branch
    }
  }
  # treeeList <- makeAlphaMono(treeeList = treeeList) # no need to update alpha, since it is not changed
  return(treeeList)
}


dropNodes <- function(treeeList){
  finalNodeIdx <- sort(getTerminalNodes(treeeList = treeeList, currentIdx = 1, keepNonTerminal = TRUE))
  treeeList <- treeeList[finalNodeIdx]
  class(treeeList) <- "SingleTreee"
  for(i in seq_along(treeeList)){
    treeeList[[i]]$currentIndex <- i # re-assign the currentIndex
    if(!is.null(treeeList[[i]]$children)) {
      treeeList[[i]]$children <- sapply(treeeList[[i]]$children, function(x) which(finalNodeIdx == x))
    }else{
      if(treeeList[[i]]$stopFlag == 0) treeeList[[i]]$stopFlag = 3 # due to pruning
    }
    treeeList[[i]]$parent <- which(finalNodeIdx == treeeList[[i]]$parent)
  }
  return(treeeList)
}

# Get tree depth ----------------------------------------------------------

getDepth <- function(treee){
  if(class(treee) == "Treee") treee = treee$treee
  depthAll <- numeric(length(treee))
  updateList <- seq_along(treee)[-1]
  if(length(updateList) != 0){
    for(i in updateList){
      currentNode <- treee[[i]]
      depthAll[i] <- depthAll[currentNode$parent] + 1
    }
  }
  return(depthAll)
}


# New Pruning (previous prune.R) ---------------------------------------------------

#> Usage
# treeList <- makeAlphaMono(treeeList = treeList, trainErrorCap = trainErrorCap)
# treeList <- pruneTreee(treeeList = treeList, trainErrorCap = trainErrorCap, alpha = -0.5)
# treeList <- dropNodes(treeList)

makeAlphaMono <- function(treeeList,
                          trainErrorCap){
  #> Purpose: Calculate alpha, and make it monotonic
  #> assumption: node index of the tree is larger than its children's indices

  for(i in rev(seq_along(treeeList))){
    # Only loop over the unpruned nodes
    if(is.null(treeeList[[i]]$pruned)){
      if(is.null(treeeList[[i]]$children)){
        # for terminal nodes
        treeeList[[i]]$alpha <- -Inf
      }else{
        # for intermediate nodes

        # current alpha
        numOfChildren <- length(treeeList[[i]]$children) #
        trainErrorCapNow <- ifelse(trainErrorCap == "zero", 0, numOfChildren - 1)

        currentAlpha <- (treeeList[[i]]$currentLoss - sum(sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$currentLoss)) - trainErrorCapNow) / (numOfChildren - 1)
        childrenAlpha <- sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$alpha)
        treeeList[[i]]$alpha <- max(c(-Inf, currentAlpha, unlist(childrenAlpha)), na.rm = TRUE)
      }
    }
  }
  return(treeeList)
}



# Splitting ---------------------------------------------------------------

### Sum of Pillai's trace ###

getSplitFunTrace <- function(datX, response, modelLDA){
  x <- getDesignMatrix(modelLDA = modelLDA, data = datX) # get the scaled x
  # if(any(attr(x, "scaled:scale") == 0)){
  #   x[, which(attr(x, "scaled:scale") == 0)] <- 0
  #   attr(x, "scaled:scale")[attr(x, "scaled:scale") == 0] <- 1
  # }

  M <- tapply(c(x), list(rep(response, dim(x)[2]), col(x)), function(o_o) mean(o_o, na.rm = TRUE))
  n <- nrow(x); nJ <- as.vector(table(response))

  getSumPillai <- function(a){ # a is the vector of interest
    # a[1] is the constant, a[-1] are coefficients of the variables
    group1Flag <- as.vector(x %*% a[-1] >= a[1])
    n1 <- sum(group1Flag)
    if(n1 %in% c(0, n)) return(Inf) # If the split is outside the data range

    # Calculate the smaller part, and use subtraction to get the other
    if(n1 > n / 2){
      group1Flag <- !group1Flag
      n1 <- sum(group1Flag)
    }
    n2 <- n - n1; idx <- which(group1Flag)
    nJ1 <- as.vector(table(response[group1Flag]))

    M1 <- tapply(c(x[idx, , drop = FALSE]), list(rep(response[idx], dim(x[idx, , drop = FALSE])[2]), col(x[idx, , drop = FALSE])), function(o_o) mean(o_o, na.rm = TRUE))
    M2 <- (M * nJ - M1 * nJ1) / (nJ - nJ1)
    Xmean1 <- apply(x[idx, , drop = FALSE],2,mean)
    Xmean2 <- - n1 * Xmean1 / n2
    Sb1 <- t(M1) - Xmean1
    Sb2 <- t(M2) - Xmean2
    -(n1 * sum(Sb1^2) + n2 * sum(Sb2^2))
  }

  optimRes <- optim(c(0, rep(1, ncol(x))), getSumPillai, control = list(maxit = 200), method = "SANN")
  aScaled <- optimRes$par
  # # scale center transformation
  aFinal <- aScaled

  # aFinal[-1] <- aFinal[-1] / attr(x, "scaled:scale")
  # aFinal[1] <- aFinal[1] + sum(aFinal[-1] * attr(x, "scaled:center"))

  # # Stop the split if all points belong to one side
  # projectionOnSplit <- unname(as.vector(x %*% matrix(aScaled[-1], ncol = 1)))
  # currentList <- list(which(projectionOnSplit < aScaled[1]), which(projectionOnSplit >= aScaled[1]))
  if(optimRes$value == Inf) return(NULL)

  res <- function(datX, missingReference){
    fixedData <- getDataInShape(data = datX, missingReference = missingReference)
    x <- getDesignMatrix(modelLDA = modelLDA, data = fixedData)
    projectionOnSplit <- unname(as.vector(x %*% matrix(aScaled[-1], ncol = 1)))
    # return the relative index
    return(list(which(projectionOnSplit < aScaled[1]), which(projectionOnSplit >= aScaled[1])))
  }
  attr(res, "a") <- aFinal # record the split function's form
  return(res)
}


### Gini Split ###

getSplitFunLD1 <- function(datX, response, modelLDA){

  LDscore <- getLDscores(modelLDA = modelLDA, data = datX, nScores = 1)
  idxRowOrdered <- order(LDscore)

  #> prevent empty nodes, so the lowest rank is removed
  #> max cut: 1000
  #> potentialCut: LDscores' ranks, a subset of 1 to length(response)
  percentageCut <- 0.05
  potentialCut <- unique(quantile(rank(LDscore,ties.method = "min"),
                                  probs = seq(percentageCut, 1 - percentageCut,length.out = 1000), type = 1))
  potentialCut <- setdiff(potentialCut,1)

  if(length(potentialCut)==0) {return(NULL)} # No cut due to ties

  # For the consideration of speed
  # increment programming is carried out
  GiniObserved <- numeric(length(potentialCut))
  NjtLeft <- table(response[idxRowOrdered][seq_len(potentialCut[1] - 1)])
  NjtRight <- table(response[idxRowOrdered][seq(potentialCut[1], length(response))])
  GiniObserved[1] <- sum(NjtLeft) * sum((NjtLeft / sum(NjtLeft))^2) +
    sum(NjtRight) * sum((NjtRight / sum(NjtRight))^2)
  for(i in seq_along(potentialCut)[-1]){
    responseTrans <- table(response[seq(potentialCut[i-1],potentialCut[i]-1)])
    NjtLeft <- NjtLeft + responseTrans
    NjtRight <- NjtRight - responseTrans
    GiniObserved[i] <- sum(NjtLeft) * sum((NjtLeft / sum(NjtLeft))^2) +
      sum(NjtRight) * sum((NjtRight / sum(NjtRight))^2)
  }
  cutPoint <- which.max(GiniObserved)
  cutScore <- LDscore[idxRowOrdered][[potentialCut[cutPoint]-1]]

  res <- function(datX, missingReference){
    fixedData <- getDataInShape(data = datX, missingReference = missingReference)
    LDscore <- unname(getLDscores(modelLDA = modelLDA, data = fixedData, nScores = 1))
    # return the relative index
    return(list(which(LDscore<=cutScore), which(LDscore>cutScore)))
  }
  return(res)
}


### FACT Split Old Version ###
#> 区别在于这个版本当predict的class少于样本中的class时，会去掉这个class重新fit
#> 新版本下不会去掉，而是在posterior prob上面做文章

# getSplitFunFACT <- function(x, response, modelLDA){
#   #> make sure that predictions have all the classes
#   #> otherwise, refit the modelLDA
#   predictedOutcome <- predict(modelLDA, x)
#   while(length(unique(predictedOutcome)) != length(modelLDA$prior)){
#     if(length(unique(predictedOutcome)) == 1) return(NULL) # if only one class is left
#     subsetIdx <- which(response %in% unique(predictedOutcome))
#     x <- x[subsetIdx,, drop = FALSE]; response <- response[subsetIdx]
#     datCombined = data.frame(response = response, x)
#     modelLDA <- ldaGSVD(response~., data = datCombined, method = "step")
#     predictedOutcome <- predict(modelLDA, x)
#   }
#
#   res <- function(x, missingReference){
#     fixedData <- getDataInShape(data = x, missingReference = missingReference)
#     predictedOutcome <- predict(modelLDA, fixedData)
#     return(lapply(names(modelLDA$prior), function(name) which(name == predictedOutcome)))
#   }
# }

# Group Mean splitting --- fast version #

getSplitFunGroupMean <- function(x, response, modelLDA){
  #> This function is called only when building the tree
  #> Fixed version

  modelFrame <- model.frame(formula = ~.-1, data = x) # get all columns without intercept
  Terms <- terms(modelFrame)
  m <- model.matrix(Terms, modelFrame)
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response, , drop = FALSE]
  # rank variables by their relative Pillai's trace in reverse order (very weird)
  rankVar <- order(apply(mW^2,2,sum) / apply(m,2,var), decreasing = FALSE)
  numOfPredictors <- min(ncol(m) - 1, nrow(groupMeans) - 1) # not include intercept


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


### linear regression line of the group means ###


getSplitFunGM <- function(datX, response, modelLDA){
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

### SPORF-like splitting ###

getSplitFunRandom <- function(datX, response, modelLDA){
  #> This is a easier version, just randomly select one direction

  modelFrame <- model.frame(formula = ~.-1, data = datX) # get all columns without intercept
  Terms <- terms(modelFrame)
  m <- model.matrix(Terms, modelFrame)
  sdM <- apply(m,2,sd) # get the SD for each column
  # k <- min(ncol(m), nlevels(response)) # the dimension of the transformation matrix
  currentCandidates <- which(sdM > 1e-10)
  if(length(currentCandidates) == 0) return(NULL) # no variable is available

  # Now we select sqrt(p) of the variables and form a random direction
  selectedColIdx <- sort(currentCandidates[sample(length(currentCandidates), ceiling(sqrt(length(currentCandidates))))])
  splitCoef <- numeric(ncol(m))
  splitCoef[selectedColIdx] <- sample(c(-1,1), length(selectedColIdx), replace = T) / sdM[selectedColIdx]
  splitCoef <- c(-median(m[, selectedColIdx, drop = FALSE] %*% splitCoef[selectedColIdx]), splitCoef)

  attr(splitCoef, "dmInfo") <- list(terms = Terms, xlevels = .getXlevels(Terms, modelFrame))

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


# Get x and response from formula and data --------------------------------

extractXnResponse <- function(formula, data){
  #> droplevels is necessary, since empty response level occurs during train/test split
  #> covariates can have empty levels as well
  modelFrame <- droplevels(model.frame(formula, data, na.action = "na.pass"))

  if(anyNA(modelFrame[,1])){
    modelFrame <- modelFrame[which(!is.na(modelFrame[,1])),, drop = FALSE]
    # no NAs are allowed in the response variable
    warning("There are missing values in the response variable. Related rows are removed")
  }

  response <- as.factor(modelFrame[,1])
  x <- modelFrame[,-1, drop = FALSE]
  return(list(x = x, response = response))
}


# stopping rule in Splits -------------------------------------------------


updateLDA <- function(oldLDA, xNew, responseNew, missingReference, ldaType){
  #> stop the update if 1.too few samples 2.only one y left
  if(nrow(xNew) < nlevels(responseNew) + 1 | length(unique(responseNew)) == 1) return(oldLDA)

  xNew <- getDataInShape(xNew, missingReference = missingReference)
  datCombined = data.frame(response = responseNew, xNew)
  datCombined = datCombined[,constantColCheck(datCombined), drop = FALSE]
  newLDA <- tryCatch({ldaGSVD(formula = response~., data = datCombined, method = ldaType, varName = rownames(oldLDA$scaling))},
                     error = function(e){oldLDA})
  return(newLDA)
}

checkCurrentSplit <- function(x, response, currentNode, childNodes, ldaType){
  #> x and response are already subsets
  idxRowBootstrap <- sample(nrow(x), replace = TRUE)
  idxRowOOB <- setdiff(seq_len(nrow(x)), idxRowBootstrap)
  if(length(idxRowOOB) < 5) return(1) # stop the split

  #> update the current LDA
  currentNode$nodePredict <- updateLDA(oldLDA = currentNode$nodePredict,
                                       xNew = x[idxRowBootstrap, , drop = FALSE],
                                       responseNew = response[idxRowBootstrap],
                                       missingReference = currentNode$misReference,
                                       ldaType = ldaType)
  trainIndex <- currentNode$splitFun(x = x[idxRowBootstrap,,drop = FALSE], missingReference = currentNode$misReference)
  for(i in seq_along(childNodes)) trainIndex[[i]] <- idxRowBootstrap[trainIndex[[i]]]

  #> Update the childNodes
  for(i in seq_along(childNodes)){
    #> update the LDA if the node model is LDA and there is data left
    if(length(trainIndex[[i]]) != 0){
      if(childNodes[[i]]$nodeModel == "LDA"){
        childNodes[[i]]$nodePredict <- updateLDA(oldLDA = childNodes[[i]]$nodePredict,
                                                 xNew = x[trainIndex[[i]], , drop = FALSE],
                                                 responseNew = response[trainIndex[[i]]],
                                                 missingReference = childNodes[[i]]$misReference,
                                                 ldaType = ldaType)
      }else childNodes[[i]]$nodePredict <- getMode(response[trainIndex[[i]]])
    }
  }

  # Build a treeList for easier predictions
  treeList = structure(list(), class = "SingleTreee")
  treeList[[1]] <- currentNode
  testResBefore <- predict(treeList, x[idxRowOOB,, drop = FALSE]) == as.character(response[idxRowOOB])

  for(i in seq_along(childNodes)) treeList[[i + 1]] <- childNodes[[i]]
  treeList[[1]]$children <- 1 + seq_along(childNodes)
  testResAfter <- predict(treeList, x[idxRowOOB,, drop = FALSE]) == as.character(response[idxRowOOB])
  # tTestRes <- t.test(testResAfter, testResBefore, paired = TRUE, alternative = "greater")
  tTestRes <- tryCatch({t.test(testResAfter, testResBefore, paired = TRUE, alternative = "greater")},
                       error = function(e){list(p.value = pt(sqrt(length(testResAfter)), df = length(testResAfter) - 1, lower.tail = FALSE))})
  return(tTestRes$p.value)
  # return(tTestRes$statistic)
}


# New split checking using Bootstrap ------------------------------------------------------


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

#> Here we are using bootstrap sample to evaluate the model performance
#> and decide if we are going to split more
# nBoots <- 3
# # tTestPvalue <- median(sapply(seq_len(nBoots), function(o_o) checkCurrentSplit(x = x[treeList[[currentIdx]]$idxRow,,drop = FALSE],
# #                                                                response = response[treeList[[currentIdx]]$idxRow],
# #                                                                currentNode = treeList[[currentIdx]],
# #                                                                childNodes = childNodes,
# #                                                                ldaType = ldaType)),na.rm = T)
#
# tTestPvalue <- median(sapply(seq_len(nBoots), function(o_o) generateSplitNchildren(x = x,
#                                                                                    response = response,
#                                                                                    idxCol = treeList[[currentIdx]]$idxCol,
#                                                                                    idxRow = treeList[[currentIdx]]$idxRow,
#                                                                                    treeType = treeType,
#                                                                                    ldaType = ldaType,
#                                                                                    fastTree = fastTree,
#                                                                                    missingMethod = missingMethod,
#                                                                                    splitMethod = splitMethod,
#                                                                                    minNodeSize = minNodeSize,
#                                                                                    bootstrap = TRUE)),na.rm = T)
# if(is.na(tTestPvalue) | tTestPvalue >= 0.1) next


# Forest Prediction -------------------------------------------------------

# #' @export
predict.ForestTreee <- function(object, newdata, type = "response", ...){
  if(type == "all") type = "prob" # there is no node info in Forest

  predCurrent <- lapply(object, function(treee) predict(treee, newdata = newdata, type = "prob"))
  # sometimes there are classes that not show up in the current tree
  allClassNames <- unique(unlist(lapply(predCurrent, colnames)))

  predCurrent <- lapply(predCurrent, function(matrix) {
    for (colName in allClassNames) {
      if (!(colName %in% colnames(matrix))) matrix[, colName] <- 0
    }
    return(matrix[, allClassNames])
  })
  predCurrent <- Reduce("+", predCurrent) / length(object) # get the standardized posterior

  if(type == "response"){
    predCurrent <- allClassNames[max.col(predCurrent, ties.method = "first")]
  }
  return(predCurrent)
}


# Variable Selection ------------------------------------------------------

getChiSqStat <- function(datX, y){
  sapply(datX, function(x) getChiSqStatHelper(x, y))
}

getChiSqStatHelper <- function(x,y){
  if(getNumFlag(x)){ # numerical variable: first change to factor
    m = mean(x,na.rm = T); s = sd(x,na.rm = T)
    if(sum(!is.na(x)) >= 30 * nlevels(y)){
      splitNow = c(m - s *sqrt(3)/2, m, m + s *sqrt(3)/2)
    }else splitNow = c(m - s *sqrt(3)/3, m + s *sqrt(3)/3)

    if(length(unique(splitNow)) == 1) return(0) # No possible split
    x = cut(x, breaks = c(-Inf, splitNow, Inf), right = TRUE)
  }

  if(anyNA(x)){
    levels(x) = c(levels(x), 'newLevel')
    x[is.na(x)] <- 'newLevel'
  }
  if(length(unique(x)) == 1) return(0) # No possible split

  fit <- suppressWarnings(chisq.test(x, y))

  #> Change to 1-df wilson_hilferty chi-squared stat unless
  #> the original df = 1 and p-value is larger than 10^(-16)
  ans = unname(ifelse(fit$parameter > 1L, ifelse(fit$p.value > 10^(-16),
                                                 qchisq(1-fit$p.value, df = 1),
                                                 wilson_hilferty(fit$statistic,fit$parameter)), fit$statistic))
  return(ans)
}


wilson_hilferty = function(chi, df){ # change df = K to df = 1
  ans = max(0, (7/9 + sqrt(df) * ( (chi / df) ^ (1/3) - 1 + 2 / (9 * df) ))^3)
  return(ans)
}


# svd helper ---------------------------------------------------------------


saferSVD <- function(x, ...){
  #> Target for error code 1 from Lapack routine 'dgesdd' non-convergence error
  #> Current solution: Round the design matrix to make approximations,
  #> hopefully this will solve the problem
  #>
  #> The code is a little lengthy, since the variable assignment in tryCatch is tricky
  parList <- list(svdObject = NULL,
                  svdSuccess = FALSE,
                  errorDigits = 16,
                  x = x)
  while (!parList$svdSuccess) {
    parList <- tryCatch({
      parList$svdObject <- svd(parList$x, ...)
      parList$svdSuccess <- TRUE
      parList
    }, error = function(e) {
      if (grepl("error code 1 from Lapack routine 'dgesdd'", e$message)) {
        parList$x <- round(x, digits = parList$errorDigits)
        parList$errorDigits <- parList$errorDigits - 1
        return(parList)
      } else stop(e)
    })
  }
  return(parList$svdObject)
}
