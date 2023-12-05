
# Check for input prior ---------------------------------------------------

checkPrior <- function(prior, response){
  ## Modified from randomForest.default Line 114
  if (is.null(prior)) {
    prior <- table(response) / length(response) # Default: Estimated Prior
  } else {
    if (length(prior) != nlevels(response))
      stop("length of prior not equal to number of classes")
    if (!is.null(names(prior))){
      prior <- prior[findTargetIndex(names(prior), levels(response))]
    }
    if (any(prior < 0)) stop("prior must be non-negative")
  }
  return(prior / sum(prior))
}


# Check for input misclassification cost ---------------------------------------------------

checkMisClassCost <- function(misClassCost, response){
  if (is.null(misClassCost)) {
    misClassCost <- (1 - diag(nlevels(response))) # Default: 1-identity
    colnames(misClassCost) <- rownames(misClassCost) <- levels(response)
  } else {
    if (dim(misClassCost)[1] != dim(misClassCost)[2] | dim(misClassCost)[1] != nlevels(response))
      stop("misclassification costs matrix has wrong dimension")
    if(!all.equal(colnames(misClassCost), rownames(misClassCost)))
      stop("misClassCost: colnames should be the same as rownames")
    if (!is.null(colnames(misClassCost))){
      misClassCost <- misClassCost[findTargetIndex(colnames(misClassCost), levels(response)),
                                   findTargetIndex(colnames(misClassCost), levels(response))]
    }
  }
  return(misClassCost)
}


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


# pruning -----------------------------------------------------------------

prune <- function(oldTreee,
                  x,
                  response,
                  idxCol,
                  idxRow,
                  splitMethod,
                  maxTreeLevel,
                  minNodeSize,
                  numberOfPruning,
                  missingMethod,
                  verbose){


  # Parameter Clean Up ------------------------------------------------------

  oldTreee <- updateAlphaInTree(oldTreee)
  treeeSaved = oldTreee

  # pruning and error estimate ----------------------------------------------

  idxCV <- sample(seq_len(numberOfPruning), length(response), replace = TRUE)
  savedGrove <- sapply(seq_len(numberOfPruning), function(i) new_SingleTreee(x = x,
                                                                             response = response,
                                                                             idxCol = idxCol,
                                                                             idxRow = idxRow[idxCV!=i],
                                                                             splitMethod = splitMethod,
                                                                             maxTreeLevel = maxTreeLevel,
                                                                             minNodeSize = minNodeSize,
                                                                             missingMethod = missingMethod,
                                                                             verbose = verbose), simplify = FALSE)

  treeeForPruning <- sapply(savedGrove, updateAlphaInTree, simplify = FALSE)

  CV_Table <- data.frame()
  numOfPruning <- 0
  while(TRUE){
    nodesCount <- sum(sapply(oldTreee, function(treeeNode) is.null(treeeNode$pruned)))
    if(verbose) cat("There are ",nodesCount," node(s) left in the tree.\n")
    # browser()
    meanAndSE <- getMeanAndSE(treeeListList = treeeForPruning,
                              idxCV = idxCV,
                              x = x,
                              response = response)

    currentCutAlpha = getCutAlpha(treeeList = oldTreee)
    # summary output
    CV_Table <- rbind(CV_Table, c(numOfPruning, nodesCount, meanAndSE, currentCutAlpha))
    if (nodesCount == 1) {break}

    # Cut the treee
    oldTreee <- pruneTreee(treeeList = oldTreee, alpha = currentCutAlpha)
    treeeForPruning <- sapply(treeeForPruning, function(treeeList) pruneTreee(treeeList = treeeList, alpha = currentCutAlpha), simplify = FALSE)
    numOfPruning <- numOfPruning + 1
  }

  colnames(CV_Table) <- c("treeeNo", "nodeCount", "meanMSE", "seMSE", "alpha")

  # Evaluation --------------------------------------------------------------

  kSE = 0.1
  pruneThreshold <- (CV_Table$meanMSE + kSE * CV_Table$seMSE)[which.min(CV_Table$meanMSE)]
  idxFinal <- dim(CV_Table)[1] + 1 - which.max(rev(CV_Table$meanMSE <= pruneThreshold))
  for(i in seq_len(idxFinal-1)){
    treeeSaved <- pruneTreee(treeeList = treeeSaved, alpha = CV_Table$alpha[i])
  }
  treeeNew <- dropNodes(treeeSaved)

  return(list(treeeNew = treeeNew,
              CV_Table = CV_Table,
              savedGrove = savedGrove))
}

updateAlphaInTree <- function(treeeList){
  #> Purpose: Calculate alpha, and make it monotonic
  #> assumption: the tree is a pre-order Depth-First tree.
  #> so loop backward will update the alpha correctly

  for(i in rev(seq_along(treeeList))){
    if(is.null(treeeList[[i]]$pruned)){ # Only loop over the unpruned nodes
      ## Get all terminal nodes
      treeeList[[i]]$offsprings <- getTerminalNodes(currentIdx = i, treeeList = treeeList)
      ## Get re-substitution error
      treeeList[[i]]$offspringLoss <- sum(sapply(treeeList[[i]]$offsprings, function(idx) treeeList[[idx]]$currentLoss))
      ## Get alpha: terminal nodes have NaN as their alpha
      treeeList[[i]]$alpha <- (treeeList[[i]]$currentLoss - treeeList[[i]]$offspringLoss) / (length(treeeList[[i]]$offsprings) - 1)
      ## Update alpha to be monotonic
      childrenAlpha <- sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$alpha)
      #> In case that the tree is not root, but all alpha are the same
      #> we add one to the root node alpha to separate them
      # treeeList[[i]]$alpha <- max(c(-Inf, treeeList[[i]]$alpha-1, unlist(childrenAlpha)), na.rm = TRUE) + 1
      treeeList[[i]]$alpha <- max(c(-Inf, treeeList[[i]]$alpha, unlist(childrenAlpha)), na.rm = TRUE)
    }
  }

  return(treeeList)
}

getMeanAndSE <- function(treeeListList, idxCV, x, response){

  error <- sapply(seq_along(treeeListList), function(i) sum(predict(object = treeeListList[[i]], newdata = x[idxCV==i,,drop = FALSE]) != response[idxCV==i]))

  return(c(mean(error), sd(error) / sqrt(length(treeeListList))))
}


getCutAlpha <- function(treeeList){
  # get alpha for all non-terminal nodes, NA for terminal nodes
  alphaList <- unique(sapply(treeeList, function(treeeNode) ifelse(is.null(treeeNode$children), NA, treeeNode$alpha)))
  # Geometry average
  # return(exp(mean(log(sort(alphaList)[1:2]), na.rm = TRUE)))

  # arithmetic average: since there might be negative alphas
  # return(mean(sort(alphaList)[1:2], na.rm = TRUE))
  return(min(alphaList,na.rm = TRUE))
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

pruneTreee <- function(treeeList, trainErrorCap, alpha){
  for(i in rev(seq_along(treeeList))){
    treeeNode <- treeeList[[i]]
    # not yet pruned + non-terminal node + alpha below threshold
    currentFlag <- is.null(treeeNode$pruned) & !is.null(treeeNode$children) & treeeNode$alpha <= alpha + 1e-10 # R rounding error, 1e-10 needed
    if(currentFlag){
      allChildren <- getTerminalNodes(currentIdx = i, treeeList = treeeList, keepNonTerminal = TRUE)
      for(j in setdiff(allChildren, i)) {treeeList[[j]]$pruned <- TRUE}
      treeeList[[i]]["children"] <- list(NULL) # cut the branch
    }
  }
  treeeList <- makeAlphaMono(treeeList = treeeList, trainErrorCap = trainErrorCap)
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



# Splitting ---------------------------------------------------------------

### Sum of Pillai's trace ###

getSplitFunPillai <- function(x = x, response = response, modelLDA = modelLDA){
  x <- getDesignMatrix(modelLDA = modelLDA, data = x) # get the scaled x
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

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
    x <- getDesignMatrix(modelLDA = modelLDA, data = fixedData)
    projectionOnSplit <- unname(as.vector(x %*% matrix(aScaled[-1], ncol = 1)))
    # return the relative index
    return(list(which(projectionOnSplit < aScaled[1]), which(projectionOnSplit >= aScaled[1])))
  }
  attr(res, "a") <- aFinal # record the split function's form
  return(res)
}


### Gini Split ###

getSplitFunLDscores <- function(x, response, modelLDA){

  LDscore <- getLDscores(modelLDA = modelLDA, data = x, nScores = 1)
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

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
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



# fixNewLevel -------------------------------------------------------------


fixNewLevel <- function(datTest, datTrain){
  #> change the shape of test data to the training data
  #> and make sure that the dimension of the data is the same as missingRefernce

  nameVarIdx <- match(colnames(datTrain), colnames(datTest))
  if(anyNA(nameVarIdx)){
    #> New columns fix (or Flags): If there are less columns than it should be,
    #> add columns with NA
    datTest[,colnames(datTrain)[which(is.na(nameVarIdx))]] <- NA
    nameVarIdx <- match(colnames(datTrain), colnames(datTest))
  }
  datTest <- datTest[,nameVarIdx, drop = FALSE]
  #> The columns are the same, now fix new levels

  #> change all characters to factors
  idxC <- which(sapply(datTrain, class) == "character")
  if(length(idxC) != 0){
    for(i in idxC){
      datTrain[,i] <- as.factor(datTrain[,i])
    }
  }

  #> New levels fix
  levelReference <- sapply(datTrain, levels, simplify = FALSE)
  for(i in which(!sapply(levelReference,is.null))){
    # sapply would lose factor property but left character, why?
    datTest[,i] <- factor(datTest[,i], levels = levelReference[[i]])
  }

  return(datTest)
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
