
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
