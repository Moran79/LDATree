
# pruning -----------------------------------------------------------------

prune <- function(oldTreee,
                  numberOfPruning,
                  datX,
                  response,
                  ldaType,
                  nodeModel,
                  missingMethod,
                  prior,
                  maxTreeLevel,
                  minNodeSize,
                  pThreshold,
                  verbose){


  # Parameter Clean Up ------------------------------------------------------

  oldTreee <- updateAlphaInTree(oldTreee) # make the alpha monotone
  treeeSaved = oldTreee

  # pruning and error estimate ----------------------------------------------

  idxCV <- sample(seq_len(numberOfPruning), length(response), replace = TRUE)
  treeeForPruning <- vector(mode = "list", length = numberOfPruning)

  if(verbose){
    cat('\nPruning: Grow the trees...\n')
    pbGrowTree <- utils::txtProgressBar(min = 0, max = numberOfPruning, style = 3)
  }

  for(i in seq_len(numberOfPruning)){
    treeeForPruning[[i]] <- predict(new_SingleTreee(datX = datX[idxCV!=i, , drop = FALSE],
                                                    response = response[idxCV!=i],
                                                    ldaType = ldaType,
                                                    nodeModel = nodeModel,
                                                    missingMethod = missingMethod,
                                                    prior = prior,
                                                    maxTreeLevel = maxTreeLevel,
                                                    minNodeSize = minNodeSize,
                                                    pThreshold = pThreshold,
                                                    verbose = FALSE),
                                    newdata = datX[idxCV ==i, , drop = FALSE],
                                    insideCV = TRUE,
                                    newY = response[idxCV==i])
    if(verbose) utils::setTxtProgressBar(pbGrowTree, i)
  }

  treeeForPruning <- lapply(treeeForPruning, updateAlphaInTree)

  CV_Table <- data.frame()
  numOfPruning <- 0
  if(verbose){
    cat('\nPruning: Prune the trees...\n')
    pbPruneTree <- utils::txtProgressBar(min = 0, max = length(treeeSaved), style = 3)
  }

  while(TRUE){
    nodesCount <- sum(sapply(oldTreee, function(treeeNode) is.null(treeeNode$pruned)))
    if(verbose) utils::setTxtProgressBar(pbPruneTree, length(treeeSaved) - nodesCount)

    meanAndSE <- getMeanAndSE(treeeListList = treeeForPruning)

    currentCutAlpha = getCutAlpha(treeeList = oldTreee)

    # summary output
    CV_Table <- rbind(CV_Table, c(numOfPruning, nodesCount, meanAndSE, currentCutAlpha))
    if (nodesCount == 1) {
      if(verbose) utils::setTxtProgressBar(pbPruneTree, length(treeeSaved))
      break
    }

    # Cut the treee
    oldTreee <- pruneTreee(treeeList = oldTreee, alpha = currentCutAlpha)
    treeeForPruning <- lapply(treeeForPruning, function(treeeList) pruneTreee(treeeList = treeeList, alpha = currentCutAlpha))
    numOfPruning <- numOfPruning + 1
  }

  colnames(CV_Table) <- c("treeeNo", "nodeCount", "meanMSE", "seMSE", "alpha")

  # Evaluation --------------------------------------------------------------

  kSE = 0.5
  pruneThreshold <- (CV_Table$meanMSE + kSE * CV_Table$seMSE)[which.min(CV_Table$meanMSE)]
  idxFinal <- dim(CV_Table)[1] + 1 - which.max(rev(CV_Table$meanMSE <= pruneThreshold))
  for(i in seq_len(idxFinal-1)){
    treeeSaved <- pruneTreee(treeeList = treeeSaved, alpha = CV_Table$alpha[i])
  }
  treeeNew <- dropNodes(treeeSaved)

  return(list(treeeNew = treeeNew,
              CV_Table = CV_Table))
}

updateAlphaInTree <- function(treeeList){
  #> Purpose: Calculate alpha, and make it monotonic
  #> assumption: the index of the children must be larger than its parent's

  for(i in rev(seq_along(treeeList))){
    if(is.null(treeeList[[i]]$pruned)){ # Only loop over the unpruned nodes

      ## Get all terminal nodes
      treeeList[[i]]$offsprings <- getTerminalNodes(currentIdx = i, treeeList = treeeList)
      ## Get re-substitution error
      treeeList[[i]]$offspringLoss <- sum(sapply(treeeList[[i]]$offsprings, function(idx) treeeList[[idx]]$currentLoss))
      # treeeList[[i]]$alpha <- (treeeList[[i]]$currentLoss - treeeList[[i]]$offspringLoss) / (length(treeeList[[i]]$offsprings) - 1)
      treeeList[[i]]$alpha <- getOneSidedPvalue(N = length(treeeList[[i]]$idxRow),
                                                lossBefore = treeeList[[i]]$currentLoss,
                                                lossAfter = treeeList[[i]]$offspringLoss)
      ## Update alpha to be monotonic
      childrenAlpha <- sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$alpha)
      #> In case that the tree is not root, but all alpha are the same
      #> we add one to the root node alpha to separate them
      treeeList[[i]]$alpha <- min(c(Inf, treeeList[[i]]$alpha, unlist(childrenAlpha)), na.rm = TRUE)
    }
  }

  return(treeeList)
}


getMeanAndSE <- function(treeeListList){
  error <- sapply(treeeListList, getMeanAndSEhelper)
  return(c(mean(error), sd(error) / sqrt(length(treeeListList))))
}

getMeanAndSEhelper <- function(treeeList){
  terminalIdx <- getTerminalNodes(currentIdx = 1, treeeList = treeeList)
  return(sum(do.call(c, lapply(terminalIdx, function(i) treeeList[[i]]$CVerror))))
}


getCutAlpha <- function(treeeList){
  # get alpha for all non-terminal nodes, NA for terminal nodes
  alphaList <- unique(sapply(treeeList, function(treeeNode) ifelse(is.null(treeeNode$children), NA, treeeNode$alpha)))

  #> Geometry average
  #> When no alpha is available, use .Machine$double.eps instead
  alphaCandidates <- na.omit(sort(alphaList, decreasing = TRUE)[1:2])
  if(length(alphaCandidates) == 0) alphaCandidates <- 0
  return(exp(mean(log(alphaCandidates))))
}


pruneTreee <- function(treeeList, alpha){
  for(i in rev(seq_along(treeeList))){
    treeeNode <- treeeList[[i]]
    # not yet pruned + non-terminal node + alpha below threshold
    currentFlag <- is.null(treeeNode$pruned) & !is.null(treeeNode$children) & treeeNode$alpha >= alpha - 1e-10 # R rounding error, 1e-10 needed

    if(currentFlag){
      allChildren <- getTerminalNodes(currentIdx = i, treeeList = treeeList, keepNonTerminal = TRUE)
      for(j in setdiff(allChildren, i)) {treeeList[[j]]$pruned <- TRUE}
      treeeList[[i]]["children"] <- list(NULL) # cut the branch
    }
  }
  treeeList <- updateAlphaInTree(treeeList = treeeList)
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
