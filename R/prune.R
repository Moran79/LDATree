#' Prune a Decision Tree
#'
#' This function performs pruning on an existing decision tree by using
#' cross-validation to assess the best level of tree complexity, balancing the
#' trade-off between tree size and predictive performance.
#'
#' @noRd
prune <- function(oldTreee,
                  datX,
                  response,
                  ldaType,
                  nodeModel,
                  numberOfPruning,
                  maxTreeLevel,
                  minNodeSize,
                  pThreshold,
                  prior,
                  misClassCost,
                  missingMethod,
                  kSample,
                  verbose){

  treeeSaved = oldTreee # saved a copy for cutting
  idxCV <- sample(seq_len(numberOfPruning), length(response), replace = TRUE)
  treeeForPruning <- vector(mode = "list", length = numberOfPruning)

  if(verbose){
    cat('\nPruning: Grow the trees...\n')
    pbGrowTree <- utils::txtProgressBar(min = 0, max = numberOfPruning, style = 3)
  }

  for(i in seq_len(numberOfPruning)){
    treeeForPruning[[i]] <- predict(new_SingleTreee(datX = datX[idxCV != i, , drop = FALSE],
                                                    response = response[idxCV != i],
                                                    ldaType = ldaType,
                                                    nodeModel = nodeModel,
                                                    maxTreeLevel = maxTreeLevel,
                                                    minNodeSize = minNodeSize,
                                                    pThreshold = pThreshold,
                                                    prior = prior,
                                                    misClassCost = misClassCost,
                                                    missingMethod = missingMethod,
                                                    kSample = kSample,
                                                    verbose = FALSE),
                                    newdata = datX[idxCV == i, , drop = FALSE],
                                    insideCV = TRUE,
                                    obsY = response[idxCV == i])
    if(verbose) utils::setTxtProgressBar(pbGrowTree, i)
  }

  if(verbose){
    cat('\nPruning: Prune the trees...\n')
    pbPruneTree <- utils::txtProgressBar(min = 0, max = length(treeeSaved), style = 3)
  }

  CV_Table <- data.frame(); pruningCounter <- 0

  while(TRUE){
    nodesCount <- sum(sapply(oldTreee, function(treeeNode) is.null(treeeNode$pruned)))
    if(verbose) utils::setTxtProgressBar(pbPruneTree, length(treeeSaved) - nodesCount)

    meanAndSE <- getMeanAndSE(treeeListList = treeeForPruning)
    currentCutAlpha = getCutAlpha(treeeList = oldTreee)
    CV_Table <- rbind(CV_Table, c(pruningCounter, nodesCount, meanAndSE, currentCutAlpha))
    if (nodesCount == 1) {
      if(verbose) utils::setTxtProgressBar(pbPruneTree, length(treeeSaved))
      break
    }

    oldTreee <- pruneTreee(treeeList = oldTreee, alpha = currentCutAlpha)
    treeeForPruning <- lapply(treeeForPruning, function(treeeList) pruneTreee(treeeList = treeeList, alpha = currentCutAlpha))
    pruningCounter <- pruningCounter + 1
  }

  colnames(CV_Table) <- c("treeeNo", "nodeCount", "meanMSE", "seMSE", "alpha")

  # Evaluation --------------------------------------------------------------

  kSE = 0.25
  pruneThreshold <- (CV_Table$meanMSE + kSE * CV_Table$seMSE)[which.min(CV_Table$meanMSE)]
  idxFinal <- dim(CV_Table)[1] + 1 - which.max(rev(CV_Table$meanMSE <= pruneThreshold))
  for(i in seq_len(idxFinal-1)) treeeSaved <- pruneTreee(treeeList = treeeSaved, alpha = CV_Table$alpha[i])

  return(list(treeeNew = dropNodes(treeeSaved), CV_Table = CV_Table))
}


#' Update Alpha Values in a Tree
#'
#' This function calculates and updates the alpha values for each node in the
#' decision tree. Alpha is a measure (p-value) of the improvement in loss after
#' pruning a subtree.
#'
#' @noRd
updateAlphaInTree <- function(treeeList){
  #> assumption: the index of the children must be larger than its parent's
  for(i in rev(seq_along(treeeList))){
    if(is.null(treeeList[[i]]$pruned)){ # Only loop over the unpruned nodes

      ## Get all terminal nodes
      treeeList[[i]]$childrenTerminal <- getTerminalNodes(currentIdx = i, treeeList = treeeList)

      ## Get re-substitution error
      treeeList[[i]]$childrenTerminalLoss <- sum(sapply(treeeList[[i]]$childrenTerminal, function(idx) treeeList[[idx]]$currentLoss))
      treeeList[[i]]$alpha <- getOneSidedPvalue(N = length(treeeList[[i]]$idxRow),
                                                lossBefore = treeeList[[i]]$currentLoss,
                                                lossAfter = treeeList[[i]]$childrenTerminalLoss)
    }
  }
  return(treeeList)
}


#' Calculate Mean and Standard Error of Tree Errors
#'
#' This function computes the mean and standard error of the errors from a list
#' of decision trees.
#'
#' @noRd
#'
#' @param treeeListList A list of decision tree objects, where each tree
#'   contains information about its error.
#'
getMeanAndSE <- function(treeeListList){
  error <- sapply(treeeListList, getMeanAndSEhelper)
  return(c(mean(error), stats::sd(error) / sqrt(length(treeeListList))))
}

getMeanAndSEhelper <- function(treeeList){
  terminalIdx <- getTerminalNodes(currentIdx = 1, treeeList = treeeList)
  return(sum(do.call(c, lapply(terminalIdx, function(i) treeeList[[i]]$CVerror))))
}


#' Calculate the Cut-off Alpha for Pruning
#'
#' This function calculates the cut-off alpha value for pruning the decision
#' tree. The alpha value is computed from the internal (non-terminal) nodes of
#' the tree, representing the improvement in loss after pruning. The geometric
#' mean of the two largest alpha values is returned as the cut-off.
#'
#' @noRd
getCutAlpha <- function(treeeList){
  internalIdx <- setdiff(getTerminalNodes(1, treeeList, keepNonTerminal = T),
                         getTerminalNodes(1, treeeList, keepNonTerminal = F))
  alphaList <- unique(do.call(c, lapply(internalIdx, function(i) treeeList[[i]]$alpha)))
  alphaCandidates <- stats::na.omit(sort(alphaList, decreasing = TRUE)[1:2])

  if(length(alphaCandidates) == 0) alphaCandidates <- 1
  return(exp(mean(log(alphaCandidates))))
}


#' Prune a Decision Tree Based on Alpha
#'
#' This function prunes the decision tree by removing branches whose alpha
#' values are greater than or equal to the given threshold. Pruning is performed
#' on non-terminal nodes that meet the alpha condition, and their child nodes
#' are marked as pruned. After pruning, the alpha values of the remaining nodes
#' are updated.
#'
#' @noRd
pruneTreee <- function(treeeList, alpha){
  for(i in rev(seq_along(treeeList))){
    treeeNode <- treeeList[[i]]
    # not yet pruned + non-terminal node + alpha below threshold
    currentFlag <- is.null(treeeNode$pruned) && !is.null(treeeNode$children) && treeeNode$alpha >= alpha - 1e-10 # R rounding error, 1e-10 needed

    if(currentFlag){
      allChildren <- getTerminalNodes(currentIdx = i, treeeList = treeeList, keepNonTerminal = TRUE)
      for(j in setdiff(allChildren, i)) {treeeList[[j]]$pruned <- TRUE}
      treeeList[[i]]["children"] <- list(NULL) # cut the branch
    }
  }
  treeeList <- updateAlphaInTree(treeeList = treeeList)
  return(treeeList)
}


#' Drop Pruned Nodes from a Decision Tree
#'
#' This function removes pruned nodes from the decision tree, keeping only
#' terminal and relevant internal nodes. It reassigns the indices of the
#' remaining nodes and updates their children and parent node references
#' accordingly.
#'
#' @noRd
dropNodes <- function(treeeList){
  finalNodeIdx <- sort(getTerminalNodes(treeeList = treeeList, currentIdx = 1, keepNonTerminal = TRUE))
  treeeList <- treeeList[finalNodeIdx]
  class(treeeList) <- "Treee"
  for(i in seq_along(treeeList)){
    treeeList[[i]]$currentIndex <- i # re-assign the currentIndex
    if(!is.null(treeeList[[i]]$children)) {
      treeeList[[i]]$children <- sapply(treeeList[[i]]$children, function(x) which(finalNodeIdx == x))
    }else{
      if(treeeList[[i]]$stopInfo == "Normal") treeeList[[i]]$stopInfo = "Pruned"
    }
    treeeList[[i]]$parent <- which(finalNodeIdx == treeeList[[i]]$parent)
  }
  return(treeeList)
}


#' Retrieve Terminal Nodes in a Decision Tree
#'
#' This function retrieves all terminal (leaf) nodes that are descendants of a
#' specified node in the decision tree. Optionally, intermediate nodes can also
#' be included in the result.
#'
#' @noRd
getTerminalNodes <- function(currentIdx, treeeList, keepNonTerminal = FALSE){
  #> Get all terminal nodes that are offsprings from currentIdx
  treeeNode <- treeeList[[currentIdx]]
  if(is.null(treeeNode$children)){
    return(currentIdx)
  }else{
    nonTerminalNodes <- if(keepNonTerminal) currentIdx
    terminalNodes <- do.call(c, lapply(treeeNode$children,function(x) getTerminalNodes(currentIdx = x, treeeList = treeeList, keepNonTerminal = keepNonTerminal)))
    return(c(nonTerminalNodes, terminalNodes))
  }
}


pruneTreeeByLevel <- function(treeeList, K){
  #> Internal exploration only
  for(i in rev(seq_along(treeeList))){
    treeeNode <- treeeList[[i]]
    # not yet pruned + non-terminal node + alpha below threshold
    currentFlag <- is.null(treeeNode$pruned) && !is.null(treeeNode$children) && treeeNode$currentLevel >= K # R rounding error, 1e-10 needed

    if(currentFlag){
      allChildren <- getTerminalNodes(currentIdx = i, treeeList = treeeList, keepNonTerminal = TRUE)
      for(j in setdiff(allChildren, i)) {treeeList[[j]]$pruned <- TRUE}
      treeeList[[i]]["children"] <- list(NULL) # cut the branch
    }
  }
  treeeList <- dropNodes(treeeList = treeeList)
  return(treeeList)
}
