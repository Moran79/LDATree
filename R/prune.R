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
                                                                           missingMethod = missingMethod), simplify = FALSE)

  treeeForPruning <- sapply(savedGrove, updateAlphaInTree, simplify = FALSE)

  CV_Table <- data.frame()
  numOfPruning <- 0
  while(TRUE){
    nodesCount <- sum(sapply(oldTreee, function(treeeNode) is.null(treeeNode$pruned)))
    if(verbose) cat("There are ",nodesCount," node(s) left in the tree.\n")
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
      treeeList[[i]]$alpha <- max(c(-1, treeeList[[i]]$alpha-1, unlist(childrenAlpha)), na.rm = TRUE) + 1
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

getMeanAndSE <- function(treeeListList, idxCV, x, response){

  error <- sapply(seq_along(treeeListList), function(i) sum(predict(object = treeeListList[[i]], newdata = x[idxCV==i,,drop = FALSE]) != response[idxCV==i]))

  return(c(mean(error), sd(error) / sqrt(length(treeeListList))))
}


getCutAlpha <- function(treeeList){
  # get alpha for all terminal nodes
  alphaList <- unique(sapply(treeeList, function(treeeNode) ifelse(is.null(treeeNode$children), NA, treeeNode$alpha)))
  # Geometry average
  return(exp(mean(log(sort(alphaList)[1:2]), na.rm = TRUE)))
}

pruneTreee <- function(treeeList, alpha){
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
  treeeList <- updateAlphaInTree(treeeList = treeeList)
  return(treeeList)
}


dropNodes <- function(treeeList){
  finalNodeIdx <- sort(getTerminalNodes(treeeList = treeeList, currentIdx = 1, keepNonTerminal = TRUE))
  treeeList <- treeeList[finalNodeIdx]
  class(treeeList) <- "SingleTreee"
  for(i in seq_along(treeeList)){
    treeeList[[i]]$currentIndex <- i
    if(!is.null(treeeList[[i]]$children)) {
      treeeList[[i]]$children <- sapply(treeeList[[i]]$children, function(x) which(finalNodeIdx == x))
    }
    treeeList[[i]]$parent <- which(finalNodeIdx == treeeList[[i]]$parent)
  }
  return(treeeList)
}


