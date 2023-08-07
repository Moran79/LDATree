prune <- function(oldTreee,
                  x,
                  response,
                  idxCol,
                  idxRow,
                  splitMethod,
                  prior,
                  weights,
                  maxTreeLevel,
                  minNodeSize,
                  numberOfPrune,
                  misClassCost,
                  missingMethod,
                  randomSeed){


  # Parameter Clean Up ------------------------------------------------------

  if (!is.null(randomSeed)) {set.seed(randomSeed)} # user input seed
  oldTreee <- updateAlphaInTree(oldTreee)
  treeeSaved = oldTreee

  # pruning and error estimate ----------------------------------------------

  idxCV <- sample(seq_len(numberOfPrune), length(response), replace = TRUE)
  savedGrove <- sapply(seq_len(numberOfPrune), function(i) new_SingleTreee(x = x,
                                                                           response = response,
                                                                           idxCol = idxCol,
                                                                           idxRow = idxRow[idxCV!=i],
                                                                           splitMethod = splitMethod,
                                                                           prior = prior,
                                                                           weights = weights,
                                                                           maxTreeLevel = maxTreeLevel,
                                                                           minNodeSize = minNodeSize,
                                                                           misClassCost = misClassCost,
                                                                           missingMethod = missingMethod), simplify = FALSE)

  treeeForPruning <- sapply(savedGrove, updateAlphaInTree, simplify = FALSE)

  CV_Table <- data.frame()
  numOfPruning <- 0
  while(TRUE){
    nodesCount <- sum(sapply(oldTreee, function(treeeNode) is.null(treeeNode$pruned)))
    message("There are ",nodesCount," node(s) left in the tree.")
    meanAndSE <- getMeanAndSE(treeeListList = treeeForPruning,
                              idxCV = idxCV,
                              x = x,
                              response = response,
                              misClassCost = misClassCost)

    currentCutAlpha = getCutAlpha(treeeList = oldTreee)
    # summary output
    CV_Table <- rbind(CV_Table, c(numOfPruning, nodesCount, meanAndSE, currentCutAlpha))
    if (nodesCount == 1) {break}

    # Cut the treee
    oldTreee <- pruneTreee(treeeList = oldTreee, alpha = currentCutAlpha)
    for(i in seq_len(numberOfPrune)){
      treeeForPruning[[i]] <- pruneTreee(treeeList = treeeForPruning[[i]], alpha = currentCutAlpha)
    }

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
  #> assumption: the tree is complete in order
  #> so loop backward will update the alpha correctly

  for(i in rev(seq_along(treeeList))){
    if(is.null(treeeList[[i]]$pruned)){ # Only loop over the unpruned nodes
      ## Get all terminal nodes
      treeeList[[i]]$offsprings <- getTerminalNodes(currentIdx = i, treeeList = treeeList)
      ## Get re-substitution error
      treeeList[[i]]$childrenLoss <- sum(sapply(treeeList[[i]]$offsprings, function(idx) treeeList[[idx]]$currentLoss))
      ## Get alpha: terminal nodes have NaN as their alpha
      treeeList[[i]]$alpha <- (treeeList[[i]]$currentLoss - treeeList[[i]]$childrenLoss) / (length(treeeList[[i]]$offsprings) - 1)
      ## Update alpha to be monotonic
      childrenAlpha <- sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$alpha)
      treeeList[[i]]$alpha <- max(c(-1, treeeList[[i]]$alpha-1, unlist(childrenAlpha)), na.rm = TRUE) + 1
    }
  }

  #> In case that the tree is not root, but all alpha are the same
  #> we add one to the root node alpha to separate them

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

getAlpha <- function(treeeNode){
  currentFlag <- !is.null(treeeNode$children)
  return(ifelse(currentFlag, treeeNode$alpha, NA))
}

getCutAlpha <- function(treeeList){
  alphaList <- unique(sapply(treeeList, getAlpha))
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

getMeanAndSE <- function(treeeListList, idxCV, x, response, misClassCost){

  error <- sapply(seq_along(treeeListList), function(i) sum(table(predict(object = treeeListList[[i]], newdata = x[idxCV==i,])$pred, response[idxCV==i]) * misClassCost))

  return(c(mean(error), sd(error) / sqrt(length(treeeListList))))
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


