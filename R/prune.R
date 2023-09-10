makeAlphaMono <- function(treeeList){
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
        currentAlpha <- (treeeList[[i]]$currentLoss - sum(sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$currentLoss))) / (numOfChildren - 1)
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
  treeeList <- makeAlphaMono(treeeList = treeeList)
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
    }
    treeeList[[i]]$parent <- which(finalNodeIdx == treeeList[[i]]$parent)
  }
  return(treeeList)
}


