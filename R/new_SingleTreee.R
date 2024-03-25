new_SingleTreee <- function(datX,
                            response,
                            ldaType,
                            nodeModel,
                            missingMethod,
                            prior,
                            maxTreeLevel,
                            minNodeSize,
                            pThreshold,
                            verbose){

  treeeList = structure(list(), class = "SingleTreee") # save the tree

  ### Initialize the first Node ###

  nodeStack <- c(1)
  treeeList[[1]] <- new_TreeeNode(datX = datX,
                                 response = response,
                                 idxCol = seq_len(ncol(datX)),
                                 idxRow = seq_len(nrow(datX)),
                                 ldaType = ldaType,
                                 nodeModel = nodeModel,
                                 missingMethod = missingMethod,
                                 prior = prior,
                                 maxTreeLevel = maxTreeLevel,
                                 minNodeSize = minNodeSize,
                                 currentLevel = 0,
                                 parentIndex = 0)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[1]; nodeStack <- nodeStack[-1] # pop the first element
    if(verbose) cat("The current index is:", currentIdx, "\n")

    if(treeeList[[currentIdx]]$stopFlag == 0){ # if it has (potential) child nodes

      trainIndex <- attr(treeeList[[currentIdx]]$splitFun, "splitRes") # distribute the training set

      ### Get child nodes ###

      childNodes <- lapply(seq_along(trainIndex),
                           function(i) new_TreeeNode(datX = datX,
                                                     response = response,
                                                     idxCol = treeeList[[currentIdx]]$idxCol,
                                                     idxRow = treeeList[[currentIdx]]$idxRow[trainIndex[[i]]],
                                                     ldaType = ldaType,
                                                     nodeModel = nodeModel,
                                                     missingMethod = missingMethod,
                                                     prior = prior,
                                                     maxTreeLevel = maxTreeLevel,
                                                     minNodeSize = minNodeSize,
                                                     currentLevel = treeeList[[currentIdx]]$currentLevel + 1,
                                                     parentIndex = currentIdx))


      ### Stopping & pruning ###
      #> 1. update the p-value for loss drop
      lossBefore <- treeeList[[currentIdx]]$currentLoss
      lossAfter <- do.call(sum,lapply(childNodes, function(node) node$currentLoss))
      treeeList[[currentIdx]]$alpha <-  getOneSidedPvalue(N = length(treeeList[[currentIdx]]$idxRow),
                                                         lossBefore = lossBefore,
                                                         lossAfter = lossAfter)

      #> 2. pre-stopping
      if(treeeList[[currentIdx]]$alpha >= pThreshold){
        treeeList[[currentIdx]]$stopFlag = 6
        next
      }

      ### Put child nodes in the tree ###

      childIdx <- seq_along(childNodes) + length(treeeList)
      treeeList[[currentIdx]]$children <- childIdx
      nodeStack <- c(nodeStack, childIdx)
      for(i in seq_along(childIdx)) treeeList[[childIdx[i]]] <- childNodes[[i]]
    }
  }

  for(i in seq_along(treeeList)) treeeList[[i]]$currentIndex <- i # assign the currentIndex

  return(treeeList)
}
