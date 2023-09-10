new_SingleTreee <- function(x,
                            response,
                            xValidation,
                            responseValidation,
                            missingMethod,
                            splitMethod,
                            maxTreeLevel,
                            minNodeSize,
                            kStepAhead,
                            verbose){

  treeList = structure(list(), class = "SingleTreee") # save the tree

  # initialize the first Node
  nodeStack <- c(1)
  treeList[[1]] <- new_TreeeNode(x = x,
                                 response = response,
                                 xValidation = xValidation,
                                 responseValidation = responseValidation,
                                 idxCol = seq_len(ncol(x)),
                                 idxRow = seq_len(nrow(x)),
                                 idxRowValidation = seq_len(nrow(xValidation)),
                                 missingMethod = missingMethod,
                                 splitMethod = splitMethod,
                                 maxTreeLevel = maxTreeLevel,
                                 minNodeSize = minNodeSize,
                                 currentLevel = 0,
                                 parentIndex = 0)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[1]; nodeStack <- nodeStack[-1] # pop the first element

    if(treeList[[currentIdx]]$stopFlag == 0){ # if splitting goes on

      # distribute the training and validation set
      validIndex <- treeList[[currentIdx]]$splitFun(x = xValidation[treeList[[currentIdx]]$idxRowValidation,,drop = FALSE], missingReference = treeList[[currentIdx]]$misReference)
      trainIndex <- treeList[[currentIdx]]$splitFun(x = x[treeList[[currentIdx]]$idxRow,,drop = FALSE], missingReference = treeList[[currentIdx]]$misReference)

      # get child nodes
      childNodes <- lapply(seq_along(trainIndex), function(i) new_TreeeNode(x = x,
                                                                      response = response,
                                                                      xValidation = xValidation,
                                                                      responseValidation = responseValidation,
                                                                      idxCol = treeList[[currentIdx]]$idxCol,
                                                                      idxRow = treeList[[currentIdx]]$idxRow[trainIndex[[i]]],
                                                                      idxRowValidation = treeList[[currentIdx]]$idxRowValidation[validIndex[[i]]],
                                                                      missingMethod = missingMethod,
                                                                      splitMethod = splitMethod,
                                                                      maxTreeLevel = maxTreeLevel,
                                                                      minNodeSize = minNodeSize,
                                                                      currentLevel = treeList[[currentIdx]]$currentLevel + 1,
                                                                      parentIndex = currentIdx))

      # update alpha & lag

      treeList[[currentIdx]]$alpha <- (treeList[[currentIdx]]$currentLoss - do.call(sum,lapply(childNodes, function(node) node$currentLoss))) / (length(childNodes) - 1)
      if(treeList[[currentIdx]]$alpha >= 0){ # if non-negative alpha, refresh the counter
        treeList[[currentIdx]]$lag = 0
      }else{ # if negative alpha, lag += 1
        treeList[[currentIdx]]$lag = treeList[[currentIdx]]$lag + 1
      }

      for(i in seq_along(childNodes)) childNodes[[i]]$lag <- treeList[[currentIdx]]$lag
      if(treeList[[currentIdx]]$lag > kStepAhead){
        treeList[[currentIdx]]$stopFlag = 5
        next
      }

      # Put child nodes in the tree
      childIdx <- seq_along(childNodes) + length(treeList)
      treeList[[currentIdx]]$children <- childIdx
      nodeStack <- c(nodeStack, childIdx)
      for(i in seq_along(childIdx)) treeList[[childIdx[i]]] <- childNodes[[i]]
    }
  }
  return(treeList)
}
