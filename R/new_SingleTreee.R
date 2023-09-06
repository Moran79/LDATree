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

  #> Some notes to clarify the parameters...

  treeList = structure(list(), class = "SingleTreee") # save the tree

  # first Node
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

      if(is.null(treeList[[currentIdx]]$splitGini)) next # No cut due to ties
      treeList[[currentIdx]]$splitCut <- treeList[[currentIdx]]$splitGini$cut

      # distribute the validation set
      fixedData <- getDataInShape(data = xValidation[treeList[[currentIdx]]$idxRowValidation,,drop = FALSE], missingReference = treeList[[currentIdx]]$misReference)
      currentScore <- getLDscores(modelLDA = treeList[[currentIdx]]$nodePredict, data = fixedData, nScores = 1)
      leftIdx <- (currentScore <= treeList[[currentIdx]]$splitCut)

      # browser()
      # get child nodes
      leftNode <- new_TreeeNode(x = x,
                                response = response,
                                xValidation = xValidation,
                                responseValidation = responseValidation,
                                idxCol = treeList[[currentIdx]]$idxCol,
                                idxRow = treeList[[currentIdx]]$splitGini$left,
                                idxRowValidation = treeList[[currentIdx]]$idxRowValidation[leftIdx],
                                missingMethod = missingMethod,
                                splitMethod = splitMethod,
                                maxTreeLevel = maxTreeLevel,
                                minNodeSize = minNodeSize,
                                currentLevel = treeList[[currentIdx]]$currentLevel + 1,
                                parentIndex = currentIdx)
      rightNode <- new_TreeeNode(x = x,
                                 response = response,
                                 xValidation = xValidation,
                                 responseValidation = responseValidation,
                                 idxCol = treeList[[currentIdx]]$idxCol,
                                 idxRow = treeList[[currentIdx]]$splitGini$right,
                                 idxRowValidation = treeList[[currentIdx]]$idxRowValidation[!leftIdx],
                                 missingMethod = missingMethod,
                                 splitMethod = splitMethod,
                                 maxTreeLevel = maxTreeLevel,
                                 minNodeSize = minNodeSize,
                                 currentLevel = treeList[[currentIdx]]$currentLevel + 1,
                                 parentIndex = currentIdx)

      # update alpha & lag
      treeList[[currentIdx]]$alpha <- treeList[[currentIdx]]$currentLoss - leftNode$currentLoss - rightNode$currentLoss
      if(treeList[[currentIdx]]$alpha >= 0){ # if non-negative alpha, refresh the counter
        treeList[[currentIdx]]$lag = 0
      }else{ # if negative alpha, lag += 1
        treeList[[currentIdx]]$lag = treeList[[currentIdx]]$lag + 1
      }

      leftNode$lag <- rightNode$lag <- treeList[[currentIdx]]$lag
      if(treeList[[currentIdx]]$lag > kStepAhead){
        treeList[[currentIdx]]$stopFlag = 4
        next
      }

      # Put child nodes in the tree
      leftIdx <- length(treeList) + 1
      rightIdx <- leftIdx + 1
      treeList[[currentIdx]]$children <- c(leftIdx, rightIdx)
      nodeStack <- c(nodeStack, leftIdx, rightIdx)
      treeList[[leftIdx]] <- leftNode
      treeList[[rightIdx]] <- rightNode
    }
  }
  return(treeList)
}
