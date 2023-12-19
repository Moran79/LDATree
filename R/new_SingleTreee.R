new_SingleTreee <- function(datX,
                            response,
                            treeType,
                            splitMethod,
                            ldaType,
                            fastTree,
                            nodeModel,
                            missingMethod,
                            maxTreeLevel,
                            minNodeSize,
                            trainErrorCap,
                            verbose){

  treeList = structure(list(), class = "SingleTreee") # save the tree

  #> Forest: Bootstrap samples for each single tree
  if(treeType == "forest") idxRowRoot <- sample(length(response), length(response), replace = TRUE)
  else idxRowRoot <- seq_len(nrow(datX))

  #> initialize the first Node
  nodeStack <- c(1)
  treeList[[1]] <- new_TreeeNode(datX = datX,
                                 response = response,
                                 idxCol = seq_len(ncol(datX)),
                                 idxRow = idxRowRoot,
                                 treeType = treeType,
                                 splitMethod = splitMethod,
                                 ldaType = ldaType,
                                 fastTree = fastTree,
                                 nodeModel = nodeModel,
                                 missingMethod = missingMethod,
                                 maxTreeLevel = maxTreeLevel,
                                 minNodeSize = minNodeSize,
                                 currentLevel = 0,
                                 parentIndex = 0)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[1]; nodeStack <- nodeStack[-1] # pop the first element
    cat("The current index is:", currentIdx, "\n")

    if(treeList[[currentIdx]]$stopFlag == 0){ # if it has (potential) child nodes

      # distribute the training set
      trainIndex <- treeList[[currentIdx]]$splitFun(datX = datX[treeList[[currentIdx]]$idxRow,,drop = FALSE], missingReference = treeList[[currentIdx]]$misReference)

      # get child nodes
      childNodes <- lapply(seq_along(trainIndex), function(i) new_TreeeNode(datX = datX,
                                                                      response = response,
                                                                      idxCol = treeList[[currentIdx]]$idxCol,
                                                                      idxRow = treeList[[currentIdx]]$idxRow[trainIndex[[i]]],
                                                                      treeType = treeType,
                                                                      splitMethod = splitMethod,
                                                                      ldaType = ldaType,
                                                                      fastTree = fastTree,
                                                                      nodeModel = nodeModel,
                                                                      missingMethod = missingMethod,
                                                                      maxTreeLevel = maxTreeLevel,
                                                                      minNodeSize = minNodeSize,
                                                                      currentLevel = treeList[[currentIdx]]$currentLevel + 1,
                                                                      parentIndex = currentIdx))

      #> If we generate K more nodes, then the number of correctly classified sample
      #> should increase by at least K-1
      if(trainErrorCap != "none"){
        trainErrorCapNow <- ifelse(trainErrorCap == "zero", 0, length(childNodes) - 1)
        treeList[[currentIdx]]$alpha <- (treeList[[currentIdx]]$currentLoss - do.call(sum,lapply(childNodes, function(node) node$currentLoss)) - trainErrorCapNow) / (length(childNodes) - 1)
        if(treeList[[currentIdx]]$alpha < 0){
          treeList[[currentIdx]]$stopFlag = 5
          next
        }
      }
      treeList[[currentIdx]]$alpha <-  getOneSidedPvalue(N = length(treeList[[currentIdx]]$idxRow),
                                                         xBefore = treeList[[currentIdx]]$currentLoss,
                                                         xAfter = do.call(sum,lapply(childNodes, function(node) node$currentLoss)))


      #> Here we are using bootstrap sample to evaluate the model performance
      #> and decide if we are going to split more
      # nBoots <- 3
      # # tTestPvalue <- median(sapply(seq_len(nBoots), function(o_o) checkCurrentSplit(x = x[treeList[[currentIdx]]$idxRow,,drop = FALSE],
      # #                                                                response = response[treeList[[currentIdx]]$idxRow],
      # #                                                                currentNode = treeList[[currentIdx]],
      # #                                                                childNodes = childNodes,
      # #                                                                ldaType = ldaType)),na.rm = T)
      #
      # tTestPvalue <- median(sapply(seq_len(nBoots), function(o_o) generateSplitNchildren(x = x,
      #                                                                                    response = response,
      #                                                                                    idxCol = treeList[[currentIdx]]$idxCol,
      #                                                                                    idxRow = treeList[[currentIdx]]$idxRow,
      #                                                                                    treeType = treeType,
      #                                                                                    ldaType = ldaType,
      #                                                                                    fastTree = fastTree,
      #                                                                                    missingMethod = missingMethod,
      #                                                                                    splitMethod = splitMethod,
      #                                                                                    minNodeSize = minNodeSize,
      #                                                                                    bootstrap = TRUE)),na.rm = T)
      # if(is.na(tTestPvalue) | tTestPvalue >= 0.1) next

      # Put child nodes in the tree
      childIdx <- seq_along(childNodes) + length(treeList)
      treeList[[currentIdx]]$children <- childIdx
      nodeStack <- c(nodeStack, childIdx)
      for(i in seq_along(childIdx)) treeList[[childIdx[i]]] <- childNodes[[i]]

    }
  }

  for(i in seq_along(treeList)) treeList[[i]]$currentIndex <- i # assign the currentIndex

  return(treeList)
}
