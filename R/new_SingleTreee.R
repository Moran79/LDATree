new_SingleTreee <- function(x,
                            response,
                            treeType,
                            ldaType,
                            forest,
                            missingMethod,
                            splitMethod,
                            maxTreeLevel,
                            minNodeSize,
                            trainErrorCap,
                            verbose){

  treeList = structure(list(), class = "SingleTreee") # save the tree

  ### Debug ###
  splitInfo <- data.frame(splitIdx = 0,
                          numOfNodes = 0,
                          trainAcc = 0,
                          testAcc = 0)

  # Bootstrap sample in ensemble method
  if(treeType == "forest") idxRowRoot <- sample(length(response), length(response), replace = TRUE)
  else idxRowRoot <- seq_len(nrow(x))

  # initialize the first Node
  nodeStack <- c(1)
  treeList[[1]] <- new_TreeeNode(x = x,
                                 response = response,
                                 idxCol = seq_len(ncol(x)),
                                 idxRow = idxRowRoot,
                                 treeType = treeType,
                                 ldaType = ldaType,
                                 forest = forest,
                                 missingMethod = missingMethod,
                                 splitMethod = splitMethod,
                                 maxTreeLevel = maxTreeLevel,
                                 minNodeSize = minNodeSize,
                                 currentLevel = 0,
                                 parentIndex = 0)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[1]; nodeStack <- nodeStack[-1] # pop the first element
    if(treeList[[currentIdx]]$stopFlag == 0){ # if splitting goes on

      # distribute the training set
      trainIndex <- treeList[[currentIdx]]$splitFun(x = x[treeList[[currentIdx]]$idxRow,,drop = FALSE], missingReference = treeList[[currentIdx]]$misReference)

      # get child nodes
      childNodes <- lapply(seq_along(trainIndex), function(i) new_TreeeNode(x = x,
                                                                      response = response,
                                                                      idxCol = treeList[[currentIdx]]$idxCol,
                                                                      idxRow = treeList[[currentIdx]]$idxRow[trainIndex[[i]]],
                                                                      treeType = treeType,
                                                                      ldaType = ldaType,
                                                                      forest = forest,
                                                                      missingMethod = missingMethod,
                                                                      splitMethod = splitMethod,
                                                                      maxTreeLevel = maxTreeLevel,
                                                                      minNodeSize = minNodeSize,
                                                                      currentLevel = treeList[[currentIdx]]$currentLevel + 1,
                                                                      parentIndex = currentIdx))

      #> If we generate K more nodes, then the number of correctly classified sample
      #> should increase by at least K-1
      if(trainErrorCap != "none"){
        trainErrorCapNow <- ifelse(trainErrorCap == "zero", 0, length(childNodes) - 1)
        treeList[[currentIdx]]$alpha <- (treeList[[currentIdx]]$currentLoss - do.call(sum,lapply(childNodes, function(node) node$currentLoss)) - trainErrorCapNow) / (length(childNodes) - 1)
        if(treeList[[currentIdx]]$alpha <= 0){
          treeList[[currentIdx]]$stopFlag = 5
          next
        }
      }


      ### Debug ###
      # splitInfo <- rbind(splitInfo, c(nrow(splitInfo), length(treeList),
      #                                 mean(predict(treeList, dat_trainC) == dat_trainC[,1]),
      #                                 mean(predict(treeList, dat_testC) == dat_testC[,1])))
      #############

      # Put child nodes in the tree
      childIdx <- seq_along(childNodes) + length(treeList)
      treeList[[currentIdx]]$children <- childIdx
      nodeStack <- c(nodeStack, childIdx)
      for(i in seq_along(childIdx)) treeList[[childIdx[i]]] <- childNodes[[i]]
    }
  }

  # Update the tree, remove those splits which include the training errors
  if(trainErrorCap != "none"){
    treeList <- makeAlphaMono(treeeList = treeList, trainErrorCap = trainErrorCap)
    treeList <- pruneTreee(treeeList = treeList, trainErrorCap = trainErrorCap, alpha = -0.5)
    treeList <- dropNodes(treeList)
  }else{
    for(i in seq_along(treeList)){
      treeList[[i]]$currentIndex <- i # re-assign the currentIndex
    }
  }

  ### DEBUG ###

  attr(treeList, "summary") <- splitInfo

  return(treeList)
}
