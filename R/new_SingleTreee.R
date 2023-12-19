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
                            trainErrorDropCap,
                            verbose){

  treeList = structure(list(), class = "SingleTreee") # save the tree

  #> Forest: Bootstrap samples for each single tree
  if(treeType == "forest"){
    idxRowRoot <- sample(length(response), length(response), replace = TRUE)
  } else idxRowRoot <- seq_len(nrow(datX))


  ### Initialize the first Node ###

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
    if(verbose) cat("The current index is:", currentIdx, "\n")

    if(treeList[[currentIdx]]$stopFlag == 0){ # if it has (potential) child nodes

      trainIndex <- attr(treeList[[currentIdx]]$splitFun, "splitRes") # distribute the training set


      ### Get child nodes ###

      childNodes <- lapply(seq_along(trainIndex),
                           function(i) new_TreeeNode(datX = datX,
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


      ### Stopping & pruning ###

      lossBefore <- treeList[[currentIdx]]$currentLoss
      lossAfter <- do.call(sum,lapply(childNodes, function(node) node$currentLoss))
      treeList[[currentIdx]]$alpha <-  getOneSidedPvalue(N = length(treeList[[currentIdx]]$idxRow),
                                                         lossBefore = lossBefore,
                                                         lossAfter = lossAfter)
      if(trainErrorDropCap != "none"){
        trainErrorDropCapNow <- ifelse(trainErrorDropCap == "zero", 0, length(childNodes) - 1)
        if(lossBefore - lossAfter < trainErrorDropCapNow){
          treeList[[currentIdx]]$stopFlag = 5
          next
        }
      }


      ### Put child nodes in the tree ###

      childIdx <- seq_along(childNodes) + length(treeList)
      treeList[[currentIdx]]$children <- childIdx
      nodeStack <- c(nodeStack, childIdx)
      for(i in seq_along(childIdx)) treeList[[childIdx[i]]] <- childNodes[[i]]
    }
  }

  for(i in seq_along(treeList)) treeList[[i]]$currentIndex <- i # assign the currentIndex

  return(treeList)
}
