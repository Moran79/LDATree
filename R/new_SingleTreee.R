new_SingleTreee <- function(datX,
                            response,
                            treeType,
                            splitMethod,
                            ldaType,
                            fastTree,
                            nodeModel,
                            missingMethod,
                            pruneMethod,
                            maxTreeLevel,
                            minNodeSize,
                            pThreshold,
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
      #> 1. check the p-value for loss drop
      lossBefore <- treeList[[currentIdx]]$currentLoss
      lossAfter <- do.call(sum,lapply(childNodes, function(node) node$currentLoss))
      treeList[[currentIdx]]$alpha <-  getOneSidedPvalue(N = length(treeList[[currentIdx]]$idxRow),
                                                         lossBefore = lossBefore,
                                                         lossAfter = lossAfter)


      #> 2. check the absolute change in drop significance
      if(trainErrorDropCap == "none"){
        trainErrorDropCapNow = -Inf
      } else if(trainErrorDropCap == "zero") {
        trainErrorDropCapNow = 0
      } else trainErrorDropCapNow = length(childNodes) - 1

      if(lossBefore - lossAfter < trainErrorDropCapNow){
        treeList[[currentIdx]]$stopFlag = 5
        next
      }


      #> 3. pre-stopping in case the tree size is going wild
      if(pruneMethod == "pre"){ # two steps ahead: two consecutive non-sig p-values are permitted
        if(treeList[[currentIdx]]$alpha < pThreshold){
          treeList[[currentIdx]]$lag = 0
        } else{ # non stat sig
          parentLag <- ifelse(treeList[[currentIdx]]$parent == 0, 0, treeList[[treeList[[currentIdx]]$parent]]$lag)
          treeList[[currentIdx]]$lag <- parentLag + 1
          if(treeList[[currentIdx]]$lag > 0){ # the number 0 could be changed to feature kStepAhead
            treeList[[currentIdx]]$stopFlag = 6
            next
          }
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
