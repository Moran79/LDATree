new_SingleTreee <- function(x,
                            response,
                            treeType,
                            ldaType,
                            missingMethod,
                            splitMethod,
                            maxTreeLevel,
                            minNodeSize,
                            trainErrorCap,
                            verbose,
                            datTest){

  treeList = structure(list(), class = "SingleTreee") # save the tree

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
                                 missingMethod = missingMethod,
                                 splitMethod = splitMethod,
                                 maxTreeLevel = maxTreeLevel,
                                 minNodeSize = minNodeSize,
                                 currentLevel = 0,
                                 parentIndex = 0)

  ### Debug ###
  # currentTestAcc <- mean(predict(treeList, datTest) == as.character(datTest[,1]))
  # splitInfo <- data.frame(splitIdx = 0,
  #                         numOfNodes = 1,
  #                         # currentTrainAcc = treeList[[1]]$accuracy,
  #                         trainAcc = treeList[[1]]$accuracy,
  #                         changeInTrainAcc = 0,
  #                         testAcc = currentTestAcc,
  #                         changeInTestAcc = 0,
  #                         # currentPvalue = treeList[[1]]$nodePredict$pValue,
  #                         # currentPillai = treeList[[1]]$nodePredict$statPillai,
  #                         # changeInOverAllPillai = treeList[[1]]$nodePredict$statPillai * length(treeList[[1]]$idxRow),
  #                         # overallPillai = treeList[[1]]$nodePredict$statPillai * length(treeList[[1]]$idxRow),
  #                         tTestStat = 0,
  #                         tTestPvalue = 0,
  #                         pointSize = 0,
  #                         accBefore = 0,
  #                         mockTtestPvalue = 0,
  #                         mockTtestPvalue2 = 0,
  #                         mockTtestPvalue3 = 0)
  ###############

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[1]; nodeStack <- nodeStack[-1] # pop the first element
    cat("The current index is:", currentIdx, "\n")
    if(treeList[[currentIdx]]$stopFlag == 0){ # if it has (potential) child nodes

      # distribute the training set
      trainIndex <- treeList[[currentIdx]]$splitFun(x = x[treeList[[currentIdx]]$idxRow,,drop = FALSE], missingReference = treeList[[currentIdx]]$misReference)

      # get child nodes
      childNodes <- lapply(seq_along(trainIndex), function(i) new_TreeeNode(x = x,
                                                                      response = response,
                                                                      idxCol = treeList[[currentIdx]]$idxCol,
                                                                      idxRow = treeList[[currentIdx]]$idxRow[trainIndex[[i]]],
                                                                      treeType = treeType,
                                                                      ldaType = ldaType,
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

      #> Here we are using bootstrap sample to evaluate the model performance
      #> and decide if we are going to split more
      nBoots <- 3
      tTestPvalue <- median(sapply(1:3, function(o_o) checkCurrentSplit(x = x[treeList[[currentIdx]]$idxRow,,drop = FALSE],
                                                                     response = response[treeList[[currentIdx]]$idxRow],
                                                                     currentNode = treeList[[currentIdx]],
                                                                     childNodes = childNodes,
                                                                     ldaType = ldaType)),na.rm = T)
      if(is.na(tTestPvalue) | tTestPvalue >= 0.1) next

      ### Debug ###
      # currentTestPredBefore <- predict(treeList, datTest, type = "all")
      # currentTestResBefore <- (currentTestPredBefore$response == as.character(datTest[,1]))[currentTestPredBefore$node == currentIdx]
      #############

      # Put child nodes in the tree
      childIdx <- seq_along(childNodes) + length(treeList)
      treeList[[currentIdx]]$children <- childIdx
      nodeStack <- c(nodeStack, childIdx)
      for(i in seq_along(childIdx)) treeList[[childIdx[i]]] <- childNodes[[i]]

      ### Debug ###
      # # lastPillai <- splitInfo$overallPillai[nrow(splitInfo)]
      # # childrenPillai <- sum(sapply(treeList[[currentIdx]]$children, function(i) length(treeList[[i]]$idxRow) * ifelse(treeList[[i]]$nodeModel == "LDA", treeList[[i]]$nodePredict$statPillai, 0)))
      # # currentPillai <- treeList[[currentIdx]]$nodePredict$statPillai * length(treeList[[currentIdx]]$idxRow)
      # currentTrainAcc <- mean(predict(treeList, x) == as.character(response))
      # currentTestPred <- predict(treeList, datTest, type = "all")
      # currentTestResAfter <- (currentTestPred$response == as.character(datTest[,1]))[currentTestPred$node %in% treeList[[currentIdx]]$children]
      # currentTestAcc <- mean(currentTestPred$response == as.character(datTest[,1]))
      # tTestRes <- t.test(currentTestResAfter, currentTestResBefore, paired = TRUE, alternative = "greater")
      #
      # splitInfo <- rbind(splitInfo, c(currentIdx,
      #                                 length(treeList),
      #                                 # treeList[[currentIdx]]$accuracy,
      #                                 currentTrainAcc,
      #                                 currentTrainAcc - splitInfo$trainAcc[nrow(splitInfo)],
      #                                 currentTestAcc,
      #                                 currentTestAcc - splitInfo$testAcc[nrow(splitInfo)],
      #                                 # treeList[[currentIdx]]$nodePredict$pValue,
      #                                 # treeList[[currentIdx]]$nodePredict$statPillai,
      #                                 # childrenPillai - currentPillai,
      #                                 # lastPillai + childrenPillai - currentPillai,
      #                                 tTestRes$statistic,
      #                                 tTestRes$p.value,
      #                                 length(currentTestResAfter),
      #                                 mean(currentTestResAfter),
      #                                 tTestPvalue,
      #                                 tTestPvalue2,
      #                                 tTestPvalue3))
      #############

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
  # attr(treeList, "summary") <- splitInfo

  return(treeList)
}
