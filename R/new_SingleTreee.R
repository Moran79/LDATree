new_SingleTreee <- function(x,
                            response,
                            idxCol,
                            idxRow,
                            splitMethod,
                            prior,
                            weights,
                            maxTreeLevel,
                            minNodeSize,
                            misClassCost,
                            missingMethod,
                            currentLevel = 0,
                            parentIndex = 0,
                            treeList = structure(list(), class = "SingleTreee")){


  # Data Cleaning -----------------------------------------------------------

  #> in-place modify xCurrentNode, pass x to its children
  #> There is no flag variables in x, but in xCurrentNode
  xCurrent <- x[idxRow, idxCol, drop = FALSE]

  #> NOTICE: If a column is constant, then it will be constant in all its subsets,
  #> so we delete those columns in its descendents. But there is no way to pass the
  #> information for constant in groups, we just leave those columns as they are

  idxCurrColKeep <- constantColCheck(data = xCurrent)
  xCurrent <- xCurrent[,idxCurrColKeep, drop = FALSE]

  # change the candidates in its children nodes, there are FLAG variable as well
  idxCol <- idxCol[idxCurrColKeep[idxCurrColKeep <= length(idxCol)]]

  # Missing Values
  imputedSummary <- missingFix(data = xCurrent, method = missingMethod)
  xCurrent <- imputedSummary$data


  # Build the Treee ---------------------------------------------------------

  currentIndex <- length(treeList) + 1 # current tree node number
  cat('The current node index is', currentIndex, '\n')

  #> check stopping
  stopFlag <- stopCheck(responseCurrent = response[idxRow],
                        idxCol = idxCol,
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel) #  # 0/1/2: Normal/Stop+Median/Stop+LDA

  # if(stopFlag != 1){
  #   #> This chunk of code is moved here, since
  #   #> 1. sometimes there are no variables left
  #   #> 2. sometimes minNodeSize is not satisfied
  #   xCurrent <- fixConstantGroup(data = xCurrent, response = response[idxRow])
  #   imputedSummary$ref <- fixReferenceWithData(data = xCurrent, ref = imputedSummary$ref)
  # }


  # Build current node
  treeList[[currentIndex]] <- new_TreeeNode(xCurrent = xCurrent,
                                            response = response,
                                            idxCol = idxCol,
                                            idxRow = idxRow,
                                            prior = prior,
                                            weights = weights,
                                            misClassCost = misClassCost,
                                            currentLevel = currentLevel,
                                            currentIndex = currentIndex,
                                            parentIndex = parentIndex,
                                            misReference = imputedSummary$ref,
                                            nodeModel = ifelse(stopFlag == 1, "mode", "LDA"))

  if (treeList[[currentIndex]]$currentLoss == 0) {stopFlag = 1} # LDA has no error

  if (stopFlag == 0) {
    splitGini <- GiniSplitScreening(xCurrent = xCurrent,
                                    response = response,
                                    idxRow = idxRow,
                                    prior = prior,
                                    minNodeSize = minNodeSize,
                                    misClassCost = misClassCost,
                                    modelLDA = treeList[[currentIndex]]$nodePredict)

    if(is.null(splitGini)) {return(treeList)} # No cut due to ties

    leftIndex <- length(treeList) + 1
    treeList <- new_SingleTreee(x = x,
                                response = response,
                                idxCol = idxCol,
                                idxRow = splitGini$left,
                                splitMethod = splitMethod,
                                prior = prior,
                                weights = weights,
                                maxTreeLevel = maxTreeLevel,
                                minNodeSize = minNodeSize,
                                misClassCost = misClassCost,
                                missingMethod = missingMethod,
                                currentLevel = currentLevel + 1,
                                parentIndex = currentIndex,
                                treeList = treeList)

    rightIndex <- length(treeList) + 1 # Right branch
    treeList <- new_SingleTreee(x = x,
                                response = response,
                                idxCol = idxCol,
                                idxRow = splitGini$right,
                                splitMethod = splitMethod,
                                prior = prior,
                                weights = weights,
                                maxTreeLevel = maxTreeLevel,
                                minNodeSize = minNodeSize,
                                misClassCost = misClassCost,
                                missingMethod = missingMethod,
                                currentLevel = currentLevel + 1,
                                parentIndex = currentIndex,
                                treeList = treeList)
    treeList[[currentIndex]]$splitCut <- splitGini$cut
    treeList[[currentIndex]]$children <- c(leftIndex, rightIndex)
  }
  return(treeList)
}
