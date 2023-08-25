new_SingleTreee <- function(x,
                            response,
                            idxCol,
                            idxRow,
                            missingMethod,
                            splitMethod,
                            maxTreeLevel,
                            minNodeSize,
                            verbose,
                            currentLevel = 0,
                            parentIndex = 0,
                            treeList = structure(list(), class = "SingleTreee")){

  #> Some notes to clarify the parameters...

  # Data Cleaning -----------------------------------------------------------

  # Remove empty levels due to partition
  xCurrent <- droplevels(x[idxRow, idxCol, drop = FALSE])
  responseCurrent <- droplevels(response[idxRow])
  #> Notes: xCurrent / responseCurrent are ephemeral and
  #> will be removed once the current node is completed

  # Fix the missing values
  imputedSummary <- missingFix(data = xCurrent, missingMethod = missingMethod)
  xCurrent <- imputedSummary$data

  #> NOTICE: If a column is constant, then it will be constant in all its subsets,
  #> so we delete those columns in its descendents.

  idxCurrColKeep <- constantColCheck(data = xCurrent)
  xCurrent <- xCurrent[,idxCurrColKeep, drop = FALSE]

  # change the candidates in its children nodes, there are FLAG variable as well
  idxCol <- idxCol[idxCurrColKeep[idxCurrColKeep <= length(idxCol)]]


  # Build the Treee ---------------------------------------------------------

  currentIndex <- length(treeList) + 1 # current tree node number
  # cat('The current node index is', currentIndex, '\n')

  #> check stopping
  stopFlag <- stopCheck(responseCurrent = responseCurrent,
                        idxCol = idxCol,
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel) #  # 0/1/2: Normal/Stop+Median/Stop+LDA

  # Build current node
  treeList[[currentIndex]] <- new_TreeeNode(xCurrent = xCurrent,
                                            responseCurrent = responseCurrent,
                                            idxCol = idxCol,
                                            idxRow = idxRow,
                                            currentLevel = currentLevel,
                                            currentIndex = currentIndex,
                                            parentIndex = parentIndex,
                                            misReference = imputedSummary$ref,
                                            nodeModel = ifelse(stopFlag == 1, "mode", "LDA"))

  # if (treeList[[currentIndex]]$currentLoss == 0) stopFlag = 1 # LDA has no error

  if (stopFlag == 0) {
    splitGini <- GiniSplitScreening(xCurrent = xCurrent,
                                    responseCurrent = responseCurrent,
                                    idxRow = idxRow,
                                    minNodeSize = minNodeSize,
                                    modelLDA = treeList[[currentIndex]]$nodePredict)

    if(is.null(splitGini)) return(treeList) # No cut due to ties

    leftIndex <- length(treeList) + 1
    treeList <- new_SingleTreee(x = x,
                                response = response,
                                idxCol = idxCol,
                                idxRow = splitGini$left,
                                splitMethod = splitMethod,
                                maxTreeLevel = maxTreeLevel,
                                minNodeSize = minNodeSize,
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
                                maxTreeLevel = maxTreeLevel,
                                minNodeSize = minNodeSize,
                                missingMethod = missingMethod,
                                currentLevel = currentLevel + 1,
                                parentIndex = currentIndex,
                                treeList = treeList)
    treeList[[currentIndex]]$splitCut <- splitGini$cut
    treeList[[currentIndex]]$children <- c(leftIndex, rightIndex)
  }
  return(treeList)
}
