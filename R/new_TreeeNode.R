new_TreeeNode <- function(x,
                          response,
                          xValidation,
                          responseValidation,
                          idxCol,
                          idxRow,
                          idxRowValidation,
                          missingMethod,
                          splitMethod,
                          maxTreeLevel,
                          minNodeSize,
                          currentLevel,
                          parentIndex) {


  # Data Cleaning -----------------------------------------------------------

  # Remove empty levels due to partition
  xCurrent <- droplevels(x[idxRow, idxCol, drop = FALSE])
  responseCurrent <- droplevels(response[idxRow])

  # Fix the missing values
  imputedSummary <- missingFix(data = xCurrent, missingMethod = missingMethod)
  xCurrent <- imputedSummary$data

  #> NOTICE: If a column is constant, then it will be constant in all its subsets,
  #> so we delete those columns in its descendents.

  idxCurrColKeep <- constantColCheck(data = xCurrent)
  xCurrent <- xCurrent[,idxCurrColKeep, drop = FALSE]

  # change the candidates in its children nodes, there are FLAG variable as well
  idxCol <- idxCol[idxCurrColKeep[idxCurrColKeep <= length(idxCol)]]

  #> check stopping
  stopFlag <- stopCheck(responseCurrent = responseCurrent,
                        idxCol = idxCol,
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel,
                        validSize = length(idxRowValidation)) #  # 0/1/2: Normal/Stop+Median/Stop+LDA
  nodeModel = ifelse(stopFlag == 1, "mode", "LDA")

  if (nodeModel == "mode") {
    nodePredict <- getMode(responseCurrent)
    resubPredict <- rep(nodePredict, length(responseCurrent))
    validPredict <- rep(nodePredict, length(idxRowValidation))
  } else if (nodeModel == "LDA") {
    #> Empty response level can not be dropped if prior exists
    datCombined = data.frame(response = responseCurrent, xCurrent)
    nodePredict <- ldaGSVD(response~., data = datCombined)
    resubPredict <- predict(object = nodePredict, newdata = datCombined)

    if(length(idxRowValidation) != 0){ # validation set has at least one obs
      fixedData <- getDataInShape(data = xValidation[idxRowValidation,,drop = FALSE], missingReference = imputedSummary$ref)
      validPredict <- predict(object = nodePredict, newdata = fixedData)
    }else{validPredict <- numeric()}
  }

  # validation no error, stop.
  currentLoss = sum(validPredict != responseValidation[idxRowValidation])
  if(stopFlag == 0 & currentLoss == 0){
    stopFlag = 3
  }


  if(stopFlag == 0){ # if splitting goes on
    # find the splits
    splitGini <- GiniSplitScreening(xCurrent = xCurrent,
                                    responseCurrent = responseCurrent,
                                    idxRow = idxRow,
                                    minNodeSize = minNodeSize,
                                    modelLDA = nodePredict)
  }else{splitGini <- NULL}


  currentTreeeNode <- list(
    # currentIndex = currentIndex,
    currentLevel = currentLevel,
    idxCol = idxCol,
    idxRow = idxRow,
    idxRowValidation = idxRowValidation,
    currentLoss = currentLoss, # this loss should account for sample size
    accuracy = mean(resubPredict == responseCurrent),
    lag = 0, # count how many negative alphas in its ancestors
    stopFlag = stopFlag,
    proportions = table(responseCurrent, dnn = NULL), # remove the name of the table
    parent = parentIndex,
    children = c(), # is.null to check terminal nodes
    misReference = imputedSummary$ref,
    splitGini = splitGini, # save the Gini Split
    splitCut = NA, # Splitting criteria
    # offsprings = c(), # all terminal nodes
    # alpha = NA, # for model selection
    # offspringLoss = NA, # sum of currentLoss of all its offsprings
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )

  # Set the name for the class
  class(currentTreeeNode) <- "TreeeNode"
  return(currentTreeeNode)
}
