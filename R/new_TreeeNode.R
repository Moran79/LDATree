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
  idxCol <- idxCol[idxCurrColKeep[idxCurrColKeep <= length(idxCol)]] # there are FLAGs


  # Model Fitting -----------------------------------------------------------


  #> check stopping
  stopFlag <- stopCheck(responseCurrent = responseCurrent,
                        idxCol = idxCol,
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel,
                        validSize = length(idxRowValidation)) #  # 0/1/2: Normal/Stop+Median/Stop+LDA


  #> Generate the model in the current node
  nodeModel = ifelse(stopFlag == 1, "mode", "LDA")
  if (nodeModel == "mode") {
    nodePredict <- getMode(responseCurrent)
    resubPredict <- rep(nodePredict, length(responseCurrent))
  } else if (nodeModel == "LDA") {
    #> Empty response level can not be dropped if prior exists
    datCombined = data.frame(response = responseCurrent, xCurrent)
    nodePredict <- ldaGSVD(response~., data = datCombined)
    resubPredict <- predict(object = nodePredict, newdata = datCombined)
  }


  #> Update the validation error
  if(length(idxRowValidation) == 0){ # no validation left -> too few data and we stop here
    currentLoss = 0
  }else{ # calculate the validation error
    if(nodeModel == "mode") validPredict <- rep(nodePredict, length(idxRowValidation))
    else{
      fixedData <- getDataInShape(data = xValidation[idxRowValidation,,drop = FALSE], missingReference = imputedSummary$ref)
      validPredict <- predict(object = nodePredict, newdata = fixedData)
    }
    currentLoss = sum(validPredict != responseValidation[idxRowValidation])
  }
  if(currentLoss == 0) stopFlag = 3 # validation no error, stop.


  # Splits Generating -----------------------------------------------------------

  #> Generate the splits
  if(stopFlag != 0) splitFun <- NULL
  else{ # if splitting goes on, find the splits
    splitFun <- getSplitFun(x = xCurrent,
                            response = responseCurrent,
                            method = splitMethod,
                            modelLDA = nodePredict)
    if(is.null(splitFun)) stopFlag <- 4 # no splits
  }


  # Final Results -----------------------------------------------------------


  currentTreeeNode <- list(
    # currentIndex = currentIndex, # will be updated in dropNodes()
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
    splitFun = splitFun, # save the splitting rules
    # alpha = NA, # for model selection
    # pruned = NULL, # for model selection
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )
  class(currentTreeeNode) <- "TreeeNode" # Set the name for the class
  return(currentTreeeNode)
}
