new_TreeeNode <- function(datX,
                          response,
                          idxCol,
                          idxRow,
                          ldaType,
                          nodeModel,
                          missingMethod,
                          prior,
                          maxTreeLevel,
                          minNodeSize,
                          currentLevel,
                          parentIndex) {


  # Data Cleaning -----------------------------------------------------------

  #> Remove empty levels due to partition
  xCurrent <- droplevels(datX[idxRow, idxCol, drop = FALSE])
  responseCurrent <- droplevels(response[idxRow])

  #> Fix the missing values
  imputedSummary <- missingFix(data = xCurrent, missingMethod = missingMethod)
  xCurrent <- imputedSummary$data

  #> NOTICE: The missingRef should not be subset after constant check, since there
  #> are cases when the original X are constant after imputation, but its flag is important
  idxCurrColKeep <- constantColCheck(data = xCurrent)
  xCurrent <- xCurrent[, idxCurrColKeep, drop = FALSE]


  # Model Fitting -----------------------------------------------------------

  #> check stopping
  stopFlag <- stopCheck(responseCurrent = responseCurrent,
                        numCol = ncol(xCurrent),
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel) #  # 0/1/2: Normal/Stop+Mode/Stop+LDA

  #> Based on the node model, decide whether we should fit LDA
  if(nodeModel == "LDA" | stopFlag == 0){ # LDA model, or mode model with LDA splits
    if(stopFlag == 1){ # when LDA can not be fitted
      nodeModel <- "mode"
    } else{
      #> Empty response level can not be dropped if prior exists
      splitLDA <- nodePredict <- ldaGSVD(datX = xCurrent, response = responseCurrent, method = ldaType, fixNA = FALSE, prior = prior, insideTree = TRUE)
      resubPredict <- predict(object = nodePredict, newdata = xCurrent)
    }
  }

  if(nodeModel == "mode"){
    nodePredict <- getMode(responseCurrent, prior = prior)
    resubPredict <- rep(nodePredict, length(responseCurrent))
  }
  currentLoss = sum(resubPredict != responseCurrent) # save the currentLoss for future accuracy calculation

  #> if not as good as mode, change it to mode,
  #> but the splitting goes on, since the next split might be better.
  #> The code is subjective to change if prior will be added
  if(currentLoss >= length(responseCurrent) - max(unname(table(responseCurrent)))){
    nodeModel <- "mode"
    nodePredict <- getMode(responseCurrent, prior = prior)
    resubPredict <- rep(nodePredict, length(responseCurrent))
    currentLoss = sum(resubPredict != responseCurrent)
  }


  # Splits Generating -----------------------------------------------------------

  #> Generate the splits
  if(stopFlag == 0){ # if splitting goes on, find the splits
    splitFun <- getSplitFunLDA(datX = xCurrent,
                               response = responseCurrent,
                               modelLDA = splitLDA)
    if(is.null(splitFun)) stopFlag <- 4 # no splits
  } else splitFun <- NULL


  # Final Results -----------------------------------------------------------

  currentTreeeNode <- list(
    # currentIndex = currentIndex, # will be updated in new_SingleTreee()
    currentLevel = currentLevel,
    idxCol = idxCol,
    idxRow = idxRow,
    currentLoss = currentLoss, # this loss should account for sample size
    accuracy = 1 - currentLoss / length(responseCurrent),
    stopFlag = stopFlag,
    proportions = table(responseCurrent, dnn = NULL), # remove the name of the table
    parent = parentIndex,
    children = c(), # is.null to check terminal nodes
    misReference = imputedSummary$ref,
    splitFun = splitFun, # save the splitting rules
    # alpha = NA, # p-value from t-test, to measure the split's strength for model selection
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )
  class(currentTreeeNode) <- "TreeeNode" # Set the name for the class
  return(currentTreeeNode)
}
