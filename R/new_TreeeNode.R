new_TreeeNode <- function(datX,
                          response,
                          idxCol,
                          idxRow,
                          treeType,
                          splitMethod,
                          ldaType,
                          fastTree,
                          nodeModel,
                          missingMethod,
                          maxTreeLevel,
                          minNodeSize,
                          currentLevel,
                          parentIndex) {


  # Data Cleaning -----------------------------------------------------------

  #> First, we subset the columns for forest / fastTree
  #> The only problem is that the number mtry calculated is based
  #> on the data.frame, not the design matrix
  if(treeType == "forest"){
    mtry <- min(50, ceiling(sqrt(length(idxCol))))
  } else if(fastTree){
    mtry <- min(50, length(idxCol))
  } else mtry <- length(idxCol)
  idxCol <- idxCol[sort(sample(seq_along(idxCol), mtry))]

  #> Remove empty levels due to partition
  xCurrent <- droplevels(datX[idxRow, idxCol, drop = FALSE])
  responseCurrent <- droplevels(response[idxRow])

  #> Fix the missing values
  #> [might be changed in future due to other considerations]
  #> [if the node-wise imputation is the same as global imputation]
  imputedSummary <- missingFix(data = xCurrent, missingMethod = missingMethod)
  xCurrent <- imputedSummary$data

  #> NOTICE: If a column is constant, then it will be constant in all its subsets,
  #> and we can delete those columns in its descendants ONLY when the constant column
  #> is not from imputation, since the missing flags could be useful
  idxCurrColKeep <- constantColCheck(data = xCurrent)
  # idxCol <- idxCol[idxCurrColKeep[idxCurrColKeep <= length(idxCol)]]
  xCurrent <- xCurrent[, idxCurrColKeep, drop = FALSE]
  #> NOTICE: The missingRef should not be subset after constant check, since there
  #> are cases when the original X are constant after imputation, but its flag is important


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
      datCombined = data.frame(response = responseCurrent, xCurrent)

      if(ldaType == "step"){
        splitLDA <- nodePredict <- ldaGSVD(response~., data = datCombined, method = "step")
      } else splitLDA <- nodePredict <- ldaGSVD(response~., data = datCombined, method = "all")
      resubPredict <- predict(object = nodePredict, newdata = datCombined)
    }
  }

  if(nodeModel == "mode"){
    nodePredict <- getMode(responseCurrent)
    resubPredict <- rep(nodePredict, length(responseCurrent))
  }
  currentLoss = sum(resubPredict != responseCurrent) # save the currentLoss for future accuracy calculation

  # if not as good as mode, change it to mode
  # subject to change if prior will be added
  if(currentLoss >= length(responseCurrent) - max(unname(table(responseCurrent)))){
    nodeModel <- "mode"
    nodePredict <- getMode(responseCurrent)
    resubPredict <- rep(nodePredict, length(responseCurrent))
    currentLoss = sum(resubPredict != responseCurrent)
  }


  # Splits Generating -----------------------------------------------------------

  #> Generate the splits
  if(stopFlag == 0){ # if splitting goes on, find the splits
    splitFun <- getSplitFun(datX = xCurrent,
                            response = responseCurrent,
                            method = splitMethod,
                            modelLDA = splitLDA)
    if(is.null(splitFun)) stopFlag <- 4 # no splits
  } else splitFun <- NULL


  # Final Results -----------------------------------------------------------

  currentTreeeNode <- list(
    # currentIndex = currentIndex, # will be updated in new_SingleTreee() & pruneByTrainErrDrop()
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
    # alpha = NA, # for model selection
    # lag = NA, # for two-steps ahead
    # pruned = NULL, # for model selection
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )
  class(currentTreeeNode) <- "TreeeNode" # Set the name for the class
  return(currentTreeeNode)
}
