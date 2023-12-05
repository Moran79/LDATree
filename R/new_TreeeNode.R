new_TreeeNode <- function(x,
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
  if(treeType == "forest") mtry <- min(50, sqrt(length(idxCol))) # for debug, suppose to be 50
  else if(fastTree) mtry <- min(50, length(idxCol)) # for debug, suppose to be 50
  else mtry <- length(idxCol)
  idxCol <- idxCol[sort(sample(seq_along(idxCol), mtry))]

  #> Remove empty levels due to partition
  xCurrent <- droplevels(x[idxRow, idxCol, drop = FALSE])
  responseCurrent <- droplevels(response[idxRow])

  #> Fix the missing values
  #> [might be changed in future due to other considerations]
  #> [if the node-wise imputation is the same as global imputation]
  imputedSummary <- missingFix(data = xCurrent, missingMethod = missingMethod)
  xCurrent <- imputedSummary$data

  #> NOTICE: If a column is constant, then it will be constant in all its subsets,
  #> so we can delete those columns in its descendents ONLY if the constant column
  #> is not from imputation, since the missing flags could be useful
  idxCurrColKeep <- constantColCheck(data = xCurrent)
  # idxCol <- idxCol[idxCurrColKeep[idxCurrColKeep <= length(idxCol)]]
  xCurrent <- xCurrent[,idxCurrColKeep, drop = FALSE]
  #> NOTICE: The missingRef should not be subset after constant check, since there
  #> are cases when the original X are constant after imputation, but its flag is important


  # Model Fitting -----------------------------------------------------------


  #> check stopping
  stopFlag <- stopCheck(responseCurrent = responseCurrent,
                        numCol = ncol(xCurrent),
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel) #  # 0/1/2: Normal/Stop+Mode/Stop+LDA


  #> Generate the split model in the current node
  # splitModel = ifelse(stopFlag == 1, "mode", "LDA")

  if(nodeModel == "LDA" | stopFlag == 0){
    if(stopFlag == 1) nodeModel <- "mode"
    else{
      #> Empty response level can not be dropped if prior exists
      datCombined = data.frame(response = responseCurrent, xCurrent)
      if(ldaType == "step") splitLDA <- nodePredict <- ldaGSVD(response~., data = datCombined, method = "step")
      else splitLDA <- nodePredict <- ldaGSVD(response~., data = datCombined, method = "all")
      resubPredict <- predict(object = nodePredict, newdata = datCombined)
    }
  }

  if(nodeModel == "mode"){
    nodePredict <- getMode(responseCurrent)
    resubPredict <- rep(nodePredict, length(responseCurrent))
  }
  currentLoss = sum(resubPredict != responseCurrent)

  # Splits Generating -----------------------------------------------------------

  #> Generate the splits
  if(stopFlag != 0) splitFun <- NULL
  else{ # if splitting goes on, find the splits
    splitFun <- getSplitFun(x = xCurrent,
                            response = responseCurrent,
                            method = splitMethod,
                            modelLDA = splitLDA)
    if(is.null(splitFun)) stopFlag <- 4 # no splits
  }


  # Final Results -----------------------------------------------------------


  currentTreeeNode <- list(
    # currentIndex = currentIndex, # will be updated in dropNodes()
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
    # pruned = NULL, # for model selection
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )
  class(currentTreeeNode) <- "TreeeNode" # Set the name for the class
  return(currentTreeeNode)
}
