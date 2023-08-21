new_TreeeNode <- function(xCurrent,
                          responseCurrent,
                          idxCol,
                          idxRow,
                          currentLevel,
                          currentIndex,
                          parentIndex,
                          misReference,
                          nodeModel) {
  if (nodeModel == "mode") {
    nodePredict <- getMode(responseCurrent)
    resubPredict <- rep(nodePredict, length(responseCurrent))
  } else if (nodeModel == "LDA") {
    #> Empty response level can not be dropped if prior exists
    datCombined = data.frame(response = responseCurrent, xCurrent)
    nodePredict <- ldaGSVD(response~., data = datCombined)
    resubPredict <- predict(object = nodePredict, newdata = datCombined)
  }

  currentTreeeNode <- list(
    currentIndex = currentIndex,
    currentLevel = currentLevel,
    idxRow = idxRow,
    idxCol = idxCol,
    currentLoss = sum(resubPredict != responseCurrent), # this loss should account for sample size
    accuracy = mean(resubPredict == responseCurrent),
    proportions = table(responseCurrent, dnn = NULL), # remove the name of the table
    parent = parentIndex,
    children = c(), # is.null to check terminal nodes
    misReference = misReference,
    splitCut = NA, # Splitting criteria
    # offsprings = c(), # all terminal nodes
    # alpha = NA, # for CART pruning
    # offspringLoss = NA, # sum of currentLoss of all its offsprings
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )

  # Set the name for the class
  class(currentTreeeNode) <- "TreeeNode"
  return(currentTreeeNode)
}
