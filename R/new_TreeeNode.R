#' Create a New Tree Node in the Decision Tree
#'
#' This function creates a new node for the decision tree by fitting a model
#' (such as ULDA or a mode model) based on the data at the current node. It
#' checks for stopping conditions, fits the model, and generates splits if
#' necessary.
#'
#' @noRd
new_TreeeNode <- function(datX,
                          response,
                          idxCol,
                          idxRow,
                          ldaType,
                          nodeModel,
                          maxTreeLevel,
                          minNodeSize,
                          prior,
                          misClassCost,
                          missingMethod,
                          kSample,
                          currentLevel,
                          parentIndex) {

  # Data Cleaning -----------------------------------------------------------

  #> Remove empty levels due to partition
  xCurrent <- droplevels(datX[idxRow, idxCol, drop = FALSE])
  responseCurrent <- droplevels(response[idxRow])
  priorAndMisClassCost <- updatePriorAndMisClassCost(prior = prior, misClassCost = misClassCost, response = responseCurrent, insideNode = TRUE)
  prior <- priorAndMisClassCost$prior; misClassCost <- priorAndMisClassCost$misClassCost

  # Model Fitting -----------------------------------------------------------

  stopInfo <- stopCheck(responseCurrent = responseCurrent,
                        numCol = ncol(xCurrent),
                        maxTreeLevel = maxTreeLevel,
                        minNodeSize = minNodeSize,
                        currentLevel = currentLevel) # Normal/Stop+Mode/Stop+ULDA

  #> Based on the node model, decide whether we should fit ULDA
  if(nodeModel == "ULDA" | stopInfo == "Normal"){ # ULDA model, or mode model with ULDA splits
    if(stopInfo == "Insufficient data"){ # when LDA can not be fitted
      nodeModel <- "mode"
    } else{
      splitLDA <- nodePredict <- tryCatch({
        folda::folda(datX = xCurrent,
                     response = responseCurrent,
                     subsetMethod = ldaType,
                     prior = prior,
                     misClassCost = misClassCost,
                     missingMethod = missingMethod,
                     downSampling = (kSample != -1),
                     kSample = kSample)
      }, error = function(e) {NULL})

      if(!is.null(nodePredict)){
        resubPredict <- predict(object = nodePredict, newdata = xCurrent)
        currentLoss = sum(resubPredict != responseCurrent) # save the currentLoss for future accuracy calculation
        #> if not as good as mode, change it to mode,
        #> but the splitting goes on, since the next split might be better.
        if(currentLoss >= length(responseCurrent) - max(table(responseCurrent))) nodeModel <- "mode"
      } else{
        nodeModel <- "mode"
        stopInfo <- "Insufficient data"
      }
    }
  }

  if(nodeModel == "mode"){
    nodePredict <- folda::getMode(responseCurrent, prior = prior)
    resubPredict <- rep(nodePredict, length(responseCurrent))
    currentLoss = sum(resubPredict != responseCurrent)
  }

  # Splits Generating -----------------------------------------------------------

  if(stopInfo == "Normal"){ # if splitting goes on, find the splits
    splitFun <- getSplitFunLDA(datX = xCurrent,
                               modelULDA = splitLDA)
    if(is.null(splitFun)) stopInfo <- "No feasible splits"
  } else splitFun <- NULL

  # Final Results -----------------------------------------------------------

  currentTreeeNode <- list(
    currentLevel = currentLevel,
    idxCol = idxCol,
    idxRow = idxRow,
    currentLoss = currentLoss, # this loss should account for sample size
    accuracy = 1 - currentLoss / length(responseCurrent),
    stopInfo = stopInfo,
    proportions = table(responseCurrent, dnn = NULL), # remove the name of the table
    prior = prior,
    misClassCost = misClassCost,
    parent = parentIndex,
    children = c(), # is.null to check terminal nodes
    splitFun = splitFun, # save the splitting rules
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
    # currentIndex = currentIndex, # will be updated in new_SingleTreee()
    # alpha = NA, # p-value from t-test, to measure the split's strength for model selection
    # pruned = NULL # generated during pruning
  )
  class(currentTreeeNode) <- "TreeeNode" # Set the name for the class
  return(currentTreeeNode)
}
