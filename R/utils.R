#' Check if tree construction should stop
#'
#' @noRd
#'
#' @param responseCurrent A vector of the current response values at the node.
#' @param numCol The number of covariate columns remaining.
#' @param maxTreeLevel The maximum allowed level of the tree.
#' @param minNodeSize The minimum number of observations required at a node.
#' @param currentLevel The current level of the tree being constructed.
stopCheck <- function(responseCurrent, numCol, maxTreeLevel, minNodeSize, currentLevel){
  flagNodeSize <- length(responseCurrent) <= minNodeSize # Data size too small
  flagTreeLevel <- currentLevel >= maxTreeLevel
  flagCol <- numCol == 0 # no covs left
  flagResponse <- length(unique(responseCurrent)) == 1 # only one class left

  if (flagResponse | flagCol | flagNodeSize) {return("Insufficient data")}
  if (flagTreeLevel) {return("Maximum level reached")}
  return("Normal")
}


#' Update Prior and Misclassification Cost
#'
#' This function updates the class prior probabilities and misclassification
#' cost matrix based on the observed response distribution. It adjusts the prior
#' and misclassification costs either inside or outside a node, depending on the
#' `insideNode` parameter.
#'
#' @noRd
updatePriorAndMisClassCost <- function(prior, misClassCost, response, insideNode = TRUE){
  if(!insideNode){ # Calculate the relative prior
    res <- folda::checkPriorAndMisClassCost(prior = prior, misClassCost = misClassCost, response = response)
    priorObs <- as.vector(table(response, dnn = NULL)) / length(response)
    res$prior <- res$prior / priorObs
  }else{
    priorObs <- table(response, dnn = NULL) / length(response)
    levelLeftIdx <- match(names(priorObs), names(prior))
    prior <- prior[levelLeftIdx] * priorObs
    res <- list(prior = prior / sum(prior),
                misClassCost = misClassCost[levelLeftIdx, levelLeftIdx, drop = FALSE])
  }
  return(res)
}


getOneSidedPvalue <- function(N, lossBefore, lossAfter){
  #> Get the p-value for testing the current split's performance
  #> H1: lossBefore > lossAfter. loss stands for the prediction error
  zStat <- (lossBefore - lossAfter) / sqrt((lossBefore * (N - lossBefore) + lossAfter * (N - lossAfter)) / N + 1e-16)
  stats::pnorm(zStat, lower.tail = FALSE)
}
