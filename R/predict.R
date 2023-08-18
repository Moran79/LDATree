#'Predictions from a fitted Treee object
#'
#'Prediction of test data using a fitted Treee object
#'
#'@param object a fitted model object of class `Treee`, which is assumed to be
#'  the result of the [Treee()] function.
#'@param newdata data frame containing the values at which predictions are
#'  required. Missing values are allowed.
#'@param type character string denoting the type of predicted value returned.
#'  The default is to return the predicted class (`type = 'response'`). The
#'  predicted posterior probabilities for each class will be returned if `type =
#'  'prob'`. `'all'` returns a data frame with predicted classes, posterior
#'  probabilities, and the predicted node indices. If cross-validation is
#'  carried out during the LDATree fitting, `'grove'` option is available and
#'  will output an ensemble result from `k` LDATrees where `k` is the number of
#'  cross-validation.
#' @param ... further arguments passed to or from other methods.
#'
#'@returns The function returns different values based on the `type`, if
#' * `type = 'response'`: vector of predicted responses.
#' * `type = 'prob'`: a data frame of the posterior probabilities. Each class takes a column.
#' * `type = 'all'`: a data frame contains the predicted responses, posterior probabilities, and the predicted node indices.
#' * `type = 'grove'`: vector of predicted responses using
#'  the ensemble method. Only available when cross-validation is carried out
#'  during the tree generating process.
#'
#'Note: for factor predictors, if it contains a level which is not used to
#'  grow the tree, it will be converted to missing and will be imputed according
#'  to the `missingMethod` in the fitted tree.
#'@export
#'
#' @examples
#' fit <- Treee(Species~., data = iris)
#' predict(fit,iris)
#' # output prosterior probabilities
#' predict(fit,iris,type = "prob")
predict.Treee <- function(object, newdata, type = c("response", "prob", "all", "grove"), ...){
  # input type: data.frame / matrix / vector
  # if(!inherits(object, "Treee")) stop("object not of class \"Treee\"")
  stopifnot(is.data.frame(newdata))

  type <- match.arg(type, c("response", "prob", "all", "grove"))
  if(type == "grove"){
    ensembleRes <- do.call(cbind.data.frame, sapply(object$savedGrove, function(treee) predict(treee, newdata = newdata), simplify = FALSE))
    return(apply(ensembleRes,1,getMode))
  }
  return(predict(object$treee, newdata = newdata, type = type))
}


#' @export
predict.SingleTreee <- function(object, newdata, type = "response", ...){

  cname <- names(object[[1]]$proportions)
  res <- data.frame(response = character(nrow(newdata)),
                    node = numeric(nrow(newdata)),
                    newCols = matrix(0,nrow = nrow(newdata), ncol = length(cname)))
  colnames(res)[2+seq_along(cname)] <- cname

  nodeList <- vector(mode = "list", length = length(object))
  nodeList[[1]] <- seq_len(nrow(newdata))
  nodeStack <- c(1)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[length(nodeStack)]
    nodeStack <- nodeStack[seq_len(length(nodeStack)-1)]
    currentNode <- object[[currentIdx]]
    currentObs <- nodeList[[currentIdx]]

    if(length(currentObs) == 0) {next} # If there is no observations in one node

    fixedData <- getDataInShape(data = newdata[currentObs,,drop = FALSE], missingReference = currentNode$misReference)

    if(is.null(currentNode$children)){ # terminal nodes
      res$node[currentObs] <- currentIdx
      res$response[currentObs] <- predNode(data = fixedData, treeeNode = currentNode, type = "response")
      posteriorProbs <- predNode(data = fixedData, treeeNode = currentNode, type = "prob")
      # browser()
      res[currentObs, match(colnames(posteriorProbs), colnames(res))] <- posteriorProbs
    }else{
      currentScore <- getLDscores(modelLDA = currentNode$nodePredict, data = fixedData, nScores = 1)
      leftIdx <- (currentScore <= currentNode$splitCut)
      nodeStack <- c(nodeStack, currentNode$children)
      nodeList[[currentNode$children[1]]] <- currentObs[leftIdx]
      nodeList[[currentNode$children[2]]] <- currentObs[!leftIdx]
    }
  }
  if(type == "response") return(res$response)
  if(type == "prob") return(res[2+seq_along(cname)])
  return(res)
}
