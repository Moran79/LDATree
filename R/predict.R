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
#'  probabilities, and the predicted node indices.
#' @param ... further arguments passed to or from other methods.
#'
#'@returns The function returns different values based on the `type`, if
#' * `type = 'response'`: vector of predicted responses.
#' * `type = 'prob'`: a data frame of the posterior probabilities. Each class takes a column.
#' * `type = 'all'`: a data frame contains the predicted responses, posterior probabilities, and the predicted node indices.
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
predict.Treee <- function(object, newdata, type = c("response", "prob", "all"), insideCV = FALSE, newY = NULL, ...){
  # input type: data.frame / matrix / vector
  # if(!inherits(object, "Treee")) stop("object not of class \"Treee\"")
  stopifnot(is.data.frame(newdata))

  type <- match.arg(type, c("response", "prob", "all"))

  return(predict(object$treee, newdata = newdata, type = type, insideCV = insideCV, newY = newY))
}


#' @export
predict.SingleTreee <- function(object, newdata, type = "response", insideCV = FALSE, newY = NULL, ...){
  cname <- names(object[[1]]$proportions) # find the class names
  res <- data.frame(response = character(nrow(newdata)),
                    node = numeric(nrow(newdata)),
                    newCols = matrix(0,nrow = nrow(newdata), ncol = length(cname)))
  colnames(res)[2+seq_along(cname)] <- cname # posertior probs' name
  nodeList <- vector(mode = "list", length = length(object)) # keep the testing obs index
  nodeList[[1]] <- seq_len(nrow(newdata))
  nodeStack <- c(1)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[length(nodeStack)]
    nodeStack <- nodeStack[seq_len(length(nodeStack)-1)] # last in first out
    currentNode <- object[[currentIdx]]
    currentObs <- nodeList[[currentIdx]]
    if(length(currentObs) == 0) next

    if(insideCV){
      object[[currentIdx]]$CVerror <- sum(predNode(data = newdata[currentObs,,drop = FALSE],
                                                   treeeNode = currentNode,
                                                   missingReference = currentNode$misReference,
                                                   type = "response") != newY[currentObs])
    }

    if(is.null(currentNode$children)){ # terminal nodes
      res$node[currentObs] <- currentIdx
      posteriorProbs <- predNode(data = newdata[currentObs,,drop = FALSE],
                                 treeeNode = currentNode,
                                 missingReference = currentNode$misReference,
                                 type = "prob")
      res$response[currentObs] <- colnames(posteriorProbs)[max.col(posteriorProbs, ties.method = "first")]
      res[currentObs, match(colnames(posteriorProbs), colnames(res))] <- posteriorProbs
    }else{ # internal nodes
      trainIndex <- currentNode$splitFun(datX = newdata[currentObs,,drop = FALSE], missingReference = currentNode$misReference)
      nodeStack <- c(nodeStack, currentNode$children)
      for(i in seq_along(currentNode$children)) nodeList[[currentNode$children[i]]] <- currentObs[trainIndex[[i]]]
    }
  }

  if(insideCV) return(object)
  if(type == "response") return(res$response)
  if(type == "prob") return(res[2+seq_along(cname)])
  return(res)
}





