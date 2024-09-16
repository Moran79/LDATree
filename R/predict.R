#' Predictions From a Fitted Treee Object
#'
#' Generate predictions on new data using a fitted `Treee` model.
#'
#' @param object A fitted model object of class `Treee`, typically the result of
#'   the [Treee()] function.
#' @param newdata A data frame containing the predictor variables. Missing
#'   values are allowed and will be handled according to the fitted tree's
#'   method for handling missing data.
#' @param type A character string specifying the type of prediction to return.
#'   Options are:
#'  * `'response'`: returns the predicted class for each observation (default).
#'  * `'prob'`: returns a data frame of posterior probabilities for each class.
#'  * `'all'`: returns a data frame containing predicted classes, posterior probabilities, and the predicted node indices.
#' @param ... Additional arguments passed to or from other methods.
#'
#' @return Depending on the value of `type`, the function returns:
#'  * If `type = 'response'`: A character vector of predicted class labels.
#'  * If `type = 'prob'`: A data frame of posterior probabilities, where each class has its own column.
#'  * If `type = 'all'`: A data frame containing predicted class labels, posterior probabilities, and the predicted node indices.
#'
#'   Note: For factor predictors, if a level not present in the training data is
#'   found in `newdata`, it will be treated as missing and handled according to
#'   the `missingMethod` specified in the fitted tree.
#'
#' @export
#'
#' @examples
#' fit <- Treee(datX = iris[, -5], response = iris[, 5], verbose = FALSE)
#' head(predict(fit, iris)) # Predicted classes
#' head(predict(fit, iris[, -5], type = "prob")) # Posterior probabilities
#' head(predict(fit, iris[, -5], type = "all")) # Full details
predict.Treee <- function(object, newdata, type = c("response", "prob", "all"), ...){
  if (!is.data.frame(newdata)) stop("datX must be a data.frame")

  dots <- list(...)
  insideCV <- if (!is.null(dots$insideCV)) dots$insideCV else FALSE
  obsY <- if (!is.null(dots$obsY)) dots$obsY else NULL

  type <- match.arg(type, c("response", "prob", "all"))
  cname <- names(object[[1]]$proportions) # find the class names
  res <- data.frame(response = character(nrow(newdata)),
                    node = numeric(nrow(newdata)),
                    newCols = matrix(0, nrow = nrow(newdata), ncol = length(cname)))
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
      object[[currentIdx]]$CVerror <- sum(predNode(datX = newdata[currentObs,,drop = FALSE],
                                                   treeeNode = currentNode,
                                                   type = "response") != obsY[currentObs])
    }

    if(is.null(currentNode$children)){ # terminal nodes
      res$node[currentObs] <- currentIdx
      posteriorProbs <- predNode(datX = newdata[currentObs,,drop = FALSE],
                                 treeeNode = currentNode,
                                 type = "prob")
      res$response[currentObs] <- colnames(posteriorProbs)[max.col(-posteriorProbs %*% t(currentNode$misClassCost), ties.method = "first")]
      res[currentObs, match(colnames(posteriorProbs), colnames(res))] <- posteriorProbs
    }else{ # internal nodes
      trainIndex <- currentNode$splitFun(datX = newdata[currentObs,,drop = FALSE])
      nodeStack <- c(nodeStack, currentNode$children)
      for(i in seq_along(currentNode$children)) nodeList[[currentNode$children[i]]] <- currentObs[trainIndex[[i]]]
    }
  }

  if(insideCV) return(object)
  if(type == "response") return(res$response)
  if(type == "prob") return(res[2+seq_along(cname)])
  return(res)
}


predNode <- function(datX, treeeNode, type){
  if(treeeNode$nodeModel != "ULDA"){
    if(type == "response"){
      pred <- rep(treeeNode$nodePredict, nrow(datX))
    } else{ # if type = "all", the extra response column will be added later
      pred <- matrix(0,nrow = nrow(datX), ncol = length(treeeNode$proportions), dimnames = list(c(), names(treeeNode$proportions)))
      pred[,which(treeeNode$nodePredict == colnames(pred))] <- 1
    }
  } else pred <- predict(object = treeeNode$nodePredict, newdata = datX, type = type)
  return(pred)
}




