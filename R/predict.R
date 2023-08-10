#' Title
#'
#' @param object
#' @param newdata
#' @param type
#'
#' @return
#' @export
#'
#' @examples
predict.Treee <- function(object, newdata, type = c("response", "prob", "all", "grove")){
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

#' Title
#'
#' @param object
#' @param newdata
#' @param type
#'
#' @return
#' @export
#'
#' @examples
predict.SingleTreee <- function(object, newdata, type = "response"){

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
