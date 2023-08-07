#' Title
#'
#' @param object
#' @param newdata
#' @param type
#' @param ...
#'
#' @return
#' @rawNamespace S3method(predict, Treee)
#' @rawNamespace S3method(predict, SingleTreee)
#'
#' @examples
predict.Treee <- function(object, newdata, type = c("single", "grove"), ...){
  # input type: data.frame / matrix / vector
  # if(!inherits(object, "Treee")) stop("object not of class \"Treee\"")
  stopifnot(is.data.frame(newdata))

  type <- match.arg(type, c("single", "grove"))
  if(type == "grove"){
    ensembleRes <- do.call(cbind.data.frame, sapply(object$savedGrove, function(treee) predict(treee, newdata = newdata)$pred, simplify = FALSE))
    return(apply(ensembleRes,1,getMode))
  }
  return(predict(object$treee, newdata = newdata))
}

predict.SingleTreee <- function(object, newdata, ...){

  res <- data.frame(node = numeric(nrow(newdata)),
                    pred = factor(character(nrow(newdata)), levels = names(object[[1]]$posterior)))
  nodeList <- vector(mode = "list", length = length(object))
  nodeList[[1]] <- seq_len(nrow(newdata))
  nodeStack <- c(1)

  while(length(nodeStack) != 0){
    currentIdx <- nodeStack[length(nodeStack)]
    nodeStack <- nodeStack[seq_len(length(nodeStack)-1)]
    currentNode <- object[[currentIdx]]
    currentObs <- nodeList[[currentIdx]]
    # cat(currentIdx, currentObs, "\n")
    if(length(currentObs) == 0) {next} # If there is no observations in one node

    fixedData <- getDataInShape(data = newdata[currentObs,], missingReference = currentNode$misReference)

    # browser()
    if(is.null(currentNode$children)){ # terminal nodes
      # print("123")
      res$node[currentObs] <- currentIdx
      res$pred[currentObs] <- predNode(data = fixedData, treeeNode = currentNode)
    }else{
      currentScore <- getLDScores(modelLDA = currentNode$nodePredict, data = fixedData)
      leftIdx <- (currentScore <= currentNode$splitCut)
      nodeStack <- c(nodeStack, currentNode$children)
      nodeList[[currentNode$children[1]]] <- currentObs[leftIdx]
      nodeList[[currentNode$children[2]]] <- currentObs[!leftIdx]
    }
  }
  return(res)
}
