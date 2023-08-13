#' @export
print.Treee <- function(treeeOutput){
  print(treeeOutput$treee)
  invisible(treeeOutput)
}


#' @export
print.TreeeNode <- function(treeeNode){
  cat(paste0("Node ",treeeNode$currentIndex,":\n"))
  cat("Number of observations: ", length(treeeNode$idxRow), "\n")
  # cat("Number of variables: ", length(treeeNode$idxCol), "\n")
  cat("Children:", treeeNode$children, "\n")
  cat("Node model:", treeeNode$nodeModel, "\n")
  # cat("Current Loss:", treeeNode$currentLoss, "\n")
  cat("Within-node accuracy:", treeeNode$accuracy, "\n")
  invisible(treeeNode)
}


#' @export
summary.TreeeNode <- function(treeeNode){
  # just print out everything besides some super long info
  treeeNode$idxRow <- treeeNode$idxCol <- treeeNode$nodePredict <- treeeNode$misReference <- NULL
  # class(treeeNode) <- "summary.TreeeNode"
  # return(treeeNode)
  return(unclass(treeeNode))
}
