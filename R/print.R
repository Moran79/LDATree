#' Title
#'
#' @param treeeList
#'
#' @return
#' @rawNamespace S3method(print, SingleTreee)
#' @rawNamespace S3method(print, TreeeNode)
#' @rawNamespace S3method(print, Treee)
#'
#' @examples
print.SingleTreee <- function(treeeList){
  cat("The current tree has",length(treeeList),"node(s)\n")
  for(i in seq_along(treeeList)){
    treeeNode <- treeeList[[i]]
    cat("\n")
    message(paste0("Node ",treeeNode$currentIndex,":"))
    cat("Number of observations: ", length(treeeNode$idxRow), "\n")
    cat("Number of variables: ", length(treeeNode$idxCol), "\n")
    cat("Children:", treeeNode$children, "\n")
    cat("Node model:", treeeNode$nodeModel, "\n")
    cat("Current Loss:", treeeNode$currentLoss, "\n")
    cat("Accuracy:", treeeNode$accuracy, "\n")
  }
  invisible(treeeList)
}

print.TreeeNode <- function(treeeNode){
  treeeNodePrint <- treeeNode
  # Hide some page-long information
  treeeNodePrint$idxRow <- treeeNodePrint$idxCol <- treeeNodePrint$nodePredict <- treeeNodePrint$levelReference <- NULL
  print(unclass(treeeNodePrint))
  invisible(treeeNode)
}

print.Treee <- function(treeeOutput){
  print(treeeOutput$treee)
  invisible(treeeOutput)
}
