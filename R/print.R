#' @export
print.TreeeNode <- function(x, ...){
  cat(paste0("Node ",x$currentIndex,":\n"))
  cat("Number of observations: ", length(x$idxRow), "\n")
  # cat("Number of variables: ", length(x$idxCol), "\n")
  cat("Parent:", x$parent, "\n")
  cat("Children:", x$children, "\n")
  cat("Node model:", x$nodeModel, "\n")
  # cat("Current Loss:", x$currentLoss, "\n")
  cat("Within-node accuracy:", x$accuracy, "\n")
  invisible(x)
}


#' @export
summary.TreeeNode <- function(object, ...){
  # just print out everything besides some super long info
  object$idxRow <- object$idxCol <- object$nodePredict <- object$splitFun <- NULL
  return(unclass(object))
}
