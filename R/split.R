getSplitFunLDA <- function(datX, modelULDA){
  if(modelULDA$predGini <= 0.1) modelULDA$prior[] <- 1 / length(modelULDA$prior) # change to equal prior
  return(getSplitFunLDAhelper(datX = datX, modelULDA = modelULDA))
}


#' Helper Function for LDA-based Splitting in Tree Construction
#'
#' This function generates a splitting function based on a fitted ULDA model. It
#' assigns observations to the class with the minimal classification cost, and
#' returns the corresponding split results.
#'
#' @noRd
getSplitFunLDAhelper <- function(datX, modelULDA){
  predictedOutcome <- predict(modelULDA, datX)
  if(length(unique(predictedOutcome)) == 1)  return(NULL) # a valid split needs at least two child nodes

  idxPred <- which(names(modelULDA$prior) %in% predictedOutcome) # in case some classes are not predicted
  splitRes <- lapply(idxPred, function(i) which(names(modelULDA$prior)[i] == predictedOutcome))

  res <- function(datX){
    predictedProb <- predict(modelULDA, datX, type = "prob")[, idxPred, drop = FALSE]
    predictedOutcome <- max.col(-predictedProb %*% t(modelULDA$misClassCost[idxPred, idxPred, drop = FALSE]), ties.method = "first")
    return(lapply(seq_along(idxPred), function(i) which(i == predictedOutcome)))
  }

  attr(res, "splitRes") <- splitRes # record the split function's form
  attr(res, "model") <- modelULDA
  return(res)
}








