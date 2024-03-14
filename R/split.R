
# Main function -----------------------------------------------------------


getSplitFun <- function(datX, response, method, modelLDA){
  if(method == "LDA"){
    return(getSplitFunLDA(datX = datX,
                          response = response,
                          modelLDA = modelLDA))
  } else if(method != "LDA") stop("Please wait for new features to come")
}

# mixed -------------------------------------------------------------------


getSplitFunLDA <- function(datX, response, modelLDA){
  if(modelLDA$predGini <= 0.1) modelLDA$prior[] <- 1 / length(modelLDA$prior) # change to equal prior
  return(getSplitFunLDAhelper(datX = datX,
                              response = response,
                              modelLDA = modelLDA))
}

# LDA --------------------------------------------------------------------


getSplitFunLDAhelper <- function(datX, response, modelLDA){
  #> This function is called only when building the tree

  #> If there are some classes not being predicted
  #> we will assign them to the class with the second largest posterior prob
  predictedOutcome <- predict(modelLDA, datX)
  # if(length(unique(predictedOutcome)) == 1)  return(NULL) # This will never happens, delete before next release
  idxPred <- which(names(modelLDA$prior) %in% predictedOutcome)
  splitRes <- lapply(idxPred, function(i) which(names(modelLDA$prior)[i] == predictedOutcome))

  res <- function(datX, missingReference){
    # fixedData <- getDataInShape(data = datX, missingReference = missingReference, NBmethod = modelLDA$NBmethod)
    fixedData <- getDataInShape(data = datX, missingReference = missingReference)
    # browser()
    predictedProb <- predict(modelLDA, fixedData,type = "prob")[,idxPred, drop = FALSE]
    predictedOutcome <- max.col(predictedProb, ties.method = "first")
    return(lapply(seq_along(idxPred), function(i) which(i == predictedOutcome)))
  }

  attr(res, "splitRes") <- splitRes # record the split function's form
  return(res)
}








