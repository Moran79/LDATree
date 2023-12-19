
# Missing Value Imputation ------------------------------------------------


missingFix <- function(data, missingMethod = c("meanFlag", "newLevel")){

  #> data: a data.frame
  #> missingMethod: for numerical / categorical variables, respectively

  misMethod <- misMethodHelper(missingMethod = missingMethod)

  data <- createFlagColumns(data = data, misMethod = misMethod) # create flag columns

  numOrNot <- getNumFlag(data) # num or cat
  NAcolumns <- sapply(data, anyNA)
  dataNRef <- rbind(data, NA) # add ref to the last row of data, initialize using NAs

  for(i in seq_len(ncol(dataNRef))){
    if(numOrNot[i]){
      # numerical / logical vars
      #> The function below output NaN when all entries are NA,
      #> so we add an if to prevent a vector of only NA and NaN (not constant)
      if(!all(is.na(dataNRef[,i]))) dataNRef[is.na(dataNRef[,i]),i] <- do.call(misMethod$numMethod, list(dataNRef[,i], na.rm = TRUE))
    }else{
      # categorical vars
      if(NAcolumns[i]){ # any NA
        dataNRef[,i] <- as.character(dataNRef[,i]) # for new level addition
        dataNRef[is.na(dataNRef[,i]),i] <- ifelse(misMethod$catMethod == "newLevel", "new0_0Level", as.character(getMode(dataNRef[,i])))
        dataNRef[,i] <- as.factor(dataNRef[,i]) # level information will be used in prediction
      }else{
        dataNRef[,i] <- as.factor(dataNRef[,i])
        dataNRef[nrow(dataNRef),i] <- factor(getMode(dataNRef[,i]), levels = levels(dataNRef[,i]))
      }
    }
  }

  return(list(data = dataNRef[-nrow(dataNRef),,drop = FALSE],
              ref = dataNRef[nrow(dataNRef),,drop = FALSE]))
}

createFlagColumns <- function(data, misMethod){
  #> given a data and numOrNot (from getNumFlag)
  #> output a data with added flag columns with correct 0/1
  #> We only add _FLAG to vars where NAs exist, not all columns
  #> since even we add NAflags, they will not be trained
  #> Notes: num flags are before cat flags,
  #> regardless of their relative column positions

  numOrNot <- getNumFlag(data) # num or cat
  NAcolumns <- sapply(data, anyNA)

  if(misMethod$numFlagOrNot & sum(numOrNot) > 0){
    NAcol <- which(numOrNot & NAcolumns)
    if(length(NAcol) > 0){
      #> The line below will treat missing flags as numerical variables
      #> as.factor() can be applied if we want them to be factors
      numFlagCols <- do.call(cbind, sapply(NAcol, function(colIdx) is.na(data[, colIdx])+0, simplify = FALSE))
      colnames(numFlagCols) <- paste(colnames(data)[NAcol], "FLAG", sep = "_")
      data <- cbind(data, numFlagCols)
    }
  }

  if(misMethod$catFlagOrNot & sum(!numOrNot) > 0){
    NAcol <- which((!numOrNot) & NAcolumns)
    if(length(NAcol) > 0){
      catFlagCols <- do.call(cbind, sapply(NAcol, function(colIdx) is.na(data[, colIdx])+0, simplify = FALSE))
      colnames(catFlagCols) <- paste(colnames(data)[NAcol], "FLAG", sep = "_")
      data <- cbind(data, catFlagCols)
    }
  }
  return(data)
}

misMethodHelper <- function(missingMethod){
  #> Aim: classify the missing imputation methods based on their methods and flags
  numMethod <- missingMethod[1]; catMethod <- missingMethod[2]
  numFlagOrNot <- grepl("Flag", numMethod)
  catFlagOrNot <- grepl("Flag", catMethod)
  numMethod <- ifelse(grepl("mean", numMethod), "mean", "median")
  catMethod <- ifelse(grepl("mode", catMethod), "mode", "newLevel")
  return(list(numMethod = numMethod,
              catMethod = catMethod,
              numFlagOrNot = numFlagOrNot,
              catFlagOrNot = catFlagOrNot))
}


# constant check ---------------------------------------------------------

constantColCheck <- function(data, idx, tol = 1e-8){
  if(missing(idx)) idx <- seq_len(ncol(data))  # default output columns
  #> constant columns fix: the data in this step should not contains NA

  constantColCheckHelper <- function(x, tol = 1e-8){
    if(getNumFlag(x)) x <- round(x, digits = -log(tol,base = 10))
    return(length(unique(x)) > 1)
  }

  idxNotConst <- which(sapply(data, function(x) constantColCheckHelper(x, tol)))

  return(idx[idxNotConst])
}

# Input check helper ------------------------------------------------------

findTargetIndex <- function(nameObj, nameTarget){
  #> Assume nameObj and nameTarget are of the same length
  #> check if nameObj are all in nameTarget
  #> If yes, return the corresponding index
  #> so that nameObj[idx] == nameTarget
  targetIndex <- match(nameTarget, nameObj)
  if (anyNA(targetIndex)) {
    stop("The names do not match with the response")
  }
  return(targetIndex)
}


# Get column Types --------------------------------------------------------

getNumFlag <- function(x, index = FALSE){
  #> index decides return type. e.g. return c(1,0,0,1,1,0) (FALSE) or c(1,4,5) (YES)
  #> logical has to be included, a column with all NAs has to be viewed as numeric
  if(is.null(dim(x))){return(class(x) %in% c('numeric', 'integer', 'logical'))}

  stopifnot(is.data.frame(x))

  numOrNot <- sapply(x, class) %in% c('numeric', 'integer', 'logical')

  if(index){numOrNot <- which(numOrNot)}

  return(numOrNot)
}

# Get mode --------------------------------------------------------------------

getMode <- function(v, prior, posterior = FALSE){
  #> posterior: return posterior probs (TRUE) or mode (FALSE)
  #> NA will be ignored
  v <- as.factor(v)
  if(missing(prior)) prior = rep(1,nlevels(v)) # equal prior

  summary_table <- table(v) * prior
  if(posterior){return(summary_table / sum(summary_table))}
  return(names(which.max(summary_table)))
}



# Stop check --------------------------------------------------------------

stopCheck <- function(responseCurrent, numCol, maxTreeLevel, minNodeSize, currentLevel){
  # 0: Normal
  # 1: Stop and return posterior majority
  # 2: stop and fit LDA

  flagNodeSize <- length(responseCurrent) <= minNodeSize # 数据量不够了，LDA is unstable
  flagTreeLevel <- currentLevel >= maxTreeLevel # 层数到了
  flagCol <- numCol == 0 # no covs left
  flagResponse <- length(unique(responseCurrent)) == 1 # 只有一种y

  if (flagResponse | flagCol | flagNodeSize) {return(1)}
  if (flagTreeLevel) {return(2)}
  return(0)
}



# Get LD scores -----------------------------------------------------------

getDesignMatrix <- function(modelLDA, data, scale = FALSE){
  # Output: the design matrix
  Terms <- delete.response(modelLDA$terms)
  modelX <- model.matrix(Terms, data = data, xlev = modelLDA$xlevels)

  if(scale){ # reserved for the scaling in getLDscores
    modelX <- sweep(modelX[,modelLDA$varIdx,drop = FALSE], 2, modelLDA$varCenter, "-")
    modelX <- sweep(modelX, 2, modelLDA$varSD, "/")
  }
  return(modelX)
}

getLDscores <- function(modelLDA, data, nScores = -1){

  modelX <- getDesignMatrix(modelLDA = modelLDA, data = data, scale = TRUE)
  if(nScores > 0) modelLDA$scaling <- modelLDA$scaling[, seq_len(nScores), drop = FALSE]
  LDscores <- modelX %*% modelLDA$scaling

  return(LDscores)
}


# New level fix + Missing -----------------------------------------------------------

getDataInShape <- function(data, missingReference){
  #> change the shape of test data to the training data
  #> and make sure that the dimension of the data is the same as missingRefernce

  cname <- colnames(missingReference)
  nameVarIdx <- match(cname, colnames(data))
  if(anyNA(nameVarIdx)){
    #> New columns fix (or Flags): If there are less columns than it should be,
    #> add columns with NA
    data[,cname[which(is.na(nameVarIdx))]] <- NA
    nameVarIdx <- match(cname, colnames(data))
  }

  for(currentIdx in seq_len(ncol(missingReference))){ # Main program starts
    #> The tricky part is the iterator is based on the missingReference, NOT the data.
    newIdx <- nameVarIdx[currentIdx]
    numOrNot <- getNumFlag(missingReference[, currentIdx])

    ### New-level Fix for Categorical Variable ###
    if(!numOrNot) data[, newIdx] <- factor(data[, newIdx], levels = levels(missingReference[, currentIdx])) # may generate NAs

    missingOrNot <- is.na(data[, newIdx])
    if(!any(missingOrNot)) next

    ### Flag Variable Detection ###
    #> This part is actually more complicated than expected
    #> Four combinations could happen: Ori / Flag both can have NA or complete
    currentVarName <- cname[currentIdx]

    ## Scenario 1: It has a related flag variable in the data ##
    #> Only modify those flags where the original variable is missing
    #> Keep other parts still, since there could already be imputed values
    #> in the original variable that have been taken care of
    currentFlagIdx <- which(cname == paste(currentVarName,"FLAG",sep = "_"))
    if(length(currentFlagIdx) == 1) data[which(missingOrNot), nameVarIdx[currentFlagIdx]] <- 1

    ## Scenario 2: It is a flag and it has an original variable in (or not in) the data ##
    #> Only impute those NAs in the flag, but keep the values that are already in the flag
    if(grepl("_FLAG$", currentVarName)){
      orginalVarName <- sub("_FLAG$", "", currentVarName)
      orginalVarIdx <- which(cname == orginalVarName)
      if(length(orginalVarIdx) == 1) data[which(missingOrNot), newIdx] <- is.na(data[which(missingOrNot), nameVarIdx[orginalVarIdx]]) + 0
      else data[, newIdx] <- 1 # The original data is NOT found
      next
    }

    ### For numerical & categorical variables ###
    data[which(missingOrNot), newIdx] <- missingReference[1, currentIdx]

    # data[, newIdx] <- as.character(data[, newIdx])
    # data[which(missingOrNot), newIdx] <- as.character(missingReference[1, currentIdx])
    # data[, newIdx] <- factor(data[, newIdx], levels = levels(missingReference[, currentIdx]))

  }

  return(data[,nameVarIdx, drop = FALSE])
}




# Prediction in terminal Nodes --------------------------------------------

predNode <- function(data, treeeNode, type){
  #> data is a data.frame
  if(treeeNode$nodeModel == "LDA"){
    return(predict(object = treeeNode$nodePredict, newdata = data, type = type))
  }else{
    if(type == "response") return(rep(treeeNode$nodePredict, dim(data)[1]))
    else{ # if type = "all", the extra response column will be added later
      pred <- matrix(0,nrow = nrow(data), ncol = length(treeeNode$proportions), dimnames = list(c(), names(treeeNode$proportions)))
      pred[,which(treeeNode$nodePredict == colnames(pred))] <- 1
      return(pred)
    }
  }
}



# Get the p-value for testing the current nodes' performance --------------

getOneSidedPvalue <- function(N, xBefore, xAfter){
  #> H1: xBefore > xAfter. x stands for the prediction error
  zStat <- (xBefore - xAfter) / sqrt((xBefore * (N - xBefore) + xAfter * (N - xAfter)) / N + 1e-16)
  pnorm(zStat, lower.tail = FALSE)
}



# Pruning -----------------------------------------------------------------

#> In the future, I might have to change the name of alpha to\
#> something like <drop in training error>

pruneByTrainErrDrop <- function(treeeList, pThreshold = 0.01, verbose = TRUE){
  stopifnot(class(treeeList) == "SingleTreee")
  lenBefore <- length(treeeList)
  treeeList <- makeAlphaMono(treeeList = treeeList)
  treeeList <- cutNodes(treeeList = treeeList, pThreshold = pThreshold)
  treeeList <- dropNodes(treeeList = treeeList)
  lenAfter <- length(treeeList)
  if(verbose) cat(paste(lenBefore - lenAfter, 'nodes are being pruned.\n'))
  return(treeeList)
}

makeAlphaMono <- function(treeeList){
  #> Purpose: make the alpha monotonic

  for(i in rev(seq_along(treeeList))){
    if(is.null(treeeList[[i]]$pruned)){ # Only loop over the unpruned nodes
      if(is.null(treeeList[[i]]$children)) treeeList[[i]]$alpha <- 1 # for terminal nodes
      else{ # for intermediate nodes
        childrenAlpha <- sapply(treeeList[[i]]$children, function(idx) treeeList[[idx]]$alpha)
        treeeList[[i]]$alpha <- min(c(1, treeeList[[i]]$alpha, unlist(childrenAlpha)), na.rm = TRUE)
      }
    }
  }
  return(treeeList)
}

getTerminalNodes <- function(currentIdx, treeeList, keepNonTerminal = FALSE){
  #> Get all terminal nodes that are offsprings from currentIdx
  treeeNode <- treeeList[[currentIdx]]
  if(is.null(treeeNode$children)){
    return(currentIdx)
  }else{
    nonTerminalNodes <- if(keepNonTerminal) currentIdx
    terminalNodes <- do.call(c, sapply(treeeNode$children,function(x) getTerminalNodes(currentIdx = x, treeeList = treeeList, keepNonTerminal = keepNonTerminal), simplify = FALSE))
    return(c(nonTerminalNodes, terminalNodes))
  }
}

cutNodes <- function(treeeList, pThreshold){
  for(i in rev(seq_along(treeeList))){
    treeeNode <- treeeList[[i]]
    # not yet pruned + non-terminal node + p-value above threshold
    cutFlag <- is.null(treeeNode$pruned) & !is.null(treeeNode$children) & treeeNode$alpha > pThreshold # R rounding error, 1e-10 needed
    if(cutFlag){
      allChildren <- getTerminalNodes(currentIdx = i, treeeList = treeeList, keepNonTerminal = TRUE)
      for(j in setdiff(allChildren, i)) {treeeList[[j]]$pruned <- TRUE}
      treeeList[[i]]["children"] <- list(NULL) # cut the branch
    }
  }
  # treeeList <- makeAlphaMono(treeeList = treeeList) # no need to update alpha, since it is not changed
  return(treeeList)
}


dropNodes <- function(treeeList){
  finalNodeIdx <- sort(getTerminalNodes(treeeList = treeeList, currentIdx = 1, keepNonTerminal = TRUE))
  treeeList <- treeeList[finalNodeIdx]
  class(treeeList) <- "SingleTreee"
  for(i in seq_along(treeeList)){
    treeeList[[i]]$currentIndex <- i # re-assign the currentIndex
    if(!is.null(treeeList[[i]]$children)) {
      treeeList[[i]]$children <- sapply(treeeList[[i]]$children, function(x) which(finalNodeIdx == x))
    }else{
      if(treeeList[[i]]$stopFlag == 0) treeeList[[i]]$stopFlag = 3 # due to pruning
    }
    treeeList[[i]]$parent <- which(finalNodeIdx == treeeList[[i]]$parent)
  }
  return(treeeList)
}








