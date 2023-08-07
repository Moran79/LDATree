
# Check for input prior ---------------------------------------------------

checkPrior <- function(prior, response){
  ## Modified from randomForest.default Line 114
  if (is.null(prior)) {
    prior <- table(response) / length(response) # Default: Estimated Prior
  } else {
    if (length(prior) != nlevels(response))
      stop("length of prior not equal to number of classes")
    if (!is.null(names(prior))){
      prior <- prior[findTargetIndex(names(prior), levels(response))]
    }
    if (any(prior < 0)) stop("prior must be non-negative")
  }
  return(prior / sum(prior))
}


# Check for input misclassification cost ---------------------------------------------------

checkMisClassCost <- function(misClassCost, response){
  if (is.null(misClassCost)) {
    misClassCost <- (1 - diag(nlevels(response))) # Default: 1-identity
    colnames(misClassCost) <- rownames(misClassCost) <- levels(response)
  } else {
    if (dim(misClassCost)[1] != dim(misClassCost)[2] | dim(misClassCost)[1] != nlevels(response))
      stop("misclassification costs matrix has wrong dimension")
    if(!all.equal(colnames(misClassCost), rownames(misClassCost)))
      stop("misClassCost: colnames should be the same as rownames")
    if (!is.null(colnames(misClassCost))){
      misClassCost <- misClassCost[findTargetIndex(colnames(misClassCost), levels(response)),
                                   findTargetIndex(colnames(misClassCost), levels(response))]
    }
  }
  return(misClassCost)
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
  prior <- checkPrior(prior = prior, response = v)

  summary_table <- table(v) * prior
  if(posterior){return(summary_table / sum(summary_table))}
  return(names(which.max(summary_table)))
}


# Missing Value Imputation ------------------------------------------------

missingFix <- function(data, method){
  #> data: a data.frame
  #> Even the training data does not contain NA, this step is necessary since testing might have NA
  #> However, in this case there will be no <new missing level> added
  #>
  #> missingMethod: for numerical / categorical variables, respectively
  #> Constant check should not be done in this step,
  #> since it violates the one-function-at-a-time principle

  misMethod <- misMethodHelper(method = method)

  data <- droplevels(data) # empty levels due to partition
  data <- createFlagColumns(data = data, method = method) # create columns

  #> 关于missing flag会被当成numerical从而填出来0.67这种情况
  #> 我记得好像在getDataInShape那个函数中得到了解决，所以此处先乱填好了

  numOrNot <- getNumFlag(data) # num or cat
  NAcolumns <- sapply(data, anyNA)
  dataNRef <- rbind(data, NA) # add ref to the last row of data, copied from the first row

  for(i in seq_len(ncol(dataNRef))){
    if(numOrNot[i]){
      # numerical / logical vars
      # This function works even when all entries are NA, output NaN
      dataNRef[is.na(dataNRef[,i]),i] <- do.call(misMethod$numMethod, list(dataNRef[,i], na.rm = TRUE))
    }else{
      # categorical vars
      if(NAcolumns[i]){ # any NA
        dataNRef[,i] <- as.character(dataNRef[,i]) # for new level addition
        dataNRef[is.na(dataNRef[,i]),i] <- ifelse(misMethod$catMethod == "newLevel", "new0_0Level", as.character(getMode(dataNRef[,i])))
        dataNRef[,i] <- as.factor(dataNRef[,i]) # level information will be used in prediction
      }else{
        dataNRef[nrow(dataNRef),i] <- getMode(dataNRef[,i])
      }
    }
  }

  return(list(data = dataNRef[-nrow(dataNRef),,drop = FALSE],
              ref = dataNRef[nrow(dataNRef),,drop = FALSE]))
}

createFlagColumns <- function(data, method){
  #> given a data and numOrNot (from getNumFlag)
  #> output a data with added flag columns with correct 0/1
  #> We only add _FLAG to vars where NAs exist, not all columns
  #> since even we add NAflags, they will not be trained

  misMethod <- misMethodHelper(method = method)
  numOrNot <- getNumFlag(data) # num or cat
  NAcolumns <- sapply(data, anyNA)

  if(misMethod$numFlagOrNot & sum(numOrNot) > 0){
    NAcol <- which(numOrNot & NAcolumns)
    if(length(NAcol) > 0){
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

misMethodHelper <- function(method){
  # Methods preparation
  numMethod <- match.arg(method[1], c("mean", "median", "meanFlag", "medianFlag"))
  catMethod <- match.arg(method[2], c("mode", "modeFlag", "newLevel"))
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

constantColCheck <- function(data, idx = NULL, tol = 1e-8){
  if(is.null(idx)) {idx <- seq_len(ncol(data))}  # default output columns
  #> constant columns fix
  #> c(NA,1,1,1,NA) should be treated as non-constant

  constantColCheckHelper <- function(x, tol = 1e-8){
    #> one less digit than constantCol is helpful
    #> since the variance calculation in constantGroup would cause non-stopping error
    if(getNumFlag(x)) {x <- round(x, digits = -1-log(tol,base = 10))}
    return(length(unique(x)) > 1)
  }

  stopifnot(is.data.frame(data))

  idxNotConst <- which(sapply(data, function(x) constantColCheckHelper(x, tol)))

  idx <- idx[idxNotConst]

  return(idx)
}



# Stop check --------------------------------------------------------------

stopCheck <- function(responseCurrent, idxCol, maxTreeLevel, minNodeSize, currentLevel){
  # 0: Normal
  # 1: Stop and return posterior majority
  # 2: stop and fit LDA

  flagNodeSize <- length(responseCurrent) <= 2 * minNodeSize
  flagTreeLevel <- currentLevel >= maxTreeLevel # LDA is possible
  flagCol <- length(idxCol) == 0
  flagResponse <- length(unique(responseCurrent)) == 1

  if (flagNodeSize | flagResponse | flagCol) {return(1)}
  if (flagTreeLevel) {return(2)}
  return(0)
}



# Get LD scores -----------------------------------------------------------

getLDScores <- function(modelLDA, data, nScores = 1){
  #> data.frame: change to design matrix
  if (is.data.frame(data)) {data <- model.matrix(~., data)}
  #> one line of design matrix: change it to matrix so that colnames works
  if (is.vector(data)) {data <- matrix(data,1,dimnames = list(NULL, names(data)))}
  match(rownames(modelLDA$scaling), colnames(data))
  colIdx <- findTargetIndex(nameObj = colnames(data), nameTarget = rownames(modelLDA$scaling))
  LDScores <- data[,colIdx] %*% modelLDA$scaling
  return(LDScores[,seq_len(nScores)])
}


# Gini Split --------------------------------------------------------------

GiniSplitScreening <- function(xCurrent, response, idxRow, prior, minNodeSize, misClassCost, modelLDA){
  Nj = table(response) # overall proportion, for prior calculation
  response <- response[idxRow]
  posterior = getMode(v = response, prior = prior / Nj, posterior = TRUE) # p(j,t) = prior * Njt / Nj

  LDScore <- getLDScores(modelLDA = modelLDA, data = xCurrent, nScores = 1)
  idxRowOrdered <- order(LDScore)

  #> prevent empty nodes, so the lowest rank is removed
  #> max cut: 1000
  #> percentage cut: to satisfy the minNode Constraints
  #> potentialCut: LDScores' ranks, a subset of 1 to length(response)
  percentageCut <- minNodeSize / length(response) # there is a previous stopCheck that size >= 2 * minNodeSize
  stopifnot(percentageCut <= 0.5)
  potentialCut <- unique(quantile(rank(LDScore,ties.method = "min"),
                                  probs = seq(percentageCut, 1 - percentageCut,length.out = 1000), type = 1))
  potentialCut <- setdiff(potentialCut,1)

  if(length(potentialCut)==0) {return(NULL)} # No cut due to ties

  # For the consideration of speed
  # increment programming is carried out
  GiniObserved <- numeric(length(potentialCut))
  NjtLeft <- table(response[idxRowOrdered][seq_len(potentialCut[1] - 1)])
  NjtRight <- table(response[idxRowOrdered][seq(potentialCut[1], length(response))])
  posteriorLeft <- prior * NjtLeft / Nj / sum(prior * NjtLeft / Nj)
  posteriorRight <- prior * NjtRight / Nj / sum(prior * NjtRight / Nj)
  GiniObserved[1] <- sum(NjtLeft) * sum(misClassCost * posteriorLeft %*% t(posteriorLeft)) +
    sum(NjtRight) * sum(misClassCost * posteriorRight %*% t(posteriorRight))
  for(i in seq_along(potentialCut)[-1]){
    responseTrans <- table(response[seq(potentialCut[i-1],potentialCut[i]-1)])
    NjtLeft <- NjtLeft + responseTrans
    NjtRight <- NjtRight - responseTrans
    posteriorLeft <- prior * NjtLeft / Nj / sum(prior * NjtLeft / Nj)
    posteriorRight <- prior * NjtRight / Nj / sum(prior * NjtRight / Nj)
    GiniObserved[i] <- sum(NjtLeft) * sum(misClassCost * posteriorLeft %*% t(posteriorLeft)) +
      sum(NjtRight) * sum(misClassCost * posteriorRight %*% t(posteriorRight))
  }
  cutPoint <- which.min(GiniObserved)
  return(list(left = idxRow[idxRowOrdered][seq_len(potentialCut[cutPoint] - 1)],
              right = idxRow[idxRowOrdered][seq(potentialCut[cutPoint], length(response))],
              cut = LDScore[idxRowOrdered][potentialCut[cutPoint]-1]))
}


# New level fix + Missing -----------------------------------------------------------

getDataInShape <- function(data, missingReference){
  #> change the shape of test data to the training data
  #> and make sure that the dimension of the data is the same as missingRefernce

  nameVarIdx <- match(colnames(missingReference), colnames(data))
  if(anyNA(nameVarIdx)){
    #> New columns fix (or Flags): If there are less columns than it should be,
    #> add columns with NA
    data[,colnames(missingReference)[which(is.na(nameVarIdx))]] <- NA
    nameVarIdx <- match(colnames(missingReference), colnames(data))
  }
  data <- data[,nameVarIdx, drop = FALSE]

  #> New levels fix
  data <- setLevelWithReference(data = data, reference = missingReference,
                                keepNA = TRUE)

  #> Missing Value fix
  #> Assumption: missingMethod is same for both training and test
  for(i in seq_len(ncol(data))){
    missingIndicator <- is.na(data[,i])

    #> If flag is present, find its index & impute NA flags if flag has NAs
    missingFlagIdx <- which(colnames(missingReference) == paste(colnames(data)[i],"FLAG",sep = "_"))
    if(length(missingFlagIdx) == 1 & anyNA(data[, missingFlagIdx])){
      NA_InFlagIdx <- which(is.na(data[, missingFlagIdx]))
      data[NA_InFlagIdx, missingFlagIdx] <- missingIndicator[NA_InFlagIdx] + 0
    }

    if(any(missingIndicator)){
      # separate methods for cat / num
      if(getNumFlag(data[,i])){
        data[which(missingIndicator), i] <- missingReference[1,i]
      }else{
        data[,i] <- as.character(data[,i])
        data[which(missingIndicator), i] <- as.character(missingReference[1,i])
        data[,i] <- factor(data[,i], levels = levels(missingReference[,i]))
      }
    }
  }

  return(data)
}


# Prediction in terminal Nodes --------------------------------------------

predNode <- function(data, treeeNode){
  #> data is a data.frame
  if(treeeNode$nodeModel == "LDA"){
    return(predict(treeeNode$nodePredict, data))
  }else{
    return(rep(treeeNode$nodePredict, dim(data)[1]))
  }
}


# Level solver with reference -----------------------------------------

setLevelWithReference <- function(data, reference, keepNA = FALSE){
  #> reference: a data.frame 1*p
  #> data: a data.frame
  #> assumption: data and reference has same colnames

  stopifnot(all(colnames(data) == colnames(reference)))

  #> each column is a number, or a character with levels
  #> 1. if keepNA = FALSE, all NAs / new levels in categorical vars
  #> will be replaced by the reference, but numerical variables stay untouched
  #> 2. if keepNA = TRUE, the new level will be changed to NA,
  #> and numerical variables stay untouched
  levelReference <- sapply(reference, levels)
  catIdx <- which(!sapply(levelReference,is.null))

  for(i in catIdx){
    # sapply would lose factor property but left character, why?
    data[,i] <- factor(data[,i], levels = levelReference[[i]])
    if(anyNA(data[,i]) & !keepNA){
      data[is.na(data[,i]),i] <- as.character(reference[1,i])
    }
  }
  return(data)
}

# Rewrite the LDA ---------------------------------------------------------

ldaGSVD <- function(formula, data, prior){
  # response <- as.factor(data[[all.vars(formula)[1]]])
  modelFrame <- model.frame(formula, data, na.action = "na.fail")
  Terms <- terms(modelFrame)
  response <- droplevels(as.factor(modelFrame[,1])) # some levels are branched out
  # if(missing(prior)) {prior <- table(response) / length(response)} # not used for now
  prior <- table(response) / length(response) # rewrite the prior = abandon prior
  m <- model.matrix(formula, data)
  cnames <- colnames(m)

  # Step 1: SVD on the combined matrix H
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  grandMeansJ <- matrix(colSums(m) / nrow(m), nrow = nlevels(response), ncol = ncol(m), byrow = TRUE)
  Hb <- sqrt(tabulate(response)) * (groupMeans - grandMeansJ)
  fitSVD <- svd(rbind(Hb, m - groupMeans[response,]))
  rankT <- sum(fitSVD$d >= max(dim(fitSVD$u),dim(fitSVD$v)) * .Machine$double.eps * fitSVD$d[1])

  # Step 2: SVD on the P matrix
  fitSVDp <- svd(fitSVD$u[seq_len(nlevels(response)), seq_len(rankT), drop = FALSE], nu = 0L)
  rankAll <- min(nlevels(response)-1, rankT) # This is not optimal, but rank(Hb) takes time
  # Fix the variance part
  unitSD <- pmin(diag(sqrt((length(response) - nlevels(response)) / abs(1 - fitSVDp$d^2)), nrow = rankAll),1e15) # Scale to unit var
  scalingFinal <- (fitSVD$v[,seq_len(rankT), drop = FALSE] %*% diag(1 / fitSVD$d[seq_len(rankT)], nrow = rankT) %*% fitSVDp$v)[,seq_len(rankAll), drop = FALSE] %*% unitSD
  rownames(scalingFinal) <- cnames

  groupMeans <- groupMeans %*% scalingFinal
  rownames(groupMeans) <- levels(response)

  res <- list(scaling = scalingFinal, formula = formula, terms = Terms, prior = prior,
              groupMeans = groupMeans, xlevels = .getXlevels(Terms, modelFrame))
  class(res) <- "ldaGSVD"
  return(res)
}

predict.ldaGSVD <- function(object, newdata){
  # add one extra check for levels of the predictors
  # browser()
  Terms <- delete.response(object$terms)
  modelX <- model.matrix(Terms, data = newdata, xlev = object$xlevels)
  LDscores <- modelX %*% object$scaling
  # browser()
  loglikelihood <- LDscores %*% t(object$groupMeans) + matrix(log(object$prior) - 0.5 * rowSums(object$groupMeans^2), nrow(LDscores), length(object$prior), byrow = TRUE)
  # Computation Optimization 2: Prevent a very large likelihood due to exponential
  likelihood <- exp(loglikelihood - apply(loglikelihood, 1, max))
  posterior <- likelihood / apply(likelihood, 1, sum)
  return(rownames(object$groupMeans)[max.col(posterior)])
}












