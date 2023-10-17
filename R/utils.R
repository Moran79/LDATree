
# Get x and response from formula and data --------------------------------

extractXnResponse <- function(formula, data){
  #> droplevels is necessary, since empty response level occurs during train/test split
  #> covariates can have empty levels as well
  modelFrame <- droplevels(model.frame(formula, data, na.action = "na.pass"))
  stopifnot(!anyNA(modelFrame[,1])) # no NAs are allowed in the response variable

  response <- as.factor(modelFrame[,1])
  x <- modelFrame[,-1, drop = FALSE]
  return(list(x = x, response = response))
}

# Missing Value Imputation ------------------------------------------------

missingFix <- function(data, missingMethod){

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
  numMethod <- match.arg(missingMethod[1], c("mean", "median", "meanFlag", "medianFlag"))
  catMethod <- match.arg(missingMethod[2], c("mode", "modeFlag", "newLevel"))
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

stopCheck <- function(responseCurrent, idxCol, maxTreeLevel, minNodeSize, currentLevel){
  # 0: Normal
  # 1: Stop and return posterior majority
  # 2: stop and fit LDA

  flagNodeSize <- length(responseCurrent) <= minNodeSize # 数据量不够了，LDA is unstable
  flagTreeLevel <- currentLevel >= maxTreeLevel # 层数到了
  flagCol <- length(idxCol) == 0 # no covs left
  flagResponse <- length(unique(responseCurrent)) == 1 # 只有一种y

  if (flagResponse | flagCol | flagNodeSize) {return(1)}
  if (flagTreeLevel) {return(2)}
  return(0)
}



# Get LD scores -----------------------------------------------------------

getLDscores <- function(modelLDA, data, nScores = -1){

  Terms <- delete.response(modelLDA$terms)
  modelX <- model.matrix(Terms, data = data, xlev = modelLDA$xlevels)
  modelX <- sweep(modelX[,modelLDA$varIdx,drop = FALSE], 2, modelLDA$varCenter, "-")
  modelX <- sweep(modelX, 2, modelLDA$varSD, "/")

  if(nScores > 0) modelLDA$scaling <- modelLDA$scaling[, seq_len(nScores), drop = FALSE]
  LDscores <- modelX %*% modelLDA$scaling

  return(LDscores)
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
  #> The columns are the same, now fix missing values and new levels


  #> New levels fix
  levelReference <- sapply(missingReference, levels, simplify = FALSE)
  for(i in which(!sapply(levelReference,is.null))){
    # sapply would lose factor property but left character, why?
    data[,i] <- factor(data[,i], levels = levelReference[[i]])
  }


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


fixNewLevel <- function(datTest, datTrain){
  #> change the shape of test data to the training data
  #> and make sure that the dimension of the data is the same as missingRefernce

  nameVarIdx <- match(colnames(datTrain), colnames(datTest))
  if(anyNA(nameVarIdx)){
    #> New columns fix (or Flags): If there are less columns than it should be,
    #> add columns with NA
    datTest[,colnames(datTrain)[which(is.na(nameVarIdx))]] <- NA
    nameVarIdx <- match(colnames(datTrain), colnames(datTest))
  }
  datTest <- datTest[,nameVarIdx, drop = FALSE]
  #> The columns are the same, now fix new levels

  #> change all characters to factors
  idxC <- which(sapply(datTrain, class) == "character")
  if(length(idxC) != 0){
    for(i in idxC){
      datTrain[,i] <- as.factor(datTrain[,i])
    }
  }

  #> New levels fix
  levelReference <- sapply(datTrain, levels, simplify = FALSE)
  for(i in which(!sapply(levelReference,is.null))){
    # sapply would lose factor property but left character, why?
    datTest[,i] <- factor(datTest[,i], levels = levelReference[[i]])
  }

  return(datTest)
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


# update Current Loss using Validation Set --------------------------------

updateCurrentLoss <- function(){

}














