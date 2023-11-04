
# Main function -----------------------------------------------------------

getSplitFun <- function(x, response, method, modelLDA){
  if(method == "LDscores") return(getSplitFunLDscores(x = x,
                                                      response = response,
                                                      modelLDA = modelLDA))
  else if(method == "FACT") return(getSplitFunFACT(x = x,
                                                   response = response,
                                                   modelLDA = modelLDA))
  else if(method == "Pillai") return(getSplitFunPillai(x = x,
                                                   response = response,
                                                   modelLDA = modelLDA))
}


# Gini Split --------------------------------------------------------------

getSplitFunLDscores <- function(x, response, modelLDA){

  LDscore <- getLDscores(modelLDA = modelLDA, data = x, nScores = 1)
  idxRowOrdered <- order(LDscore)

  #> prevent empty nodes, so the lowest rank is removed
  #> max cut: 1000
  #> potentialCut: LDscores' ranks, a subset of 1 to length(response)
  percentageCut <- 0.05
  potentialCut <- unique(quantile(rank(LDscore,ties.method = "min"),
                                  probs = seq(percentageCut, 1 - percentageCut,length.out = 1000), type = 1))
  potentialCut <- setdiff(potentialCut,1)

  if(length(potentialCut)==0) {return(NULL)} # No cut due to ties

  # For the consideration of speed
  # increment programming is carried out
  GiniObserved <- numeric(length(potentialCut))
  NjtLeft <- table(response[idxRowOrdered][seq_len(potentialCut[1] - 1)])
  NjtRight <- table(response[idxRowOrdered][seq(potentialCut[1], length(response))])
  GiniObserved[1] <- sum(NjtLeft) * sum((NjtLeft / sum(NjtLeft))^2) +
    sum(NjtRight) * sum((NjtRight / sum(NjtRight))^2)
  for(i in seq_along(potentialCut)[-1]){
    responseTrans <- table(response[seq(potentialCut[i-1],potentialCut[i]-1)])
    NjtLeft <- NjtLeft + responseTrans
    NjtRight <- NjtRight - responseTrans
    GiniObserved[i] <- sum(NjtLeft) * sum((NjtLeft / sum(NjtLeft))^2) +
      sum(NjtRight) * sum((NjtRight / sum(NjtRight))^2)
  }
  cutPoint <- which.max(GiniObserved)
  cutScore <- LDscore[idxRowOrdered][[potentialCut[cutPoint]-1]]

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
    LDscore <- unname(getLDscores(modelLDA = modelLDA, data = fixedData, nScores = 1))
    # return the relative index
    return(list(which(LDscore<=cutScore), which(LDscore>cutScore)))
  }
  return(res)
}

# FACT --------------------------------------------------------------------

getSplitFunFACT <- function(x, response, modelLDA){
  #> This function is called only when building the tree
  predictedOutcome <- predict(modelLDA, x)
  if(length(unique(predictedOutcome)) == 1) return(NULL)

  #> If there are some classes not being predicted
  #> we will assign them to the class with the second largest posterior prob
  idxPred <- which(names(modelLDA$prior) %in% predictedOutcome)

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
    predictedProb <- predict(modelLDA, fixedData,type = "prob")[,idxPred, drop = FALSE]
    predictedOutcome <- max.col(predictedProb, ties.method = "first")
    lapply(seq_along(idxPred), function(i) which(i == predictedOutcome))
    return(lapply(seq_along(idxPred), function(i) which(i == predictedOutcome)))
  }
}

# getSplitFunFACT <- function(x, response, modelLDA){
#   #> make sure that predictions have all the classes
#   #> otherwise, refit the modelLDA
#   predictedOutcome <- predict(modelLDA, x)
#   while(length(unique(predictedOutcome)) != length(modelLDA$prior)){
#     if(length(unique(predictedOutcome)) == 1) return(NULL) # if only one class is left
#     subsetIdx <- which(response %in% unique(predictedOutcome))
#     x <- x[subsetIdx,, drop = FALSE]; response <- response[subsetIdx]
#     datCombined = data.frame(response = response, x)
#     modelLDA <- ldaGSVD(response~., data = datCombined, method = "step")
#     predictedOutcome <- predict(modelLDA, x)
#   }
#
#   res <- function(x, missingReference){
#     fixedData <- getDataInShape(data = x, missingReference = missingReference)
#     predictedOutcome <- predict(modelLDA, fixedData)
#     return(lapply(names(modelLDA$prior), function(name) which(name == predictedOutcome)))
#   }
# }


# Sum of Pillai's trace ---------------------------------------------------

getSplitFunPillai <- function(x = x, response = response, modelLDA = modelLDA){
  x <- getDesignMatrix(modelLDA = modelLDA, data = x) # get the scaled x
  M <- tapply(c(x), list(rep(response, dim(x)[2]), col(x)), function(o_o) mean(o_o, na.rm = TRUE))
  n <- nrow(x); nJ <- as.vector(table(response))

  getSumPillai <- function(a){ # a is the vector of interest
    # a[1] is the constant, a[-1] are coefficients of the variables
    group1Flag <- as.vector(x %*% a[-1] >= a[1])
    n1 <- sum(group1Flag)
    if(n1 %in% c(0, n)) return(Inf) # If the split is outside the data range

    # Calculate the smaller part, and use subtraction to get the other
    if(n1 > n / 2){
      group1Flag <- !group1Flag
      n1 <- sum(group1Flag)
    }
    n2 <- n - n1; idx <- which(group1Flag)
    nJ1 <- as.vector(table(response[group1Flag]))

    M1 <- tapply(c(x[idx, , drop = FALSE]), list(rep(response[idx], dim(x[idx, , drop = FALSE])[2]), col(x[idx, , drop = FALSE])), function(o_o) mean(o_o, na.rm = TRUE))
    M2 <- (M * nJ - M1 * nJ1) / (nJ - nJ1)
    Xmean1 <- apply(x[idx, , drop = FALSE],2,mean)
    Xmean2 <- - n1 * Xmean1 / n2
    Sb1 <- t(t(M1) - Xmean1)
    Sb2 <- t(t(M2) - Xmean2)
    -(n1 * sum(Sb1^2) + n2 * sum(Sb2^2))
  }
  aScaled <- optim(c(0, rep(1, ncol(x))), getSumPillai, control = list(maxit = 100), method = "SANN")$par
  # # scale center transformation
  aFinal <- aScaled
  aFinal[-1] <- aFinal[-1] / modelLDA$varSD
  aFinal[1] <- aFinal[1] + sum(aFinal[-1] * modelLDA$varCenter)

  # Stop the split if all points belong to one side
  projectionOnSplit <- unname(as.vector(x %*% matrix(aScaled[-1], ncol = 1)))
  currentList <- list(which(projectionOnSplit < aScaled[1]), which(projectionOnSplit >= aScaled[1]))
  if(any(sapply(currentList, length) == 0)) return(NULL)

  res <- function(x, missingReference){
    fixedData <- getDataInShape(data = x, missingReference = missingReference)
    x <- getDesignMatrix(modelLDA = modelLDA, data = fixedData)
    projectionOnSplit <- unname(as.vector(x %*% matrix(aScaled[-1], ncol = 1)))
    # return the relative index
    return(list(which(projectionOnSplit < aScaled[1]), which(projectionOnSplit >= aScaled[1])))
  }
  attr(res, "a") <- aFinal # record the split function's form
  return(res)
}


