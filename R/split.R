
# Main function -----------------------------------------------------------

getSplitFun <- function(x, response, method, modelLDA){
  if(method == "LDscores") return(getSplitFunLDscores(x = x,
                                                      response = response,
                                                      modelLDA = modelLDA))
  else if(method == "FACT") return(getSplitFunFACT(x = x,
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
    LDscore <- getLDscores(modelLDA = modelLDA, data = fixedData, nScores = 1)
    # return the relative index
    return(list(which(LDscore<=cutScore), which(LDscore>cutScore)))
  }
  return(res)
}

# FACT --------------------------------------------------------------------

getSplitFunFACT <- function(){}
