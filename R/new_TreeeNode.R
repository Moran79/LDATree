new_TreeeNode <- function(xCurrent,
                          response,
                          idxCol,
                          idxRow,
                          prior,
                          weights,
                          misClassCost,
                          currentLevel,
                          currentIndex,
                          parentIndex,
                          misReference,
                          nodeModel) {
  Nj = table(response) # overall proportion, for prior calculation
  responseCurrent <- response[idxRow]
  posterior = getMode(v = responseCurrent, prior = prior / Nj, posterior = TRUE) # p(j,t) = prior * Njt / Nj

  if (nodeModel == "mode") {
    nodePredict <- as.factor(names(which.max(posterior)))
    resubPredict <- factor(rep(nodePredict, length(responseCurrent)), levels = levels(responseCurrent))
  } else if (nodeModel == "LDA") {
    #> LDA will stop if empty level exists
    #> But response level can not be dropped
    datCombined = data.frame(droplevels(xCurrent), responseCurrent)

    ### New code
    # chi_stat = apply(datCombined,2,function(x) var_select_LDATree(xCurrent,responseCurrent, length(responseCurrent), sum(table(responseCurrent) != 0)))
    # topFiveIdx <- order(chi_stat, decreasing = TRUE)[seq_len(min(length(chi_stat), 20))]
    # nodePredict <<- MASS::lda(responseCurrent~., data = datCombined[,topFiveIdx, drop = FALSE], prior = as.vector(posterior))

    # write.csv(datCombined, paste(currentIndex,".csv",sep = ""), row.names = FALSE)

    # nodePredict <<- MASS::lda(responseCurrent~., data = datCombined, prior = as.vector(posterior), tol = 1e-8)
    # resubPredict <- predict(object = nodePredict, newdata = datCombined)$class
    nodePredict <<- ldaGSVD(responseCurrent~., data = datCombined, prior = as.vector(posterior))
    resubPredict <- predict.ldaGSVD(object = nodePredict, newdata = datCombined)
    resubPredict <- factor(resubPredict, levels = levels(response))
  }

  # print(resubPredict[1:5])
  # print(responseCurrent[1:5])
  # print(idxRow)
  # print(table(resubPredict, responseCurrent))
  # print(dim(misClassCost))

  # browser()
  currentTreeeNode <- list(
    currentIndex = currentIndex,
    currentLevel = currentLevel,
    idxRow = idxRow,
    idxCol = idxCol,
    posterior = posterior,
    currentLoss = sum(table(resubPredict, responseCurrent) * misClassCost %*% diag(posterior)),
    accuracy = mean(resubPredict == responseCurrent),
    proportions = table(responseCurrent), #
    parent = parentIndex,
    children = c(), # is.null to check terminal nodes
    # offsprings = c(), # all descendent
    misReference = misReference,
    # colClass = ifelse(getNumFlag(dat),"numeric","factor"),
    # criteria = NA, # 用来打印在output tree上面
    # splitIdx = NA, # 用来记录是用哪一个变量进行split，也可以用来判断是否为叶子结点
    # LDscaling = NA,
    splitCut = NA, # Splitting criteria
    # splitMissingAction = NA, # 把NA分到左面或者右面, 1 代表左面，0代表右面。NA代表训练时没有NA
    # linear_split_trans = NA, # for linear combination split, record the reverse function for prediction
    # alpha = NA, # for CART pruning
    # childrenLoss = NA, # R(T): for CV pruning
    nodeModel = nodeModel,
    nodePredict = nodePredict # predict Function
  )

  # Set the name for the class
  class(currentTreeeNode) <- "TreeeNode"
  return(currentTreeeNode)
}
