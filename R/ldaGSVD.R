#' Linear Discriminant Analysis using the Generalized Singular Value
#' Decomposition
#'
#' @description `r lifecycle::badge('experimental')` Fit an LDA/GSVD model.
#'
#' @details
#'
#' Traditional Fisher's Linear Discriminant Analysis (LDA) ceases to work when
#' the within-class scatter matrix is singular. The Generalized Singular Value
#' Decomposition (GSVD) is used to address this issue. GSVD simultaneously
#' diagonalizes both the within-class and between-class scatter matrices without
#' the need to invert a singular matrix. This method is believed to be more
#' accurate than PCA-LDA (as in `MASS::lda`) because it also considers the
#' information in the between-class scatter matrix.
#'
#' @param data a data frame that contains both predictors and the response.
#'   Missing values are NOT allowed.
#' @param method default to be all
#'
#' @returns An object of class `ldaGSVD` containing the following components:
#' * `scaling`: a matrix which transforms the training data to LD scores, normalized so that the within-group scatter matrix is proportional to the identity matrix.
#' * `formula`: the formula passed to the [ldaGSVD()]
#' * `terms`: a object of class `terms` derived using the input `formula` and the training data
#' * `prior`: a `table` of the estimated prior probabilities.
#' * `groupMeans`: a matrix that records the group means of the training data on the transformed LD scores.
#' * `xlevels`: a list records the levels of the factor predictors, derived using the input `formula` and the training data
#'
#' @export
#'
#' @references Ye, J., Janardan, R., Park, C. H., & Park, H. (2004). \emph{An
#'   optimization criterion for generalized discriminant analysis on
#'   undersampled problems}. IEEE Transactions on Pattern Analysis and Machine
#'   Intelligence
#'
#'   Howland, P., Jeon, M., & Park, H. (2003). \emph{Structure preserving dimension
#'   reduction for clustered text data based on the generalized singular value
#'   decomposition}. SIAM Journal on Matrix Analysis and Applications
#'
#' @examples
#' fit <- ldaGSVD(Species~., data = iris)
#' # prediction
#' predict(fit,iris)
ldaGSVD <- function(datX,
                    response,
                    method = c("all", "step"),
                    fixNA = TRUE,
                    missingMethod = c("medianFlag", "newLevel"),
                    prior = NULL,
                    misClassCost = NULL,
                    insideTree = FALSE){


  # Pre-processing: Arguments and response ----------------------------------

  method <- match.arg(method, c("all", "step"))
  missingMethod <- c(match.arg(missingMethod[1], c("mean", "median", "meanFlag", "medianFlag")),
                     match.arg(missingMethod[2], c("mode", "modeFlag", "newLevel")))
  stopifnot(!anyNA(response)) # No NAs in the response variable
  response <- droplevels(as.factor(response)) # some levels are branched out

  #> Get the prior
  if(insideTree){
    prior <- getFinalPrior(prior = prior, response = response)
  }else  prior <- checkPriorAndMisClassCost(prior = prior, misClassCost = misClassCost, response = response, internal = FALSE)


  # Pre-processing: Variables -----------------------------------------------

  #> Variable Selection Step: for stepwise LDA only
  if(method == "step"){
    chiStat <- getChiSqStat(datX = datX, y = response)
    idxKeep <- which(chiStat >= qchisq(1 - 0.05/length(chiStat), 1)) # Bonferroni
    if(length(idxKeep) == 0) idxKeep <- seq_len(length(chiStat))
    datX <- datX[, idxKeep, drop = FALSE]
  }

  if(fixNA){
    imputedSummary <- missingFix(data = datX, missingMethod = missingMethod)
    if(anyNA(datX)) datX <- imputedSummary$data
  }

  modelFrame <- model.frame(formula = ~.-1, datX, na.action = "na.fail")
  Terms <- terms(modelFrame)
  m <- scale(model.matrix(Terms, modelFrame)) # constant cols would be changed to NaN in this step
  cnames <- colnames(m)
  currentVarList <- as.vector(which(apply(m, 2, function(x) !any(is.nan(x))))) # remove constant columns and intercept


  if(length(currentVarList) == 0) stop("All variables are constant.")

  if(method == "step"){
    #> Output: currentVarList, which contains indices of the selected variables
    #> RESPECTIVELY in the design matrix, some columns of m might be removed

    stepRes <- stepVarSelByF(m = m, response = response, currentCandidates = currentVarList)

    #> When no variable is selected, use the full model
    #> it might be more time-consuming, but it is better for future LDA split

    if(length(stepRes$currentVarList) != 0){
      currentVarList <- stepRes$currentVarList

      #> modify the design matrix to make it more compact
      #> so that only the selected variables are included in the design matrix,
      #> and eventually make the prediction faster
      selectedVarRawIdx <- unique(sort(attributes(m)$assign[currentVarList])) # MUST be from the modelFrame where the factors are not dummied
      modelFrame <- model.frame(formula = ~.-1, datX[, selectedVarRawIdx, drop = FALSE], na.action = "na.fail")
      Terms <- terms(modelFrame)
      m <- scale(model.matrix(Terms, modelFrame))

      #> select CERTAIN levels of the factor variables, not ALL
      currentVarList <- which(colnames(m) %in% stepRes$stepInfo$var)
    }
  }

  varSD <- attr(m,"scaled:scale")[currentVarList]
  varCenter <- attr(m,"scaled:center")[currentVarList]
  m <- m[,currentVarList, drop = FALSE]

  # Step 1: SVD on the combined matrix H
  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  Hb <- sqrt(tabulate(response)) * groupMeans # grandMean = 0 if scaled

  Hw <- m - groupMeans[response, , drop = FALSE]
  if(diff(dim(m)) < 0){ # More rows than columns
    qrRes <- qrEigen(Hw)
    fitSVD <- svdEigen(rbind(Hb, qrRes$R))
  }else fitSVD <- svdEigen(rbind(Hb, Hw))

  # fitSVD <- saferSVD(rbind(Hb, m - groupMeans[response, , drop = FALSE]))
  rankT <- sum(fitSVD$d >= max(dim(fitSVD$u),dim(fitSVD$v)) * .Machine$double.eps * fitSVD$d[1])

  # Step 2: SVD on the P matrix
  #> The code below can be changed to saferSVD if necessary
  # fitSVDp <- saferSVD(fitSVD$u[seq_len(nlevels(response)), seq_len(rankT), drop = FALSE], nu = 0L)
  fitSVDp <- svdEigen(fitSVD$u[seq_len(nlevels(response)), seq_len(rankT), drop = FALSE])
  rankAll <- min(nlevels(response)-1, rankT) # This is not optimal, but rank(Hb) takes time

  # Fix the variance part
  unitSD <- pmin(diag(sqrt((length(response) - nlevels(response)) / abs(1 - fitSVDp$d^2 + 1e-15)), nrow = rankAll),1e15) # Scale to unit var
  scalingFinal <- (fitSVD$v[,seq_len(rankT), drop = FALSE] %*% diag(1 / fitSVD$d[seq_len(rankT)], nrow = rankT) %*% fitSVDp$v)[,seq_len(rankAll), drop = FALSE] %*% unitSD
  rownames(scalingFinal) <- cnames[currentVarList]

  groupMeans <- groupMeans %*% scalingFinal
  rownames(groupMeans) <- levels(response)
  colnames(groupMeans) <- colnames(scalingFinal) <- paste("LD", seq_len(ncol(groupMeans)), sep = "")

  # Get the test statistics and related p-value
  statPillai <- sum(fitSVDp$d[seq_len(rankAll)]^2)
  #> s & p are changed here, since sometimes design matrix is not of full rank p
  p <- rankT; J <- nlevels(response); N <- nrow(m)
  s <- rankAll; numF <- N-J-p+s; denF <- abs(p-J+1)+s

  #> When numF is non-positive, Pillai = s & training accuracy = 100%
  #> since there always exist a dimension where we can separate every class perfectly
  # pValue <- pf(numF / denF * statPillai / (s - statPillai), df1 = s*denF, df2 = s*numF, lower.tail = F) # the same answer
  pValue <- ifelse(numF > 0, pbeta(1 - statPillai / s, shape1 = numF * s / 2, shape2 = denF * s / 2), 0)
  if(method == "step") pValue <- pValue * length(cnames) # Bonferroni correction

  res <- list(scaling = scalingFinal, terms = Terms, prior = prior,
              groupMeans = groupMeans, xlevels = .getXlevels(Terms, modelFrame),
              varIdx = currentVarList, varSD = varSD, varCenter = varCenter, statPillai = statPillai,
              pValue = pValue)
  if(fixNA) res$misReference <- imputedSummary$ref
  if(method == "step"){
    res$stepInfo = stepRes$stepInfo
    res$stopFlag <- stepRes$stopFlag
  }
  class(res) <- "ldaGSVD"

  # for LDA splitting
  currentP <- unname(table(predict(res, datX)) / dim(datX)[1])
  res$predGini <- 1 - sum(currentP^2)
  return(res)
}



#' Predictions from a fitted ldaGSVD object
#'
#' Prediction of test data using a fitted ldaGSVD object
#'
#' Unlike the original paper, which uses the k-nearest neighbor (k-NN) as the
#' classifier, we use a faster and more straightforward likelihood-based method.
#' One limitation of the traditional likelihood-based method for LDA is that it
#' ceases to work when there are Linear Discriminant (LD) directions with zero
#' variance in the within-class scatter matrix. However, when using LDA/GSVD,
#' all chosen LD directions possess non-zero variance in the between-class
#' scatter matrix. This implies that LD directions with zero variance in the
#' within-class scatter matrix will yield the highest Fisher's ratio. Therefore,
#' to get these directions higher weights, we manually adjust the zero variance
#' to `1e-15` for computational reasons.
#'
#' @param object a fitted model object of class `ldaGSVD`, which is assumed to
#'   be the result of the [ldaGSVD()] function.
#' @param newdata data frame containing the values at which predictions are
#'   required. Missing values are NOT allowed.
#' @param type character string denoting the type of predicted value returned.
#'   The default is to return the predicted class (`type` = 'response'). The
#'   predicted posterior probabilities for each class will be returned if `type`
#'   = 'prob'.
#' @param ... further arguments passed to or from other methods.
#'
#' @return The function returns different values based on the `type`, if
#' * `type = 'response'`: vector of predicted responses.
#' * `type = 'prob'`: a data frame of the posterior probabilities. Each class takes a column.
#' @export
#'
#' @references Ye, J., Janardan, R., Park, C. H., & Park, H. (2004). \emph{An
#'   optimization criterion for generalized discriminant analysis on
#'   undersampled problems}. IEEE Transactions on Pattern Analysis and Machine
#'   Intelligence
#'
#'   Howland, P., Jeon, M., & Park, H. (2003). \emph{Structure preserving dimension
#'   reduction for clustered text data based on the generalized singular value
#'   decomposition}. SIAM Journal on Matrix Analysis and Applications
#'
#' @examples
#' fit <- ldaGSVD(Species~., data = iris)
#' predict(fit,iris)
#' # output prosterior probabilities
#' predict(fit,iris,type = "prob")
predict.ldaGSVD <- function(object, newdata, type = c("response", "prob"), ...){
  type <- match.arg(type, c("response", "prob"))
  # add one extra check for levels of the predictors
  LDscores <- getLDscores(modelLDA = object, data = newdata)
  loglikelihood <- LDscores %*% t(object$groupMeans) + matrix(log(object$prior) - 0.5 * rowSums(object$groupMeans^2), nrow(LDscores), length(object$prior), byrow = TRUE)
  # Computation Optimization 2: Prevent a very large likelihood due to exponential
  likelihood <- exp(loglikelihood - apply(loglikelihood, 1, max))
  posterior <- likelihood / apply(likelihood, 1, sum)
  if(type == "prob") return(posterior)
  return(rownames(object$groupMeans)[max.col(posterior, ties.method = "first")])
}


#' @export
print.ldaGSVD <- function(x, ...){
  cat("\nObserved proportions of groups:\n")
  print(x$prior)
  cat("\n\nGroup means of LD scores:\n")
  print(x$groupMeans)
  cat("\n\nScaling (coefficients) of LD scores:\n")
  print(x$scaling)
  invisible(x)
}


# helper functions --------------------------------------------------------

getPillai <- function(Sw, St){
  #> return -1 if the St is not full rank
  #> sometimes diag has negative value due to instability of R
  #> We ignore those for now
  tryCatch({
    sum(diag(solve(St, St - Sw)))
  }, error = function(e) {-1})
}

selectVar <- function(currentVar, newVar, Sw, St, direction = "forward"){
  #> return the column index
  #> return 0 if all var makes St = 0
  if(direction == "forward"){
    lambdaAll <- sapply(newVar, function(i) getPillai(Sw[c(i,currentVar),c(i,currentVar), drop = FALSE], St[c(i,currentVar),c(i,currentVar), drop = FALSE]))
  }else{
    lambdaAll <- sapply(currentVar, function(i) getPillai(Sw[setdiff(currentVar, i),setdiff(currentVar, i), drop = FALSE], St[setdiff(currentVar, i),setdiff(currentVar, i), drop = FALSE]))
  }
  maxVarIdx <- which(lambdaAll == max(lambdaAll)) # all variables that achieve maximum
  currentVarIdx <- maxVarIdx[1]

  return(list(stopflag = (lambdaAll[currentVarIdx] == -1),
              varIdx = newVar[currentVarIdx],
              statistics = lambdaAll[currentVarIdx],
              maxVarIdx = newVar[maxVarIdx]))
}


stepVarSelByF <- function(m, response, currentCandidates){
  idxOriginal <- currentCandidates
  m <- m[,currentCandidates, drop = FALSE] # all columns should be useful

  groupMeans <- tapply(c(m), list(rep(response, dim(m)[2]), col(m)), function(x) mean(x, na.rm = TRUE))
  mW <- m - groupMeans[response, , drop = FALSE]

  # Initialize
  n = nrow(m); g = nlevels(response); p = 0; currentVarList = c()
  previousPillai <- previousDiff <- previousDiffDiff <- numeric(ncol(m)+1);
  previousDiff[1] <- Inf; diffChecker <- 0
  kRes <- 1; currentCandidates <- seq_len(ncol(m))
  Sw <- St <- matrix(NA, nrow = ncol(m), ncol = ncol(m))
  diag(Sw) <- apply(mW^2,2,sum); diag(St) <- apply(m^2,2,sum)
  stopFlag <- 0

  # Empirical: Calculate the threshold for pillaiToEnter
  pillaiThreshold <- 1 / (1 + (n-g) / (abs(g-2)+1) / qf(1 - 0.1 / (ncol(m)+1), abs(g-2)+1, n-g)) / currentCandidates^(0.25)

  stepInfo <- data.frame(var = character(2*ncol(m)),
                         pillaiToEnter = 0,
                         threhold = pillaiThreshold,
                         pillaiToRemove = 0,
                         pillai = 0)

  #> If n <= g, which means there are too few observations,
  #> we output all columns, and leave that problem to outside function
  if(anyNA(pillaiThreshold)){
    currentVarList <- currentCandidates; currentCandidates <- c()
    stepInfo$var[seq_along(currentVarList)] <- colnames(m)[currentVarList]
    stopFlag <- 4
  }

  # Stepwise selection starts!
  while(length(currentCandidates) != 0){

    nCandidates <- length(currentCandidates)
    p = p + 1
    selectVarInfo <- selectVar(currentVar = currentVarList,
                               newVar = currentCandidates,
                               Sw = Sw,
                               St = St)
    bestVar <- selectVarInfo$varIdx
    if(selectVarInfo$stopflag){ # If St = 0, stop. [Might never happens, since there are other variables to choose]
      stopFlag <- 1
      break
    }

    # get the difference in Pillai's trace
    previousDiff[p+1] <- selectVarInfo$statistics - previousPillai[p]
    previousDiffDiff[p+1] <- previousDiff[p+1] - previousDiff[p]
    diffChecker <- ifelse(abs(previousDiffDiff[p+1]) < 0.001, diffChecker + 1, 0)

    if(previousDiff[p+1] > 10 * previousDiff[p]){ # Correlated variable(s) is included
      currentCandidates <- setdiff(currentCandidates, selectVarInfo$maxVarIdx)
      p <- p - 1; next
    }

    # Check the stopping rule
    if(previousDiff[p+1] < pillaiThreshold[p]){ # If no significant variable selected, stop
      stopFlag <- 2
      break
    }
    if(diffChecker == 10){ # converge
      stopFlag <- 3
      break
    }

    # Add the variable into the model
    previousPillai[p+1] <- selectVarInfo$statistics
    currentVarList <- c(currentVarList, bestVar)
    currentCandidates <- setdiff(currentCandidates, bestVar)
    stepInfo$var[kRes] <- colnames(m)[bestVar]
    stepInfo$pillaiToEnter[kRes] <- previousDiff[p+1]
    stepInfo$pillai[kRes] <- previousPillai[p+1]
    kRes <- kRes + 1

    # Update the Sw and St on the new added column
    Sw[currentCandidates, bestVar] <- Sw[bestVar, currentCandidates] <- as.vector(t(mW[, currentCandidates, drop = FALSE]) %*% mW[,bestVar, drop = FALSE])
    St[currentCandidates, bestVar] <- St[bestVar, currentCandidates] <- as.vector(t(m[, currentCandidates, drop = FALSE]) %*% m[,bestVar, drop = FALSE])
  }

  # Remove the empty rows in the stepInfo if stepLDA does not select all variables
  stepInfo <- stepInfo[seq_along(currentVarList),]

  return(list(currentVarList = idxOriginal[currentVarList], stepInfo = stepInfo, stopFlag = stopFlag))
}



saferSVD <- function(x, ...){
  #> Target for error code 1 from Lapack routine 'dgesdd' non-convergence error
  #> Current solution: Round the design matrix to make approximations,
  #> hopefully this will solve the problem
  #>
  #> The code is a little lengthy, since the variable assignment in tryCatch is tricky
  parList <- list(svdObject = NULL,
                  svdSuccess = FALSE,
                  errorDigits = 16,
                  x = x)
  while (!parList$svdSuccess) {
    parList <- tryCatch({
      parList$svdObject <- svd(parList$x, ...)
      parList$svdSuccess <- TRUE
      parList
    }, error = function(e) {
      if (grepl("error code 1 from Lapack routine 'dgesdd'", e$message)) {
        parList$x <- round(x, digits = parList$errorDigits)
        parList$errorDigits <- parList$errorDigits - 1
        return(parList)
      } else stop(e)
    })
  }
  return(parList$svdObject)
}



