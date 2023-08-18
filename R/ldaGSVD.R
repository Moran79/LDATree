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
#' @inheritParams Treee
#' @param data a data frame that contains both predictors and the response.
#'   Missing values are NOT allowed.
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
ldaGSVD <- function(formula, data){
  # response <- as.factor(data[[all.vars(formula)[1]]])
  modelFrame <- model.frame(formula, data, na.action = "na.fail")
  Terms <- terms(modelFrame)
  response <- droplevels(as.factor(modelFrame[,1])) # some levels are branched out
  prior <- table(response, dnn = NULL) / length(response) # estimated prior
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
  unitSD <- pmin(diag(sqrt((length(response) - nlevels(response)) / abs(1 - fitSVDp$d^2 + 1e-15)), nrow = rankAll),1e15) # Scale to unit var
  scalingFinal <- (fitSVD$v[,seq_len(rankT), drop = FALSE] %*% diag(1 / fitSVD$d[seq_len(rankT)], nrow = rankT) %*% fitSVDp$v)[,seq_len(rankAll), drop = FALSE] %*% unitSD
  rownames(scalingFinal) <- cnames

  groupMeans <- groupMeans %*% scalingFinal
  rownames(groupMeans) <- levels(response)
  colnames(groupMeans) <- colnames(scalingFinal) <- paste("LD", seq_len(ncol(groupMeans)), sep = "")

  res <- list(scaling = scalingFinal, formula = formula, terms = Terms, prior = prior,
              groupMeans = groupMeans, xlevels = .getXlevels(Terms, modelFrame))
  class(res) <- "ldaGSVD"
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
  return(rownames(object$groupMeans)[max.col(posterior)])
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


