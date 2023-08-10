#' Title
#'
#' @param formula
#' @param data
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param object
#' @param newdata
#' @param type
#'
#' @return
#' @export
#'
#' @examples
predict.ldaGSVD <- function(object, newdata, type = c("response", "prob")){
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

#' Title
#'
#' @param fit
#'
#' @return
#' @export
#'
#' @examples
print.ldaGSVD <- function(fit){
  message("\nObserved proportions of groups:")
  print(fit$prior)
  message("\nGroup means of LD scores:")
  print(fit$groupMeans)
  message("\nScaling (coefficients) of LD scores:")
  print(fit$scaling)
  invisible(fit)
}


