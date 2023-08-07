#' Title
#'
#' @param formula
#' @param data
#' @param splitMethod
#' @param pruneMethod
#' @param prior
#' @param weights
#' @param maxTreeLevel
#' @param minNodeSize
#' @param numberOfPrune
#' @param misClassCost
#' @param missingMethod
#' @param randomSeed
#'
#' @return
#' @export
#'
#' @examples
Treee <- function(formula,
                  data,
                  splitMethod = 'LDScores',
                  pruneMethod = 'CV',
                  prior = NULL,
                  weights = NULL,
                  maxTreeLevel = 4,
                  minNodeSize = NULL,
                  numberOfPrune = 20,
                  misClassCost = NULL,
                  missingMethod = c("meanFlag", "newLevel"),
                  randomSeed = NULL){
  ### Arguments ###
  #> x: data matrix without response
  #> response: vector with the same length of x
  #> splitMethod: univariate / LDScores
  #> pruneMethod: CV / none
  #> misClassCost: matrix C(i|j) misclassfied j as i
  #> missingMethod: for numerical / categorical variables, respectively
  #> randomSeed: used in pruning
  #>
  #> 改写function的内容，将不重要的东西放入control

  message("Hi, my friend. The Treee is growing...")

  # Data & Parameter Pre-processing -------------------------------------------------------------------

  # 这里添加更多的parameter error trap
  # 这里添加更多的data pre-processing：比如说把categorical变成factor


  #> droplevels is necessary, since empty response level occurs during train/test split
  #> covariates can have empty levels as well
  modelFrame <- droplevels(model.frame(formula, data, na.action = "na.pass"))
  modelFrame <- modelFrame[which(!is.na(modelFrame[,1])),] # remove NAs in response

  response <- as.factor(modelFrame[,1])
  x <- modelFrame[,-1, drop = FALSE]
  if(is.null(minNodeSize)) {minNodeSize <-  min(dim(data)[1] %/% 100, nlevels(response) + 1)}

  prior <- checkPrior(prior, response)
  misClassCost <- checkMisClassCost(misClassCost, response)

  # Build Single Tree ----------------------------------------------------------------
  treeeNow = new_SingleTreee(x = x,
                             response = response,
                             idxCol = seq_len(ncol(x)),
                             idxRow = seq_len(nrow(x)),
                             splitMethod = splitMethod,
                             prior = prior,
                             weights = weights,
                             maxTreeLevel = maxTreeLevel,
                             minNodeSize = minNodeSize,
                             misClassCost = misClassCost,
                             missingMethod = missingMethod)

  message(paste('The LDA tree is completed. For now, it has', length(treeeNow), 'nodes.\n'))

  finalTreee <- structure(list(treee =  treeeNow,
                               missingMethod = missingMethod), class = "Treee")


  # Pruning -----------------------------------------------------------------
  pruneMethod <- match.arg(pruneMethod, c("CV", "none"))
  if(pruneMethod == "CV" & length(treeeNow) > 1){
    message('0.0 Now, I am about to start the CV !!!')

    pruningOutput <- prune(oldTreee = treeeNow,
                           x = x,
                           response = response,
                           idxCol = seq_len(ncol(x)),
                           idxRow = seq_len(nrow(x)),
                           splitMethod = splitMethod,
                           prior = prior,
                           weights = weights,
                           maxTreeLevel = maxTreeLevel,
                           minNodeSize = minNodeSize,
                           numberOfPrune = numberOfPrune,
                           misClassCost = misClassCost,
                           missingMethod = missingMethod,
                           randomSeed = randomSeed)

    # Add something to the finalTreee
    finalTreee$treee <- pruningOutput$treeeNew
    finalTreee$CV_Table <- pruningOutput$CV_Table
    finalTreee$savedGrove <- pruningOutput$savedGrove

    message('The pruned tree is completed. Wish you a good day :-) \n')
  }

  return(finalTreee)
}
