#' Title
#'
#' @param formula
#' @param data
#' @param missingMethod
#' @param splitMethod
#' @param pruneMethod
#' @param numberOfPruning
#' @param maxTreeLevel
#' @param minNodeSize
#'
#' @return
#' @export
#'
#' @examples
Treee <- function(formula,
                  data,
                  missingMethod = c("meanFlag", "newLevel"),
                  splitMethod = 'LDScores',
                  pruneMethod = 'none',
                  numberOfPruning = 10,
                  maxTreeLevel = 4,
                  minNodeSize = NULL){
  ### Arguments ###
  #> splitMethod: univariate / LDScores
  #> pruneMethod: CV / none
  #> missingMethod: for numerical / categorical variables, respectively
  #> minNodeSize: 1% of data / J + 1

  # Data & Parameter Pre-processing -------------------------------------------------------------------

  #> droplevels is necessary, since empty response level occurs during train/test split
  #> covariates can have empty levels as well
  modelFrame <- droplevels(model.frame(formula, data, na.action = "na.pass"))
  modelFrame <- modelFrame[which(!is.na(modelFrame[,1])), , drop = FALSE] # remove NAs in response

  response <- as.factor(modelFrame[,1])
  x <- modelFrame[,-1, drop = FALSE]

  # minNodeSize: It is too arbitrary if based on % of the sample size
  if(is.null(minNodeSize)) minNodeSize <- nlevels(response) + 1

  # Build Single Tree ----------------------------------------------------------------
  treeeNow = new_SingleTreee(x = x,
                             response = response,
                             idxCol = seq_len(ncol(x)),
                             idxRow = seq_len(nrow(x)),
                             missingMethod = missingMethod,
                             splitMethod = splitMethod,
                             maxTreeLevel = maxTreeLevel,
                             minNodeSize = minNodeSize)

  message(paste('The unpruned LDA tree is completed. For now, it has', length(treeeNow), 'nodes.\n'))

  finalTreee <- structure(list(treee =  treeeNow,
                               missingMethod = missingMethod), class = "Treee")


  # Pruning -----------------------------------------------------------------
  pruneMethod <- match.arg(pruneMethod, c("CV", "none"))
  if(pruneMethod == "CV" & length(treeeNow) > 1){
    message('Pruning has started...')

    pruningOutput <- prune(oldTreee = treeeNow,
                           x = x,
                           response = response,
                           idxCol = seq_len(ncol(x)),
                           idxRow = seq_len(nrow(x)),
                           splitMethod = splitMethod,
                           maxTreeLevel = maxTreeLevel,
                           minNodeSize = minNodeSize,
                           numberOfPruning = numberOfPruning,
                           missingMethod = missingMethod)

    # Add something to the finalTreee
    finalTreee$treee <- pruningOutput$treeeNew
    finalTreee$CV_Table <- pruningOutput$CV_Table
    finalTreee$savedGrove <- pruningOutput$savedGrove

    message('The pruned tree is completed. Wish you a good day :-) \n')
  }

  return(finalTreee)
}
