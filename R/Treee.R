#' Classification trees with Linear Discriminant Analysis terminal nodes
#'
#' Fit an LDATree model.
#'
#' Unlike other classification trees, LDATree integrates LDA throughout the
#' entire tree-growing process. Here's a breakdown of its distinctive features:
#' * The tree searches for the best binary split based on sample quantiles of the first linear discriminant score.
#'
#' * An LDA/GSVD model is fitted for each terminal node (For more details, refer to [ldaGSVD()]).
#'
#' * Missing values can be imputed using the mean, median, or mode, with optional missing flags available.
#'
#' * By default, the tree employs a direct-stopping rule. However, cross-validation using the alpha-pruning from CART is also provided.
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
  dataProcessed <- extractXnResponse(formula, data)
  x <- dataProcessed$x
  response <- dataProcessed$response
  rm(dataProcessed)

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

  finalTreee <- structure(list(formula = formula,
                               treee =  treeeNow,
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
