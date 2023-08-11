#' Classification trees with Linear Discriminant Analysis terminal nodes
#'
#' Fit an LDATree model.
#'
#' Unlike other classification trees, LDATree integrates LDA throughout the
#' entire tree-growing process. Here is a breakdown of its distinctive features:
#' * The tree searches for the best binary split based on sample quantiles of the first linear discriminant score.
#'
#' * An LDA/GSVD model is fitted for each terminal node (For more details, refer to [ldaGSVD()]).
#'
#' * Missing values can be imputed using the mean, median, or mode, with optional missing flags available.
#'
#' * By default, the tree employs a direct-stopping rule. However, cross-validation using the alpha-pruning from CART is also provided.
#'
#' @param formula an object of class [stats::formula()], which has the form
#'   `class ~ x1 + x2 + ...`
#' @param data a data frame that contains both predictors and the response.
#'   Missing values are allowed in predictors but not in the response.
#' @param missingMethod Missing value solutions for numerical variables and
#'   factor variables. 'mean', 'median', 'meanFlag', 'medianFlag' are available
#'   for numerical variables. 'mode', 'modeFlag', 'newLevel' are available for
#'   factor variables. The word 'Flag' in the methods indicates whether a
#'   missing flag is added or not. The 'newLevel' method means that all missing
#'   values are replaced with a new level rather than imputing them to another
#'   existing value.
#' @param splitMethod the splitting rule in LDATree growing process. For now,
#'   'LDScores' is the only available option.
#' @param pruneMethod the model selection method in the LDATree growing process,
#'   which controls the size of the tree. By default, it's set to 'none', which
#'   applies a direct stopping rule. Alternatively, 'CV' uses the alpha-pruning
#'   process from CART. Although 'CV' is often more accurate, it can be slower,
#'   especially with large datasets.
#' @param numberOfPruning controls the number of cross-validation in the
#'   pruning. It is 10 by default.
#' @param maxTreeLevel controls the largest tree size possible for either a
#'   direct-stopping tree or a CV-pruned tree. Adding one extra level (depth)
#'   introduces an additional layer of nodes at the bottom of the current tree.
#'   e.g., when the maximum level is 1 (or 2), the maximum tree size is 3 (or
#'   7).
#' @param minNodeSize controls the minimum node size. Think carefully before
#'   changing this value. Setting a large number might result in early stopping
#'   and reduced accuracy. By default, it's set to one plus the number of
#'   classes in the response variable.
#'
#' @returns An object of class 'Treee' containing the following components:
#'
#' * `formula`: the formula passed to the `Treee`
#'
#' * `treee`: A list of all the tree nodes, and each node is an object of class 'TreeeNode'.
#'
#' * `missingMethod`: the missingMethod passed to the `Treee`
#'
#'   An object of class 'TreeeNode' containing the following components:
#'
#' * `currentIndex`
#'
#' * `currentLevel`:
#'
#' * `idxRow`, `idxCol`:
#'
#' * `currentLoss`:
#'
#' * `accuracy`:
#'
#' * `proportions`:
#'
#' * `parent`:
#'
#' * `children`:
#'
#' * `misReference`:
#'
#' * `splitCut`
#'
#' * `nodeModel`
#'
#' * `nodePredict`
#'
#' * `offsprings`
#'
#' * `offspringLoss`
#'
#' * `alpha`
#'
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
