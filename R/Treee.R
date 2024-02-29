#' Classification trees with Linear Discriminant Analysis terminal nodes
#'
#' @description `r lifecycle::badge('experimental')` Fit an LDATree model.
#'
#' @details
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
#' @param formula an object of class [formula], which has the form `class ~ x1 +
#'   x2 + ...`
#' @param data a data frame that contains both predictors and the response.
#'   Missing values are allowed in predictors but not in the response.
#' @param missingMethod Missing value solutions for numerical variables and
#'   factor variables. `'mean'`, `'median'`, `'meanFlag'`, `'medianFlag'` are
#'   available for numerical variables. `'mode'`, `'modeFlag'`, `'newLevel'` are
#'   available for factor variables. The word `'Flag'` in the methods indicates
#'   whether a missing flag is added or not. The `'newLevel'` method means that
#'   all missing values are replaced with a new level rather than imputing them
#'   to another existing value.
#' @param splitMethod the splitting rule in LDATree growing process. For now,
#'   `'LDscores'` is the only available option.
#' @param maxTreeLevel controls the largest tree size possible for either a
#'   direct-stopping tree or a CV-pruned tree. Adding one extra level (depth)
#'   introduces an additional layer of nodes at the bottom of the current tree.
#'   e.g., when the maximum level is 1 (or 2), the maximum tree size is 3 (or
#'   7).
#' @param minNodeSize controls the minimum node size. Think carefully before
#'   changing this value. Setting a large number might result in early stopping
#'   and reduced accuracy. By default, it's set to one plus the number of
#'   classes in the response variable.
#' @param verbose a logical. If TRUE, the function provides additional
#'   diagnostic messages or detailed output about its progress or internal
#'   workings. Default is FALSE, where the function runs silently without
#'   additional output.
#'
#' @returns An object of class `Treee` containing the following components:
#' * `formula`: the formula passed to the [Treee()]
#' * `treee`: a list of all the tree nodes, and each node is an object of class `TreeeNode`.
#' * `missingMethod`: the missingMethod passed to the [Treee()]
#'
#'   An object of class `TreeeNode` containing the following components:
#' * `currentIndex`: the node index of the current node
#' * `currentLevel`: the level of the current node in the tree
#' * `idxRow`, `idxCol`: the row and column indices showing which portion of data is used in the current node
#' * `currentLoss`: ?
#' * `accuracy`: the training accuracy of the current node
#' * `stopFlag`: ?
#' * `proportions`: shows the observed frequency for each class
#' * `parent`: the node index of its parent
#' * `children`: the node indices of its direct children (not including its children's children)
#' * `misReference`: a data frame, serves as the reference for missing value imputation
#' * `splitFun`: ?
#' * `nodeModel`: one of `'mode'` or `'LDA'`. It shows the type of predictive model fitted in the current node
#' * `nodePredict`: the fitted predictive model in the current node. It is an object of class `ldaGSVD` if LDA is fitted. If `nodeModel = 'mode'`, then it is a vector of length one, showing the plurality class.
#' * `offsprings`: (available only if `pruneMethod = 'CV'`) showing all terminal descendant nodes of the current node
#' * `alpha`: (available only if `pruneMethod = 'CV'`) the alpha in alpha-pruning from CART
#' @export
#'
#' @examples
#' fit <- Treee(Species~., data = iris)
#' # Use cross-validation to prune the tree
#' fitCV <- Treee(Species~., data = iris, pruneMethod = "CV")
#' # prediction
#' predict(fit,iris)
#' # plot the overall tree
#' plot(fit)
#' # plot a certain node
#' plot(fit, iris, node = 1)
Treee <- function(datX,
                  response,
                  treeType = c("single", "forest"),
                  splitMethod = c("LDA"),
                  ldaType = c("step", "all"),
                  fastTree = FALSE,
                  nodeModel = c("LDA", "mode"),
                  missingMethod = c("meanFlag", "newLevel"),
                  pruneMethod = c("pre", "post", "none"),
                  nTree = 100,
                  maxTreeLevel = 20,
                  minNodeSize = NULL,
                  pThreshold = 0.01,
                  trainErrorDropCap = c("none", "zero", "numOfChildren"),
                  verbose = TRUE,
                  ...){

  # Standardize the Arguments -----------------------------------------------

  response <- as.factor(response) # make it a factor
  treeType <- match.arg(treeType, c("single", "forest"))
  #> The splitMethod option might be deleted in the future if no more methods are implemented
  splitMethod <- match.arg(splitMethod, c("LDA"))
  ldaType <- match.arg(ldaType, c("step", "all"))
  nodeModel <- match.arg(nodeModel, c("LDA", "mode"))
  missingMethod <- c(match.arg(missingMethod[1], c("mean", "median", "meanFlag", "medianFlag")),
                     match.arg(missingMethod[2], c("mode", "modeFlag", "newLevel")))
  pruneMethod <- match.arg(pruneMethod, c("pre", "post", "none"))
  if(is.null(minNodeSize)) minNodeSize <- nlevels(response) + 1 # minNodeSize: If not specified, set to J+1
  trainErrorDropCap <- match.arg(trainErrorDropCap, c("none", "zero", "numOfChildren"))


  # Build Different Trees ---------------------------------------------------

  if(treeType == "single"){
    resNow = new_SingleTreee(datX = datX,
                             response = response,
                             treeType = treeType,
                             splitMethod = splitMethod,
                             ldaType = ldaType,
                             fastTree = fastTree,
                             nodeModel = nodeModel,
                             missingMethod = missingMethod,
                             pruneMethod = pruneMethod,
                             maxTreeLevel = maxTreeLevel,
                             minNodeSize = minNodeSize,
                             pThreshold = pThreshold,
                             trainErrorDropCap = trainErrorDropCap,
                             verbose = verbose,
                             ...)
    #> If we use kStepAhead, "pre" has to be included as well
    if(pruneMethod %in% c("post")) resNow <- pruneByTrainErrDrop(treeeList = resNow,
                                                                 pThreshold = pThreshold,
                                                                 verbose = verbose)
    if(verbose) cat(paste('The LDA tree is completed. It has', length(resNow), 'nodes.\n'))
  } else if(treeType == "forest"){
    resNow <- replicate(nTree, new_SingleTreee(datX = datX,
                                               response = response,
                                               treeType = treeType,
                                               splitMethod = splitMethod,
                                               ldaType = ldaType,
                                               fastTree = fastTree,
                                               nodeModel = nodeModel,
                                               missingMethod = missingMethod,
                                               pruneMethod = "none",
                                               maxTreeLevel = maxTreeLevel,
                                               minNodeSize = minNodeSize,
                                               pThreshold = pThreshold,
                                               trainErrorDropCap = "none", # could be other options
                                               verbose = verbose,
                                               ...), simplify = FALSE)
    class(resNow) <- "ForestTreee"
  }

  finalTreee <- structure(list(resNow =  resNow,
                               treeType = treeType,
                               missingMethod = missingMethod), class = "Treee")
  names(finalTreee)[1] <- ifelse(treeType == "single", "treee", "forest")
  return(finalTreee)
}
