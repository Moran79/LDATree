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
Treee <- function(formula,
                  data,
                  treeType = c("single", "forest", "boosting"),
                  ldaType = c("step", "all"),
                  forest = FALSE,
                  missingMethod = c("meanFlag", "newLevel"),
                  splitMethod = c("FACT", "groupMean", "mixed"),
                  nTree = 20,
                  maxTreeLevel = 100,
                  minNodeSize = NULL,
                  trainErrorCap = c("numOfNodes", "none", "zero"),
                  verbose = TRUE,
                  datTest = NULL){
  ### Arguments ###
  splitMethod <- match.arg(splitMethod, c("FACT", "groupMean", "mixed"))
  ldaType <- match.arg(ldaType, c("step", "all"))
  trainErrorCap <- match.arg(trainErrorCap, c("numOfNodes", "none", "zero"))
  treeType <- match.arg(treeType, c("single", "forest", "boosting"))

  # Data & Parameter Pre-processing -------------------------------------------------------------------
  dataProcessed <- extractXnResponse(formula, data)
  x <- dataProcessed$x
  response <- dataProcessed$response
  rm(dataProcessed)

  # minNodeSize: It is too arbitrary if based on % of the sample size
  if(is.null(minNodeSize)) minNodeSize <- nlevels(response) + 1


  # Build Different Trees ---------------------------------------------------


  if(treeType == "single"){
    treeeNow = new_SingleTreee(x = x,
                               response = response,
                               treeType = treeType,
                               ldaType = ldaType,
                               forest = forest,
                               missingMethod = missingMethod,
                               splitMethod = splitMethod,
                               maxTreeLevel = maxTreeLevel,
                               minNodeSize = minNodeSize,
                               trainErrorCap = trainErrorCap,
                               verbose = verbose,
                               datTest = datTest)

    if(verbose) cat(paste('The unpruned LDA tree is completed. It has', length(treeeNow), 'nodes.\n'))

    finalTreee <- structure(list(formula = formula,
                                 treee =  treeeNow,
                                 treeType = treeType,
                                 missingMethod = missingMethod), class = "Treee")
  }else if(treeType == "forest"){
    forestNow <- replicate(nTree, new_SingleTreee(x = x,
                                                  response = response,
                                                  treeType = treeType,
                                                  ldaType = ldaType,
                                                  forest = forest,
                                                  missingMethod = missingMethod,
                                                  splitMethod = splitMethod,
                                                  maxTreeLevel = maxTreeLevel,
                                                  minNodeSize = minNodeSize,
                                                  trainErrorCap = trainErrorCap,
                                                  verbose = verbose), simplify = FALSE)
    class(forestNow) <- "ForestTreee"
    finalTreee <- structure(list(formula = formula,
                                 forest =  forestNow,
                                 treeType = treeType,
                                 missingMethod = missingMethod), class = "Treee")
  }
  return(finalTreee)
}
