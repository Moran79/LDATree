#' Classification Trees with Uncorrelated Linear Discriminant Analysis Terminal
#' Nodes
#'
#' This function fits a classification tree where each node has a Uncorrelated
#' Linear Discriminant Analysis (ULDA) model. It can also handle missing values
#' and perform downsampling. The resulting tree can be pruned either through
#' pre-pruning or post-pruning methods.
#'
#' @param datX A data frame of predictor variables.
#' @param response A vector of response values corresponding to `datX`.
#' @param ldaType A character string specifying the type of LDA to use. Options
#'   are `"forward"` for forward ULDA or `"all"` for full ULDA. Default is
#'   `"forward"`.
#' @param nodeModel A character string specifying the type of model used in each
#'   node. Options are `"ULDA"` for Uncorrelated LDA, or `"mode"` for predicting
#'   based on the most frequent class. Default is `"ULDA"`.
#' @param pruneMethod A character string specifying the pruning method. `"pre"`
#'   performs pre-pruning based on p-value thresholds, and `"post"` performs
#'   cross-validation-based post-pruning. Default is `"pre"`.
#' @param numberOfPruning An integer specifying the number of folds for
#'   cross-validation during post-pruning. Default is `10`.
#' @param maxTreeLevel An integer controlling the maximum depth of the tree.
#'   Increasing this value allows for deeper trees with more nodes. Default is
#'   `20`.
#' @param minNodeSize An integer controlling the minimum number of samples
#'   required in a node. Setting a higher value may lead to earlier stopping and
#'   smaller trees. If not specified, it defaults to one plus the number of
#'   response classes.
#' @param pThreshold A numeric value used as a threshold for pre-pruning based
#'   on p-values. Lower values result in more conservative trees. If not
#'   specified, defaults to `0.01` for pre-pruning and `0.6` for post-pruning.
#' @param prior A numeric vector of prior probabilities for each class. If
#'   `NULL`, the prior is automatically calculated from the data.
#' @param misClassCost A square matrix \eqn{C}, where each element \eqn{C_{ij}}
#'   represents the cost of classifying an observation into class \eqn{i} given
#'   that it truly belongs to class \eqn{j}. If `NULL`, a default matrix with
#'   equal misclassification costs for all class pairs is used. Default is
#'   `NULL`.
#' @param missingMethod A character string specifying how missing values should
#'   be handled. Options include `'mean'`, `'median'`, `'meanFlag'`,
#'   `'medianFlag'` for numerical variables, and `'mode'`, `'modeFlag'`,
#'   `'newLevel'` for factor variables. `'Flag'` options indicate whether a
#'   missing flag is added, while `'newLevel'` replaces missing values with a
#'   new factor level.
#' @param kSample An integer specifying the number of samples to use for
#'   downsampling during tree construction. Set to `-1` to disable downsampling.
#' @param verbose A logical value. If `TRUE`, progress messages and detailed
#'   output are printed during tree construction and pruning. Default is
#'   `FALSE`.
#'
#' @returns An object of class `Treee` containing the fitted tree, which is a
#'   list of nodes, each an object of class `TreeeNode`. Each `TreeeNode`
#'   contains:
#' * `currentIndex`: The node index in the tree.
#' * `currentLevel`: The depth of the current node in the tree.
#' * `idxRow`, `idxCol`: Row and column indices indicating which part of the original data was used for this node.
#' * `currentLoss`: The training error for this node.
#' * `accuracy`: The training accuracy for this node.
#' * `stopInfo`: Information on why the node stopped growing.
#' * `proportions`: The observed frequency of each class in this node.
#' * `prior`: The (adjusted) class prior probabilities used for ULDA or mode prediction.
#' * `misClassCost`: The misclassification cost matrix used in this node.
#' * `parent`: The index of the parent node.
#' * `children`: A vector of indices of this node’s direct children.
#' * `splitFun`: The splitting function used for this node.
#' * `nodeModel`: Indicates the model fitted at the node (`'ULDA'` or `'mode'`).
#' * `nodePredict`: The fitted model at the node, either a ULDA object or the plurality class.
#' * `alpha`: The p-value from a two-sample t-test used to evaluate the strength of the split.
#' * `childrenTerminal`: A vector of indices representing the terminal nodes that are descendants of this node.
#' * `childrenTerminalLoss`: The total training error accumulated from all nodes listed in `childrenTerminal`.
#'
#' @export
#'
#' @references Wang, S. (2024). FoLDTree: A ULDA-Based Decision Tree Framework
#'   for Efficient Oblique Splits and Feature Selection. \emph{arXiv preprint
#'   arXiv:2410.23147}. Available at \url{https://arxiv.org/abs/2410.23147}.
#'
#'   Wang, S. (2024). A New Forward Discriminant Analysis Framework Based On
#'   Pillai's Trace and ULDA. \emph{arXiv preprint arXiv:2409.03136}. Available
#'   at \url{https://arxiv.org/abs/2409.03136}.
#'
#' @examples
#' fit <- Treee(datX = iris[, -5], response = iris[, 5], verbose = FALSE)
#' # Use cross-validation to prune the tree
#' fitCV <- Treee(datX = iris[, -5], response = iris[, 5], pruneMethod = "post", verbose = FALSE)
#' head(predict(fit, iris)) # prediction
#' plot(fit) # plot the overall tree
#' plot(fit, datX = iris[, -5], response = iris[, 5], node = 1) # plot a certain node
Treee <- function(datX,
                  response,
                  ldaType = c("forward", "all"),
                  nodeModel = c("ULDA", "mode"),
                  pruneMethod = c("pre", "post"),
                  numberOfPruning = 10L,
                  maxTreeLevel = 20L,
                  minNodeSize = NULL,
                  pThreshold = NULL,
                  prior = NULL,
                  misClassCost = NULL,
                  missingMethod = c("medianFlag", "newLevel"),
                  kSample = -1,
                  verbose = TRUE){

  # Standardize the Arguments -----------------------------------------------

  datX <- data.frame(datX) # change to data.frame, remove the potential tibble attribute
  for(i in seq_along(datX)){ # remove ordered factors
    if(inherits(datX[,i], c("ordered"))) class(datX[,i]) <- "factor"
  }
  # Remove NAs in the response
  idxNonNA <- which(!is.na(response)); response <- droplevels(factor(response[idxNonNA], ordered = FALSE))
  datX <- datX[idxNonNA, , drop = FALSE]

  ldaType <- match.arg(ldaType, c("forward", "all"))
  nodeModel <- match.arg(nodeModel, c("ULDA", "mode"))
  pruneMethod <- match.arg(pruneMethod, c("pre", "post"))
  if(is.null(minNodeSize)) minNodeSize <- nlevels(response) + 1 # minNodeSize: If not specified, set to J+1
  if(is.null(pThreshold)) pThreshold <- ifelse(pruneMethod == "pre", 0.01, 0.6)
  priorAndMisClassCost <- updatePriorAndMisClassCost(prior = prior, misClassCost = misClassCost, response = response, insideNode = FALSE)
  prior <- priorAndMisClassCost$prior; misClassCost <- priorAndMisClassCost$misClassCost

  # Build Different Trees ---------------------------------------------------

  treeeNow = new_SingleTreee(datX = datX,
                             response = response,
                             ldaType = ldaType,
                             nodeModel = nodeModel,
                             maxTreeLevel = maxTreeLevel,
                             minNodeSize = minNodeSize,
                             pThreshold = pThreshold,
                             prior = prior,
                             misClassCost = misClassCost,
                             missingMethod = missingMethod,
                             kSample = kSample,
                             verbose = verbose)

  if(verbose) cat(paste('\nThe pre-pruned LDA tree is completed. It has', length(treeeNow), 'nodes.\n'))

  # Pruning -----------------------------------------------------------------

  if(pruneMethod == "post" & length(treeeNow) > 1){
    pruningOutput <- prune(oldTreee = treeeNow,
                           datX = datX,
                           response = response,
                           ldaType = ldaType,
                           nodeModel = nodeModel,
                           numberOfPruning = numberOfPruning,
                           maxTreeLevel = maxTreeLevel,
                           minNodeSize = minNodeSize,
                           pThreshold = pThreshold,
                           prior = prior,
                           misClassCost = misClassCost,
                           missingMethod = missingMethod,
                           kSample = kSample,
                           verbose = verbose)

    treeeNow <- pruningOutput$treeeNew
    attr(treeeNow, "CV_Table") <- pruningOutput$CV_Table
    if(verbose) cat(paste('\nThe post-pruned tree is completed. It has', length(treeeNow), 'nodes.\n'))
  }

  return(treeeNow)
}
