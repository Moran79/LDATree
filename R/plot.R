#' Plot a Treee object
#'
#' Provide a diagram of the whole tree structure or a scatter/density plot for a
#' specific tree node.
#'
#' @section Overall tree structure:
#'
#'   A full tree diagram (via the R package [visNetwork]) is shown if `node` is
#'   not provided (default is `-1`). The color shows the most common (plurality)
#'   class inside each node. The size of each terminal node is based on its
#'   relative sample size. Under every node, you see the plurality class, the
#'   fraction of the correctly predicted training sample vs. the node's sample
#'   size, and the node index, respectively. When you click on the node, an
#'   information panel with more details will appear.
#'
#' @section Individual plot for each node:
#'
#'   The node index and the original training data are required to return a more
#'   detailed plot within a specific node. The density plot will be provided
#'   when only two levels are left for the response variable in a node (like in
#'   a binary classification problem). Samples are projected down to their first
#'   linear discriminant scores (LD1). A scatter plot will be provided if a node
#'   contains more than two classes. Samples are projected down to their first
#'   and second linear discriminant scores.
#'
#' @param x a fitted model object of class `Treee`, which is assumed to be the
#'   result of the [Treee()] function.
#' @param data the original data you used to fit the `Treee` object if you want
#'   the individual plot for each node. Otherwise, you can leave this parameter
#'   blank if you only need the overall tree structure diagram.
#' @param node the node index that you are interested in. By default, it is set
#'   to `-1` and the overall tree structure is drawn.
#' @param ... further arguments passed to or from other methods.
#'
#' @returns For overall tree structure (`node = -1`), A figure of class
#'   `visNetwork` is drawn. Otherwise, a figure of class `ggplot` is drawn.
#'
#' @export
#'
#' @examples
#' fit <- Treee(Species~., data = iris)
#' # plot the overall tree
#' plot(fit)
#' # plot a certain node
#' plot(fit, iris, node = 1)
plot.Treee <- function(tree, datX, response, node = -1, ...){
  treeeOutput <- tree
  if(node>0){
    if(missing(datX) | missing(response)) stop("Please input the orginal training data for nodewise LDA plots")
    if(treeeOutput$treee[[node]]$nodeModel == "mode") return(paste("Every observation in this node is predicted to be", treeeOutput$treee[[node]]$nodePredict))
    # Get the data ready, impute the NAs (if any)
    response <- as.factor(response)
    newX <- getDataInShape(data = datX[treeeOutput$treee[[node]]$idxRow,,drop = FALSE], missingReference = treeeOutput$treee[[node]]$misReference)
    colorIdx <- match(names(treeeOutput$treee[[node]]$proportions), levels(response))

    plotLDA2d(ldaModel = treeeOutput$treee[[node]]$nodePredict,
              data = cbind.data.frame(response = response[treeeOutput$treee[[node]]$idxRow], newX),
              node = node,
              colorManual = scales::hue_pal()(nlevels(response))[colorIdx])
  }else{ # default overall plot
    plot(treeeOutput$treee)
  }
}


#' @export
plot.SingleTreee <- function(x, ...){
  idTransVec <- seq_along(x)
  nodes <- do.call(rbind, sapply(x, function(treeeNode) nodesHelper(treeeNode = treeeNode, idTransVec = idTransVec),simplify = FALSE))
  edges <- do.call(rbind, sapply(x, edgesHelper,simplify = FALSE))
  p <- visNetwork::visNetwork(nodes, edges, width = "100%", height = "600px")%>%
    visNetwork::visNodes(shape = 'dot')%>%
    visNetwork::visHierarchicalLayout(levelSeparation = 100)%>%
    visNetwork::visLegend(width = 0.1, position = "right", main = "Group")%>%
    visNetwork::visInteraction(dragNodes = FALSE,
                   dragView = TRUE,
                   zoomView = TRUE)

  # Change the color manual
  colorManual = scales::hue_pal()(length(x[[1]]$proportions))
  for(i in seq_along(colorManual)){
    p <- p %>% visNetwork::visGroups(groupname = names(x[[1]]$proportions)[i], color = colorManual[i])
  }

  return(p)
}

infoClickSingle <- function(treeeNode, idTransVec){
  line1 = '#### Information Panel ####'
  line2 = paste('</br>Current Node Index:', idTransVec[treeeNode$currentIndex])
  line3 = paste('</br>There are', length(treeeNode$idxRow), 'data in this node')
  # line4 = paste('</br>The proportion of', paste(names(treeeNode$proportions), collapse = ', '),'are',
  #               paste(sprintf("%.1f%%", treeeNode$proportions / length(treeeNode$idxRow) * 100), collapse = ', '))
  # line4 = paste('</br>', length(treeeNode$idxRow) - treeeNode$currentLoss, 'of them are correctly classified')
  line5 = paste('</br>The resubstitution acc is ', round(treeeNode$accuracy,3))
  line5.5 = paste('</br>Plurality class (', round(max(treeeNode$proportions) / sum(treeeNode$proportions),4)*100, '%) is ', names(sort(treeeNode$proportions, decreasing = TRUE))[1], sep = "")
  line6 = paste('</br>The model in this node is ', treeeNode$nodeModel)
  line6.3 = paste("</br>Pillai's trace is ", tryCatch({treeeNode$nodePredict$statPillai}, error = function(e) {NULL}))
  line6.5 = paste('</br>p value is ', tryCatch({treeeNode$nodePredict$pValue}, error = function(e) {NULL}))
  line6.6 = paste('</br>predGini is ', tryCatch({treeeNode$nodePredict$predGini}, error = function(e) {NULL}))
  line6.7 = paste('</br>alpha is ', tryCatch({treeeNode$alpha}, error = function(e) {NULL}))
  line7 = paste('</br>stopFlag is ', treeeNode$stopFlag)
  return(paste(line1,line2,line3,line5,line5.5,line6,line6.3,line6.5,line6.6,line6.7,line7))
}


nodesHelper <- function(treeeNode, idTransVec){
  terminalFlag <- is.null(treeeNode$children)
  id = treeeNode$currentIndex
  title = infoClickSingle(treeeNode = treeeNode, idTransVec = idTransVec) # Show when you click
  value = ifelse(terminalFlag, log(length(treeeNode$idxRow)), 2) # node size
  level = treeeNode$currentLevel
  group = names(sort(treeeNode$proportions, decreasing = TRUE))[1]
  label = paste(# group, # paste(treeeNode$proportions, collapse = ' / '),
                paste(length(treeeNode$idxRow) - treeeNode$currentLoss, length(treeeNode$idxRow), sep = ' / '),
                # treeeNode$currentLoss,
                paste('Node', idTransVec[id]),sep = "\n")
                # paste('alpha:', treeeNode$alpha)
                # paste('Tnodes:', paste(treeeNode$offsprings, collapse = "/")),
                # paste('Pruned:', treeeNode$pruned)
  return(data.frame(id, title, value, level, group, label, shadow = TRUE))
}


edgesHelper <- function(treeeNode){
  terminalFlag <- is.null(treeeNode$children)
  if(!terminalFlag){
    return(data.frame(from = treeeNode$currentIndex, to = treeeNode$children))
  }
}

plotLDA2d <- function(ldaModel, data, node, colorManual){
  LD1 <- LD2 <- response <- NULL # walk around the binding error in R CMD check
  # browser()
  if(dim(ldaModel$scaling)[2] == 1){
    # Only one LD is available, draw the histogram
    datCombined <- cbind.data.frame(response = data$response, LD1 = getLDscores(modelLDA = ldaModel, data = data, nScores = 1))
    estimatedPrior <- table(datCombined$response) / length(datCombined$response)
    estimatedPrior <- estimatedPrior[which(estimatedPrior != 0)] # some classes are not available
    datPlot <- do.call(rbind, lapply(seq_along(estimatedPrior), function(i) cbind(with(density(datCombined$LD1[datCombined$response == names(estimatedPrior)[i]]), data.frame(LD1 = x, density = y * estimatedPrior[i])), response = names(estimatedPrior)[i])))
    datPlot$response <- factor(datPlot$response, levels = names(estimatedPrior))
    p <- ggplot2::ggplot(data = datPlot)+
      ggplot2::geom_line(ggplot2::aes(x = LD1, y = density, color = response))+
      ggplot2::geom_ribbon(ggplot2::aes(x = LD1, ymin = 0, ymax = density, fill = response), alpha = 0.5)+
      ggplot2::scale_fill_manual(values = colorManual)+
      ggplot2::theme_bw()+
      ggplot2::labs(title = "Density plot of LD1", subtitle = paste("Node:",node))
  }else{
    LDscores <- getLDscores(modelLDA = ldaModel, data = data, nScores = 2)
    datPlot <- cbind.data.frame(response = data$response, LDscores)
    p <- ggplot2::ggplot(data = datPlot)+
      ggplot2::geom_point(ggplot2::aes(x = LD1, y = LD2, color = response), alpha = 0.7)+
      ggplot2::scale_color_manual(values = colorManual)+
      ggplot2::theme_bw()+
      ggplot2::labs(title = "Scatter plot by first two LDscores", subtitle = paste("Node:",node))
  }
  return(p)
}



