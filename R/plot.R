#' Plot a Treee Object
#'
#' This function visualizes either the entire decision tree or a specific node
#' within the tree. The tree is displayed as an interactive network of nodes and
#' edges, while individual nodes are scatter/density plots using `ggplot2`.
#'
#' @section Overall Tree Structure:
#'
#'   A full tree diagram is displayed using \link[visNetwork]{visNetwork} when `node` is not
#'   specified (the default is `-1`). The color represents the most common
#'   (plurality) class within each node, and the size of each terminal node
#'   reflects its relative sample size. Below each node, the fraction of
#'   correctly predicted training samples and the total sample size for that
#'   node are shown, along with the node index. Clicking on a node opens an
#'   information panel with additional details.
#'
#' @section Individual Node Plot:
#'
#'   To plot a specific node, you must provide the node index along with the
#'   original training predictors (`datX`) and responses (`response`). A scatter
#'   plot is generated if more than one discriminant score is available,
#'   otherwise, a density plot is created. Samples are projected onto their
#'   linear discriminant score(s).
#'
#' @param x A fitted model object of class `Treee`, typically the result of the
#'   [Treee()] function.
#' @param datX A data frame of predictor variables. Required for plotting
#'   individual nodes.
#' @param response A vector of response values. Required for plotting individual
#'   nodes.
#' @param node An integer specifying the node to plot. If `node = -1`, the
#'   entire tree is plotted. Default is `-1`.
#' @param ... Additional arguments passed to the plotting functions.
#'
#' @return A `visNetwork` interactive plot of the decision tree if `node = -1`,
#'   or a `ggplot2` object if a specific node is plotted.
#' @export
#'
#' @examples
#' fit <- Treee(datX = iris[, -5], response = iris[, 5], verbose = FALSE)
#' plot(fit) # plot the overall tree
#' plot(fit, datX = iris, response = iris[, 5], node = 1) # plot a specific node
plot.Treee <- function(x, datX, response, node = -1, ...){
  # Save the color manual, since some classes might be empty during branching
  colorManual = grDevices::hcl.colors(length(x[[1]]$proportions))
  names(colorManual) <- responseLevels <- names(x[[1]]$proportions)

  if(node < 0){ # Overall tree plot
    idTransVec <- seq_along(x)
    nodes <- do.call(rbind, lapply(x, function(treeeNode) nodesHelper(treeeNode = treeeNode, idTransVec = idTransVec)))
    edges <- do.call(rbind, lapply(x, edgesHelper))
    p <- visNetwork::visNetwork(nodes, edges, width = "100%", height = "600px")%>%
      visNetwork::visNodes(shape = 'dot')%>%
      visNetwork::visHierarchicalLayout(levelSeparation = 100)%>%
      visNetwork::visInteraction(dragNodes = FALSE,
                                 dragView = TRUE,
                                 zoomView = TRUE)

    for (i in seq_along(responseLevels)) {
      p <- p %>% visNetwork::visGroups(groupname = responseLevels[i], color = unname(colorManual[i]))
    }

    legend_nodes <- lapply(seq_along(responseLevels), function(i) { # add legends
      list(label = responseLevels[i],
           shape = "dot",
           color = unname(colorManual[i]))
    })
    p <- p %>% visNetwork::visLegend(addNodes = legend_nodes, width = 0.1, useGroups = FALSE, position = "right", main = "Class")
  } else{ # individual node plot
    if(x[[node]]$nodeModel == "mode") return(paste("Every observation in node", node, "is predicted to be", x[[node]]$nodePredict))
    if(missing(datX) || missing(response)) stop("Please input the training X and Y for the nodewise plot")
    response <- droplevels(as.factor(response))
    colorIdx <- match(names(x[[node]]$proportions), levels(response))
    p <- plot(x = x[[node]]$nodePredict,
              datX = datX[x[[node]]$idxRow,,drop = FALSE],
              response = response[x[[node]]$idxRow])
    p$scales$scales <- list() # remove old color palette
    p <- p +
      ggplot2::scale_color_manual(values = colorManual[colorIdx])+
      ggplot2::scale_fill_manual(values = colorManual[colorIdx])+
      ggplot2::labs(caption = paste("Node", node))
  }
  return(p)
}


infoClickSingle <- function(treeeNode, idTransVec){
  line1 = paste('#### Information Panel: Node', idTransVec[treeeNode$currentIndex], '####')
  line2 = paste('</br>There are', length(treeeNode$idxRow), 'data in this node')
  line3 = paste('</br>The resubstitution acc is ', round(treeeNode$accuracy,3))
  line4 = paste('</br>Plurality class (', round(max(treeeNode$proportions) / sum(treeeNode$proportions),4)*100, '%) is ', names(sort(treeeNode$proportions, decreasing = TRUE))[1], sep = "")
  line5 = paste('</br>The model in this node is ', treeeNode$nodeModel)
  line6 = paste('</br>stopInfo:', treeeNode$stopInfo)

  if (treeeNode$nodeModel != "mode") {
    line7 = paste("</br>Pillai's trace is ", round(treeeNode$nodePredict$statPillai, 3))
    line8 = paste('</br>MANOVA p value is ', format(treeeNode$nodePredict$pValue, scientific= TRUE, digits = 4))
    line9 = paste('</br>Gini Index is ', format(treeeNode$nodePredict$predGini, scientific= TRUE, digits = 4))
    line10 = paste('</br>Splitting p value is ', format(treeeNode$alpha, scientific= TRUE, digits = 4))
  } else line7 = line8 = line9 = line10 = ""

  return(paste(line1,line2,line3,line4,line5,line6,line7,line8,line9,line10))
}


nodesHelper <- function(treeeNode, idTransVec){
  terminalFlag <- is.null(treeeNode$children)
  id = treeeNode$currentIndex
  title = infoClickSingle(treeeNode = treeeNode, idTransVec = idTransVec) # Show when you click
  value = ifelse(terminalFlag, log(length(treeeNode$idxRow)), 2) # node size
  level = treeeNode$currentLevel
  group = names(sort(treeeNode$proportions, decreasing = TRUE))[1]
  label = paste(paste(length(treeeNode$idxRow) - treeeNode$currentLoss, length(treeeNode$idxRow), sep = ' / '),
                paste('Node', idTransVec[id]),sep = "\n")
  return(data.frame(id, title, value, level, group, label, shadow = TRUE))
}


edgesHelper <- function(treeeNode){
  terminalFlag <- is.null(treeeNode$children)
  if(!terminalFlag){
    return(data.frame(from = treeeNode$currentIndex, to = treeeNode$children))
  }
}
