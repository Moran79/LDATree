#' Title
#'
#' @param treeeOutput
#'
#' @return
#' @rawNamespace S3method(plot, Treee)
#' @rawNamespace S3method(plot, SingleTreee)
#'
#' @examples
plot.Treee <- function(treeeOutput){
  plot(treeeOutput$treee)
}

plot.SingleTreee <- function(treeeList){
  # idTransVec <- transIdForBiPlot(treeeList)
  idTransVec <- seq_along(treeeList)

  nodes <- do.call(rbind, sapply(treeeList, function(treeeNode) nodesHelper(treeeNode = treeeNode, idTransVec = idTransVec),simplify = FALSE))
  edges <- do.call(rbind, sapply(treeeList, edgesHelper,simplify = FALSE))

  p <- visNetwork(nodes, edges, width = "100%", height = "600px")%>%
    visNodes(shape = 'dot', color = list(background = "white",
                                         border = "black"))%>%
    visHierarchicalLayout(levelSeparation = 100)%>%
    visLegend(width = 0.1, position = "right", main = "Group")%>%
    visInteraction(dragNodes = FALSE,
                   dragView = TRUE,
                   zoomView = TRUE)
  return(p)
}

infoClickSingle <- function(treeeNode, idTransVec){
  line1 = '#### Information Panel ####'
  line2 = paste('</br>Current Node Index:', idTransVec[treeeNode$currentIndex])
  line3 = paste('</br>There are', length(treeeNode$idxRow), 'data in this node')
  line4 = paste('</br>The proportion of', paste(names(treeeNode$proportions), collapse = ', '),'are',
                paste(sprintf("%.1f%%", treeeNode$proportions / length(treeeNode$idxRow) * 100), collapse = ', '))
  line5 = paste('</br>The accuracy is ', treeeNode$accuracy)
  line6 = paste('</br>The model in this node is ', treeeNode$nodeModel)
  return(paste(line1,line2,line3,line4,line5,line6))
}


nodesHelper <- function(treeeNode, idTransVec){
  terminalFlag <- is.null(treeeNode$children)
  id = treeeNode$currentIndex
  title = infoClickSingle(treeeNode = treeeNode, idTransVec = idTransVec) # Show when you click
  value = ifelse(terminalFlag, log(length(treeeNode$idxRow)), 2)
  level = treeeNode$currentLevel
  group = names(sort(treeeNode$proportions, decreasing = TRUE))[1]
  label = paste(group, paste(treeeNode$proportions, collapse = ' / '),
                paste(length(treeeNode$idxRow) * c(treeeNode$accuracy, 1), collapse = ' / '),
                paste('Node', idTransVec[id]),
                paste('alpha:', treeeNode$alpha),
                paste('Tnodes:', paste(treeeNode$offsprings, collapse = "/")),
                paste('Pruned:', treeeNode$pruned),sep = "\n")
  return(data.frame(id, title,value, level, group, label, shadow = TRUE))
}


edgesHelper <- function(treeeNode){
  terminalFlag <- is.null(treeeNode$children)
  if(!terminalFlag){
    return(data.frame(from = treeeNode$currentIndex, to = treeeNode$children))
  }
}

transIdForBiPlot <- function(treeeList){
  # Keep the tree structure, but change the indexing system
  res <- (seq_along(treeeList) == 1) + 0
  for(i in seq_along(treeeList)){
    treeeNode <- treeeList[[i]]
    terminalFlag <- is.null(treeeNode$children)
    # Index for binary tree
    if(!terminalFlag){
      currentIdx <- res[i]
      res[treeeNode$children] <- 2 * currentIdx + c(0,1)
    }
    # future: index for multi-split tree
  }
  return(res)
}
