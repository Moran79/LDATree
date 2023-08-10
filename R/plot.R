#' Title
#'
#' @param treeeOutput
#' @param data
#' @param node
#'
#' @return
#' @export
#'
#' @examples
plot.Treee <- function(treeeOutput, data, node = -1){
  if(node>0){
    if(missing(data)) stop("Please input the orginal training data for nodewise LDA plots")
    if(treeeOutput$treee[[node]]$nodeModel == "mode") return(paste("Every observation in this node is predicted to be", treeeOutput$treee[[node]]$nodePredict))
    # Get the data ready, impute the NAs (if any)
    dataProcessed <- extractXnResponse(treeeOutput$formula, data)
    newX <- getDataInShape(data = dataProcessed$x[treeeOutput$treee[[node]]$idxRow,], missingReference = treeeOutput$treee[[node]]$misReference)
    colorIdx <- match(names(treeeOutput$treee[[node]]$proportions), levels(dataProcessed$response))

    plotLDA2d(ldaModel = treeeOutput$treee[[node]]$nodePredict,
              data = cbind.data.frame(response = dataProcessed$response[treeeOutput$treee[[node]]$idxRow], newX),
              node = node,
              colorManual = hue_pal()(nlevels(dataProcessed$response))[colorIdx])
  }else{ # default overall plot
    plot(treeeOutput$treee)
  }
}

#' Title
#'
#' @param treeeList
#'
#' @return
#' @export
#'
#' @examples
plot.SingleTreee <- function(treeeList){
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
  # line4 = paste('</br>The proportion of', paste(names(treeeNode$proportions), collapse = ', '),'are',
  #               paste(sprintf("%.1f%%", treeeNode$proportions / length(treeeNode$idxRow) * 100), collapse = ', '))
  line4 = paste('</br>', length(treeeNode$idxRow) - treeeNode$currentLoss, 'of them are correctly classified')
  line5 = paste('</br>The accuracy is ', round(treeeNode$accuracy,3))
  line5.5 = paste('</br>Plurality class (', round(max(treeeNode$proportions) / sum(treeeNode$proportions),4)*100, '%) is ', names(sort(treeeNode$proportions, decreasing = TRUE))[1], sep = "")
  line6 = paste('</br>The model in this node is ', treeeNode$nodeModel)
  return(paste(line1,line2,line3,line4,line5,line5.5,line6))
}


nodesHelper <- function(treeeNode, idTransVec){
  terminalFlag <- is.null(treeeNode$children)
  id = treeeNode$currentIndex
  title = infoClickSingle(treeeNode = treeeNode, idTransVec = idTransVec) # Show when you click
  value = ifelse(terminalFlag, log(length(treeeNode$idxRow)), 2) # node size
  level = treeeNode$currentLevel
  group = names(sort(treeeNode$proportions, decreasing = TRUE))[1]
  label = paste(group, # paste(treeeNode$proportions, collapse = ' / '),
                paste(length(treeeNode$idxRow) - treeeNode$currentLoss, length(treeeNode$idxRow), sep = ' / '),
                paste('Node', idTransVec[id]),sep = "\n")
                # paste('alpha:', treeeNode$alpha),
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
  if(dim(ldaModel$scaling)[2] == 1){
    # Only one LD is available, draw the histogram
    datPlot <- cbind.data.frame(response = data$response, LD1 = getLDscores(modelLDA = ldaModel, data = data, nScores = 1))
    p <- ggplot(data = datPlot)+
      geom_density(aes(x = LD1, fill = response), alpha = 0.7)+
      scale_fill_manual(values = colorManual)+
      theme_bw()+
      labs(title = "Density plot of LD1", subtitle = paste("Node:",node))
  }else{
    LDscores <- getLDscores(modelLDA = ldaModel, data = data, nScores = 2)
    datPlot <- cbind.data.frame(response = data$response, LDscores)
    p <- ggplot(data = datPlot)+
      geom_point(aes(x = LD1, y = LD2, color = response), size = 3, alpha = 0.7)+
      scale_color_manual(values = colorManual)+
      theme_bw()+
      labs(title = "Scatter plot by first two LDscores", subtitle = paste("Node:",node))
  }
  return(p)
}



