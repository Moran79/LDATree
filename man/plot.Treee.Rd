% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot.Treee}
\alias{plot.Treee}
\title{Plot a Treee object}
\usage{
\method{plot}{Treee}(tree, datX, response, node = -1, ...)
}
\arguments{
\item{node}{the node index that you are interested in. By default, it is set
to \code{-1} and the overall tree structure is drawn.}

\item{...}{further arguments passed to or from other methods.}

\item{x}{a fitted model object of class \code{Treee}, which is assumed to be the
result of the \code{\link[=Treee]{Treee()}} function.}

\item{data}{the original data you used to fit the \code{Treee} object if you want
the individual plot for each node. Otherwise, you can leave this parameter
blank if you only need the overall tree structure diagram.}
}
\value{
For overall tree structure (\code{node = -1}), A figure of class
\code{visNetwork} is drawn. Otherwise, a figure of class \code{ggplot} is drawn.
}
\description{
Provide a diagram of the whole tree structure or a scatter/density plot for a
specific tree node.
}
\section{Overall tree structure}{


A full tree diagram (via the R package \link{visNetwork}) is shown if \code{node} is
not provided (default is \code{-1}). The color shows the most common (plurality)
class inside each node. The size of each terminal node is based on its
relative sample size. Under every node, you see the plurality class, the
fraction of the correctly predicted training sample vs. the node's sample
size, and the node index, respectively. When you click on the node, an
information panel with more details will appear.
}

\section{Individual plot for each node}{


The node index and the original training data are required to return a more
detailed plot within a specific node. The density plot will be provided
when only two levels are left for the response variable in a node (like in
a binary classification problem). Samples are projected down to their first
linear discriminant scores (LD1). A scatter plot will be provided if a node
contains more than two classes. Samples are projected down to their first
and second linear discriminant scores.
}

\examples{
fit <- Treee(Species~., data = iris)
# plot the overall tree
plot(fit)
# plot a certain node
plot(fit, iris, node = 1)
}
