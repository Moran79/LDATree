% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict.Treee}
\alias{predict.Treee}
\title{Predictions From a Fitted Treee Object}
\usage{
\method{predict}{Treee}(object, newdata, type = c("response", "prob", "all"), ...)
}
\arguments{
\item{object}{A fitted model object of class \code{Treee}, typically the result of
the \code{\link[=Treee]{Treee()}} function.}

\item{newdata}{A data frame containing the predictor variables. Missing
values are allowed and will be handled according to the fitted tree's
method for handling missing data.}

\item{type}{A character string specifying the type of prediction to return.
Options are:
\itemize{
\item \code{'response'}: returns the predicted class for each observation (default).
\item \code{'prob'}: returns a data frame of posterior probabilities for each class.
\item \code{'all'}: returns a data frame containing predicted classes, posterior probabilities, and the predicted node indices.
}}

\item{...}{Additional arguments passed to or from other methods.}
}
\value{
Depending on the value of \code{type}, the function returns:
\itemize{
\item If \code{type = 'response'}: A character vector of predicted class labels.
\item If \code{type = 'prob'}: A data frame of posterior probabilities, where each class has its own column.
\item If \code{type = 'all'}: A data frame containing predicted class labels, posterior probabilities, and the predicted node indices.
}

Note: For factor predictors, if a level not present in the training data is
found in \code{newdata}, it will be treated as missing and handled according to
the \code{missingMethod} specified in the fitted tree.
}
\description{
Generate predictions on new data using a fitted \code{Treee} model.
}
\examples{
fit <- Treee(datX = iris[, -5], response = iris[, 5], verbose = FALSE)
head(predict(fit, iris)) # Predicted classes
head(predict(fit, iris[, -5], type = "prob")) # Posterior probabilities
head(predict(fit, iris[, -5], type = "all")) # Full details
}
