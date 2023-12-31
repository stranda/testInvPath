% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abcRF_model.R
\name{predRFmodel}
\alias{predRFmodel}
\title{Predict Using ABC Random Forest Model}
\usage{
predRFmodel(rfFit, refobj, modelColumn = "introModel", cores = 1)
}
\arguments{
\item{rfFit}{A fitted ABC Random Forest model.}

\item{refobj}{A reference object.}

\item{modelColumn}{A string specifying the model column. Default is "introModel".}

\item{cores}{An integer specifying the number of cores to use. Default is 1.}
}
\value{
Predictions from the ABC Random Forest model.
}
\description{
This function makes predictions using an ABC Random Forest model.
}
