% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abc_model.R
\name{ModelPostProbs}
\alias{ModelPostProbs}
\title{Calculate Posterior Probabilities for Different Models}
\usage{
ModelPostProbs(
  refobj,
  method = c("neuralnet"),
  tol = c(0.1, 0.05, 0.01, 0.005, 0.001, 5e-04),
  modelColumn = "introModel",
  cores = 1,
  ...
)
}
\arguments{
\item{refobj}{A reference object containing the observed data and reference table.}

\item{method}{A character vector indicating the regression method(s). Default is "neuralnet".}

\item{tol}{A numeric vector indicating the tolerance levels for ABC. 
Default is c(0.1,0.05,0.01,0.005,0.001,0.005).}

\item{modelColumn}{A string specifying the model column. Default is "introModel".}

\item{cores}{An integer specifying the number of cores to use for parallel computation. Default is 1.}

\item{...}{Additional arguments passed to the regression methods.}
}
\value{
A list containing the posterior probabilities for different models under 
        various method and tolerance settings.
}
\description{
This function computes the posterior probabilities of different models 
given observed data using Approximate Bayesian Computation (ABC).
}
