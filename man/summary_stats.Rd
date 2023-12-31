% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_stats.R
\name{summary_stats}
\alias{summary_stats}
\title{calculate summary stats}
\usage{
summary_stats(gin, meta, dataType, popPairwise = FALSE)
}
\arguments{
\item{gin}{is a gtypes object with genetic data}

\item{meta}{is a dataframe with metadata}

\item{dataType}{is a character that defines the type of genetic
data}

\item{popPairwise}{calculate pairwise stats among every population (default FALSE)}
}
\value{
vector of summary statistics
}
\description{
this function calls different summary stats functions
    based on the value of dataType (possibilities include
    "sequence", "snp", "microsatellite")
}
