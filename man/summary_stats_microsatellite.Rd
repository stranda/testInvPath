% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary_stats.R
\name{summary_stats_microsatellite}
\alias{summary_stats_microsatellite}
\title{Summary Statistics for Microsatellites}
\usage{
summary_stats_microsatellite(gin, meta, popPairwise = FALSE)
}
\arguments{
\item{gin}{Genetic input data in strataG gtypes format}

\item{meta}{Metadata related to the `gin` input.}
}
\value{
A vector containing various summary statistics. The output is ordered to ensure comparability 
between original and simulated statistics.
}
\description{
This function computes various summary statistics related to microsatellite data.
}
\details{
The function computes multiple summary statistics such as overall heterozygosity, 
population-specific heterozygosity, pairwise population statistics, and so on.
The method to extract population information from `gin` is also determined within the function based on 
the presence of a specific pattern in the individual names.
}
\seealso{
\code{\link{heterozygosity}}, \code{\link{overallTest}}, \code{\link{pairwiseTest}}, \code{\link{mRatio}}
}
