% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_meta.R
\name{is.meta}
\alias{is.meta}
\title{Test if meta has the correct columns}
\usage{
is.meta(meta, required = c("pop", "longpop", "source", "reg"))
}
\arguments{
\item{meta}{meta data}

\item{required}{a vector of names of columns required in meta}
}
\value{
boolean T if all required cols are present in meta
}
\description{
Test if meta has the correct columns
}
