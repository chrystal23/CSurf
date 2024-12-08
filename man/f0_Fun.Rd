% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/I-SLAMM.R
\name{f0_Fun}
\alias{f0_Fun}
\title{f0_Fun}
\usage{
f0_Fun(X, y)
}
\arguments{
\item{X}{A matrix (n x p) of the covariates. The covariates are also used as potential thresholding variables for subgroup identification.}

\item{y}{A vector (length n) of the true response.}
}
\value{
A list containing the following fields: 
\item{f0}{A vector (length n) of the generated threshold function initialization.}
\item{omeg.hat}{The estimated change plane parameters.}
}
\description{
The function to generate threshold function initialization using change plane model.
}
