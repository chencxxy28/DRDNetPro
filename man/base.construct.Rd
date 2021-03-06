% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_learning.R
\name{base.construct}
\alias{base.construct}
\title{The base generation function}
\usage{
base.construct(data_observe,data_fitted, degree=3,
len.knots=3,data_fitted_cov,agent,x_cov)
}
\arguments{
\item{data_observe}{The gene expression matrix.}

\item{data_fitted}{The matrix containing varying intercepts.}

\item{degree}{The degree in the B-spline base. The default is 3}

\item{len.knots}{The number of knots. The default is 3}

\item{data_fitted_cov}{The matrix containing varying covariate effects}

\item{agent}{The imputed disease risk.}

\item{x_cov}{The vector including covariate values (eg., smoking: 0 and 1).}
}
\value{
The list of four matrices based on varying intercepts (X_big), varying covariate effects (X_big_cov), varying
intercepts with the column containing ones (X_big_int), and varying covariates with the column containing ones (X_big_int_cov).
}
\description{
This function generates the base matrices used in the network learning.
}
