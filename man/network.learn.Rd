% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/network_learning.R
\name{network.learn}
\alias{network.learn}
\title{Network learning}
\usage{
network.learn(data_observe,x_cov,X_big_int,X_big_int_cov,
agent,degree=3,len.knots=3,cv=TRUE,nfolds=20,alpha=1)
}
\arguments{
\item{data_observe}{The gene expression matrix.}

\item{x_cov}{The vector including covariate values (eg., smoking: 0 and 1).}

\item{X_big_int}{The matrix based on varying intercepts with the column
containing ones. It can be extracted from the output of function `base.construct`.}

\item{X_big_int_cov}{(X_big_int), The matrix based on varying covariates with the column
containing ones. It can be extracted from the output of function `base.construct`.}

\item{agent}{The imputed disease risk.}

\item{degree}{The degree in the B-spline base. The default is 3}

\item{len.knots}{The number of knots. The default is 3}

\item{cv}{Use cross-validation or not. The default is TRUE. If cv=FALSE,
then BIC method will be applied.}

\item{nfolds}{The number of folds used in cross-validation.
Only workable when cv=TRUE.}

\item{alpha}{0 or 1 values: 1: LASSO; 0: ridge regression.}
}
\value{
The list of network properties: self-node size for baseline, self-node size for covariate effect,
overall self-node size, gene interaction effects from the baseline, gene interaction effects from the covariate,
and overall gene interaction effects
}
\description{
This function can be used to learn the networks related to disease prograssion.
}
