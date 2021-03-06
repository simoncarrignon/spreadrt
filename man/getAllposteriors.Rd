% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ABC_analyse.R
\name{getAllposteriors}
\alias{getAllposteriors}
\title{Return posterior given a priori distribution, score and a number of particle to accept}
\usage{
getAllposteriors(allscores, allparameters.dataframe, n = 500)
}
\arguments{
\item{allscores}{a list of list of score}

\item{allparameters.dataframe}{a dataframe with all parameters}

\item{n}{an integer}
}
\value{
a list of list, with for each score and each metrics a dataframe with subset of the parameter for with score of the simulation is in top `n`
}
\description{
Return posterior given a priori distribution, score and a number of particle to accept
}
