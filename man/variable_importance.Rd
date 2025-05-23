% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/variable_importance.R
\name{variable_importance}
\alias{variable_importance}
\title{Get the variable importance of a variable based on the Meta Heuristic run.}
\usage{
variable_importance(
  results,
  column_name,
  col_binary = "Position",
  var_names = NULL
)
}
\arguments{
\item{results}{data.table with the results of the meta heuristic run.
The results should be in line with the transformation of \code{result_2_dt}}

\item{column_name}{column name of \code{results} with the results of the Meta heuristic run}

\item{col_binary}{column name of \code{results} where the binary coding is of the results}

\item{var_names}{vector of labels for the binary coding of the positions.
Length of \code{var_names} should be the same as the binary coding}
}
\value{
list, with names if \code{var_names} is given, with order as the importance of the variable.
}
\description{
Get the variable importance of a variable based on the Meta Heuristic run.
}
