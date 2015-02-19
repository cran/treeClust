\name{cramer}
\alias{cramer}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute Cramer's V for a two-way table
}
\description{
This function computes the value of Cramer's V for a two-way table.
}
\usage{
cramer(tbl)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{tbl}{
Two-way table, or matrix, of counts.
}
}
\details{
If X^2 is the usual chi-squared measure of association in a two-way table,
Cramer's V is sqrt (X^2 / (n * (k-1))), where n is the total number of 
observations in the table, and k is min (nrow(table), ncol(table)).
}
\value{
Numeric value of Cramer's V, with name "X-squared".
}
\references{
Agresti, "Categorical Data Analysis," p. 75, where V^2 is used.
}
\author{
Sam Buttrey
}
