\name{make.leaf.paths}
\alias{make.leaf.paths}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make matrix of leaf paths
}
\description{
It is helpful to know the parent nodes for each tree node. This
function creates a matrix with that information.
}
\usage{
make.leaf.paths(up.to = 2047)
}
\arguments{
  \item{up.to}{
Number of rows for which to compute leaf.paths.
}
}
\details{
The ith row of the resulting matrix lists all the leaves, including i,
 that would be traversed from the root to leaf i. Unneeded columns have zeros.
}
\value{
Numeric matrix with "up.to" rows. If 2^j <= up.to < 2^(j+1), j columns.
}
