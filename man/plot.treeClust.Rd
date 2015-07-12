\name{plot.treeClust}
\alias{plot.treeClust}
\title{
Plot treeClust object
}
\description{
Plot a picture of a treeClust object. This picture shows the deviance ratio
on the vertical axis, scaled to have maximum 1, and the tree 
index on the horizontal. Each point is shown by a digit (or digits) giving
the size of the tree.
}
\usage{
\method{plot}{treeClust} (x, extended, \ldots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Object of class treeClust
}
  \item{extended}{
Logical. If TRUE, include all variables, even those whose trees were dropped.
Otherwise only include variables whose trees were kept. Default TRUE.
}
  \item{\dots}{
Other arguments to be passed to the plot function.
}
}
\value{
None. The side effect is that the plot is produced on the current device.
}
