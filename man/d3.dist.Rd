\name{d3.dist}
\alias{d3.dist}
\title{
D3-style dissimilarity for a single tree
}
\description{
Compute the set of pairwise dissimilarities across all observations in a tree.
Each dissimilarity measures the extent to which observations are "far apart"
in the tree: the dissimilarity is 0 if the pair land in the same leaf, 1 if
they land on leaves that have only the root as common ancestors, and otherwise
something intermediate.
}
\usage{
d3.dist(mytree, return.pd = FALSE)
}
\arguments{
  \item{mytree}{
Output from "tree"
}
  \item{return.pd}{
If TRUE return the matrix of pairwise distances among leaves. Useful
for debugging. Default FALSE.
}
}
\details{
Two observations have
 distance 0 if they fall in the same leaf; otherwise, the distance measures
 the ratio of the deviance of a tree trimmed so that they do fall in the same
 leaf to the deviance of the original tree.
}
\value{
Item of class "dist" giving inter-point distances.
}
\author{
Sam Buttrey
}

\seealso{
\code{\link{treeClust}}
}
