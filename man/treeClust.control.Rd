\name{treeClust.control}
\alias{treeClust.control}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Parameters describing the output from a treeClust fit
}
\description{
This function produces a list that is used as input to \code{\link{treeClust}}
to determine which items are preserved in the output.
}
\usage{
treeClust.control(return.trees = FALSE, return.mat = TRUE, return.dists = FALSE, 
 return.newdata = FALSE, cluster.only = FALSE, serule = 0, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{return.trees}{
If TRUE, all the trees that go into the object are returned. This can make
the treeClust object very large. Default FALSE.
}
  \item{return.mat}{
If TRUE, return a matrix describing leaf membership. Default TRUE.
}
  \item{return.dists}{
If TRUE, return an object of class 'dissimilarity' giving all pairwise distances
between observations. This can be very large for large datasets. Default FALSE.
}
  \item{cluster.only}{
If TRUE, return only the clustering vector, which names the cluster into which
each observation is places. Default FALSE.
}
  \item{return.newdata}{
If TRUE, return a numeric matrix describing leaf membership and/or inter-point
distance (see "Details"). Default FALSE.
}
  \item{serule}{
Describes how to prune the rpart trees. By default, each tree is pruned to the
minimum error size. With serule > 0, each tree is pruned to the smallest size 
for which the cross-validated error is less than (min error) + (serule * sds).
}
  \item{\dots}{
Other arguments, passed onto the output.
}
}
\details{
The "newdata" item is a numeric matrix that gives inter-point distances whose form depends
on the "d.num" argument to treeClust(). When d.num = 1, each tree contributes a set of 0-1 dummy
variables that serve as leaf membership indicators, and with d.num = 2, each tree's indicators
are multiplied by that tree's "strength." With d.num = 3, a tree with k leaves contributes
k-choose-2 columns, with the distances between distinct rows matching the d3 distances, and
likewise with d.num = 4, a tree with k leaves produced k-choose-2 columns that have been
weighted by tree strength.
}
\value{
list, with all the input arguments and their supplied or default values.
}
\author{
Sam Buttrey, buttrey@nps.edu
}

\seealso{
\code{\link{treeClust}}
}
