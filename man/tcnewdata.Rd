\name{tcnewdata}
\alias{tcnewdata}
\title{
Create all-numeric data to mimic the inter-point distances from treeClust
}
\description{
treeClust produces a vector of dissimilarities, but these objects are
large. This function produces a data frame of data whose inter-point
distances are related to the treeClust ones, for use in, for example,
k-means.
}
\usage{
tcnewdata(obj, d.num = 1, tbl, mat, trees)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
Output from a call to \code{\link{treeClust}}.
}
  \item{d.num}{
Integer, 1-4, describing dissimilarity algorithm. See \code{\link{treeClust}}.
}
  \item{tbl}{
Matrix of tree deviances and sizes, if not present in \code{obj}.
}
  \item{mat}{
Matrix of leaf memberships, if not present in \code{obj}.
}
  \item{trees}{
List of trees, if not present in \code{obj} (needed for d.num = 3 or 4),
}
}
\details{
See the paper by Buttrey and Whitaker. These points have attributes that
match the treeClust distances, but the inter-point distances are not
identical. This is experimental.
}
\value{
Numeric matrix of data whose inter-point distances match the d1 distances
computed by treeClust, and which may be useful for d2-d4 as well.
}
\references{
Buttrey and Whitaker, in submission.
}
\author{
Sam Buttrey, buttrey@nps.edu
}

\seealso{
\code{\link{treeClust}}
}
