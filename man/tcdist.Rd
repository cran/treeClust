\name{tcdist}
\alias{tcdist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Compute treeClust dissimilarities
}
\description{
Given a treeClust object, or the necessary components, compute all
 pairwise dissimilarities for input to a clustering algorithm
}
\usage{
tcdist(obj, d.num = 1, tbl, mat, trees, verbose=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{obj}{
Object of class treeClust
}
  \item{d.num}{
Method of dissimilarities computation. See "Details".
}
  \item{tbl}{
Two-column of information about trees. Always included in a treeClust
 object, but may be supplied separately. Required if d.num = 2 or 4.
}
  \item{mat}{
Matrix of leaf-membership factors, if not supplied in "obj". 
}
  \item{trees}{
List of trees, if not supplied in obj.
}
  \item{verbose}{
If > 0, print some information useful for debugging.
}
}
\details{
There are four ways to compute inter-point dissimilarities from a 
 treeClust object. If d.num = 1, two points differ by the number of trees
 in which they land in different leaves. "Mat" is required. If d.num = 2,
 the computation for d.num = 1 is used, but each tree gets a different weight.
 "Mat" and "tbl" are required.tbl" are required.

 The computation for d.num = 3 requires that the set of trees be supplied. With
 this approach two observations differ, on a particular tree, according to how 
 far apart they are on that tree. For d.num = 4, both tree and "tbl" 
 are required; this is a weighted version of the d.num = 3 dissimilarity.
}
\value{
Object of class "dist" giving pairwise distances for the original data
 used to build the treeClust object.
}
\author{
Sam Buttrey
}

\seealso{
\code{\link{treeClust}}
}
