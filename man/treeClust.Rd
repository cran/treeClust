\name{treeClust}
\alias{treeClust}
\title{
Build a tree-based dissimilarity for clustering, and optionally
 perform the clustering
}
\description{
This function uses a set of classification or regression trees to build
 an inter-point dissimilarity in which two points are similar when they
 tend to fall in the same leaves of trees. The user can pass in a clustering
 algorithm and/or ask for the dissimilarities or the set of trees.
}
\usage{
treeClust(dfx, d.num = 1, col.range = 1:ncol(dfx), verbose = F, 
  final.algorithm, k, control = treeClust.control(), rcontrol = rpart.control(), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dfx}{
Input data frame. Columns may be numeric or categorical. Missing values
 are permitted.
}
  \item{d.num}{
Integer: Dissimilarity specifier. When d.num = 1, the dissimilarity between 
 two observations is the proportion of trees where they disagree. 
 With d.num = 2,
 those counts are weighted according to tree quality. In d.num = 3, 
 dissimilarities are variable with trees, reflecting the belief that some pairs
 of leaves are closer together than others. With d.num = 4, those 
 dissimilarities are weighted by tree quality.
}
  \item{col.range}{
Integer: the indices of the columns used. Defaults to all.
}
  \item{verbose}{
If non-zero, print degugging messages to the screen.
}
  \item{final.algorithm}{
Final algorithm, to be used to cluster the computed distances. This may
 be "pam", "agnes", "clara" or "kmeans".
}
  \item{k}{
If final.algorithm is supplied, the number of clusters is required.
}
  \item{control}{
List of the sort produced by \code{\link{treeClust.control}}, giving
 specifications for the fitting routine.
}
  \item{rcontrol}{
List of the sort produced by \code{\link[rpart]{rpart.control}}, giving
 arguments for the rpart routine.
}
  \item{\dots}{
Other arguments, to be passed to the final clustering algorithm if specified.
}
}
\details{
The treeClust approach builds a set of classification or regresion trees,
one for each variable. Trees are pruned, and those that are pruned to the
root are discarded. For each remaining tree, an observation's leaf membership
serves as the starting point for a dissimilarity measurement.
}
\value{
If control$cluster.only is TRUE, a vector of cluster assignments, as produced
 by the final algorthm. Otherwise, a list with these items:
\item{call }{The call that produced the object}
\item{d.num }{d.num, as supplied}
\item{tbl }{Two-column matrix with one row for each tree retained, giving
  size and deviance ratio}
\item{extended.tbl }{Two-column matrix like tbl, but with one row for every
  variable, giving size and deviance ratio (these will be 1 and 0 for variables
  whose trees were discarded}
\item{final.algorithm }{final.algorithm, as supplied}
\item{final.clust }{If final.algorithm is supplied, the output from the final 
 clustering algorithm; otherwise, NULL}
\item{additional.args }{Any additional arguments specified}
\item{tree }{If control$return.trees is TRUE, a list holding all the
 retained trees. This can make the resulting object very large.}
\item{dists }{If control$return.dists is TRUE, an object of class dist
 with the set of pairwise inter-point dissimilarities}
\item{mat }{If control$return.mat is TRUE, a data frame. If final.algorithm is
 "pam" or "agnes" this contains leaf assignment indices. Otherwise this 
  holds a dataset useful as input to k-means or clara. Experimental.}
 
}
\references{
Buttrey and Whitaker, "treeClust: An R Package for Tree-Based Clustering
Dissimilarities," The R Journal, 7/2, 2015.
}
\author{
Sam Buttrey, buttrey@nps.edu
}
\seealso{
\code{\link{treeClust.control}}
}
\examples{
iris.km6 <- treeClust (iris[,-5], d.num = 2, final.algorithm = "kmeans", k=6)
table (iris.km6$final.clust$cluster, iris$Species)
}
