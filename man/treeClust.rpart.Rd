\name{treeClust.rpart}
\alias{treeClust.rpart}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Build an rpart tree as part of treeClust}

\description{
This function builds one tree, as part of a treeClust analysis. It will
not normally be called by users.
}
\usage{
treeClust.rpart(i, dfx, d.num, control, rcontrol)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{i}{
Index of column number (in dfx) of response variable.
}
  \item{dfx}{
Data set used to build tree
}
  \item{d.num}{
Distance number, 1-4, describing measurement for clustering.
}
  \item{control}{
List of controls for treeClust, often output of treeClust.control().
}
  \item{rcontrol}{
List of controls for rpart, often output of rpart.control().
}
}
\details{
It is useful to encapsulate some of the tree-building code so that it can be
used either in a loop or in parallel.
}
\value{
List containing some of these elements (below). Size and DevRatio are always present.
\item{DevRat }{Deviance ratio (decrease in dev. / original dev.) for this tree; always present }
\item{Size }{Size of pruned tree. If no tree is grown, Size is 1.}
\item{tree }{The pruned tree, if needed}
\item{leaf.where }{Vector of leaf membership indices, if Size > 1}
}

\author{
Sam Buttrey
}

\seealso{
\code{\link{treeClust}}
}
