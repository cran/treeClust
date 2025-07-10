\name{treeClust.dist}
\alias{treeClust.dist}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Built treeClust distance
}
\description{
This function uses treeClust to build a distance. It is intended to act analagously
to \code{\link[cluster]{daisy}} and \code{\link{dist}}.
}
\usage{
treeClust.dist(x, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Data set from which to compute distances via \code{treeClust}.
}
  \item{\dots}{
Other argments to be passed to \code{treeClust}.
}
}
\details{
The \code{treeClust} function's first argument is named dfx. This calls the same code, but by naming
the first argument \code{x} it allows users to employ this function interchangeably with \code{dist} 
and \code{daisy}, which expect arguments named \code{x}. This function also sets the \code{return.dists}
flag and extract the distance object so that that is the only thing returned.
}
\value{
An object of class \code{dissimilarity}. 
}

\author{
Sam Buttrey
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{treeClust}}
}
