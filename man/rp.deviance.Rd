\name{rp.deviance}
\alias{rp.deviance}
\title{
Compute deviance within nodes of classification trees
}
\description{
An rpart regression tree carries the deviance around (in the frame$dev
 element). This function computes the deviance for classification trees.
}
\usage{
rp.deviance(x, ...)
}
\arguments{
  \item{x}{
An object of class rpart
}
  \item{\dots}{
Other arguments, currently unused.
}
}
\details{
For a vector of leaf counts n whose sum is N, the deviance is (-2) times 
 the sum of n log (n/N), taking 0 log 0 as 0.
 
}
\value{
Vector of deviances for every row in the tree's frame.
}
\author{
Sam Buttrey
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
rpart
}
