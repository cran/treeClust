\name{leaf.numbers}
\alias{leaf.numbers}
\title{
Convert "where" entry of tree frame into leaf numbers
}
\description{
The "where" entry of a tree object denotes leaves by row numbers in the
 "frame" object. This converts those to actual leaf numbers.
}
\usage{
leaf.numbers(tree)
}
\arguments{
  \item{tree}{
Item of class "tree".
}
}
\value{
Vector, the same length as tree$where, giving leaf numbers.
}
\author{
Sam Buttrey
}
