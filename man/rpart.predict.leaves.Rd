\name{rpart.predict.leaves}
\alias{rpart.predict.leaves}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Return the leaf into which observations are predicted to fall
}
\description{
The "where" element of an rpart object gives the leaf into which each
 observation used building the tree falls. This produces the equivalent for
 new data.
}
\usage{
rpart.predict.leaves(rp, newdata, type = "where")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rp}{
Object of class rpart.
}
  \item{newdata}{
New data frame, with the columns used in the rpart model.
}
  \item{type}{
Style of leaf identification: "where" or "leaf"
}
}
\details{
There are two ways to identify the leaf into which an observation falls.
 The way used in the "where" element of an rpart object is to give the row
 number of the leaf within the object's "frame" element. That is the approach
 used here when type = "where". When type = "leaf" the actual leaf number
 is returned. For example, in a tree where node 2 is a terminal node and node 3
 splits into terminal nodes 6 and 7, type = "leaf" will return a vector with
 values 2, 6 and 7. Type = "where" will return a vector with values 2, 4 and 5, 
 since rows 2, 4 and 5 of the tree's "frame" element are leaves.

}
\value{
If type = "where", numeric vector of row numbers describing leaves in the
 tree's "frame" component. If type = "leaf," character vector of leaf numbers.
}
\author{
Sam Buttrey
}
\seealso{
\code{\link[rpart]{rpart}}
}
