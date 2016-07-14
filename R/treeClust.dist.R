treeClust.dist <- function (x, ...) 
{
#
# This function runs treeClust using its arguments and produces a "dist" object.
#
arg.list <- list (...)
arg.list$dfx <- x
if (any (names (arg.list) == "control"))
    arg.list$control$return.dists <- TRUE
else
    arg.list$control <- treeClust.control (return.dists = TRUE)
do.call ("treeClust", arg.list)$dists
}
