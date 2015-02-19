treeClust.control <- function (return.trees = FALSE, return.mat = TRUE, 
   return.dists = FALSE, return.newdata = FALSE, cluster.only = FALSE, 
       serule = 0, ...) 
{
    list(return.trees = return.trees, return.mat = return.mat, 
        return.dists = return.dists, cluster.only = cluster.only, 
        return.newdata = return.newdata, serule = serule, ...)
}
