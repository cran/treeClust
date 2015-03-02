treeClust <- function(dfx, d.num = 1, col.range = 1:ncol (dfx), verbose = F, 
         final.algorithm, k, control = treeClust.control(), ...)
{
if(!is.data.frame(dfx))
    stop("This function requires a data frame")

leaf.matrix <- as.data.frame(matrix(0., nrow(dfx), ncol(dfx)))
#
# Fail if any columns have names with embedded spaces.
#
dimnames(leaf.matrix) <- dimnames(dfx)
nm <- names(dfx)
if (length (grep (" ", nm) > 0))
    stop ("Some columns have embedded spaces in their names, and that's not good.")
#
# Final.algorithm has to be agnes, pam, clara, or kmeans (or missing). For everything except
# agnes, a k is required.
#
additional.args <- list (...)
if (!missing (final.algorithm)) {
    if (!is.element (final.algorithm, c("agnes", "pam", "clara", "kmeans")))
        stop ("Unrecognized final algorithm")
#
# For pam, clara or k-means, "k" must be present. For agnes, set it to -1 if 
# it's missing.
#
    if (missing (k))
        if (is.element (final.algorithm, c("kmeans", "pam", "clara")))
            stop ("With kmeans, pam or clara, specify the number of clusters 'k'")
        else k <- -1
        if (k == -1 & control$cluster.only == T) stop ("Cluster.only TRUE requires k")

    if (d.num >= 3 && is.element (final.algorithm, c("clara", "kmeans")))
        stop ("D.num (3 or 4) + (kmeans or clara) not yet supported.")
}
#
# Pairwise distance vector, for d3 (computed later for d1, d2 and d4).
#
if (d.num == 3)
    dists <- numeric (nrow(dfx) * (nrow(dfx) - 1)/2)
#
assign ("dfx", dfx, pos=sys.frame(1))
df.name <- deparse(substitute(dfx))
results <- matrix(0., nrow = ncol(dfx), ncol = 2.)
dimnames(results) <- list(dimnames(dfx)[[2.]], c("DevRat", "Size"))
#
# Set up the big list of trees if asked -- except that we will need
# that list no matter what when d.num == 4.
#
if (control$return.trees || d.num == 4)
    big.list.of.trees <- vector("list", ncol(dfx))
#
# Start the big loop. If the response only has one value, skip it.
#
for(i in col.range) {
    if(verbose > 0)
        cat("Creating rpart tree with column", i, "\n")
    if (length (unique (dfx[,i])) == 1) {
        results[i, "DevRat"] <- 0.
        results[i, "Size"] <- 1.
        next
    }
    if (any (is.na (dfx[,i])))
        response.had.NAs <- TRUE
    else
        response.had.NAs <- FALSE
    str <- paste("rpart (", names(dfx)[i], " ~ ., data = dfx)", sep = "")
#
# Build the tree. If it has one leaf, drop this tree.
#
    mytree <- eval(parse(text = str))
    if (nrow (mytree$frame) == 1) {
        results[i, "DevRat"] <- 0.
        results[i, "Size"] <- 1.
        next
    }
#
# Extract the CP table. Compute the cutoff value based on the one-se
# rule (or whatever value serule has). Then find the smallest row whose
# xerror is smaller than that value. If that's row 1, drop this tree.
#
    cptbl <- mytree$cptable
    min.cp.dex <- which (cptbl[,"xerror"] == min(cptbl[,"xerror"]))[1]
    serule.value <-         cptbl[min.cp.dex,"xerror"] + 
           control$serule * cptbl[min.cp.dex,"xstd"]
    best.row <- min(which (cptbl[,"xerror"] <= serule.value))
    if(best.row == 1.) {
        results[i, "DevRat"] <- 0.
        results[i, "Size"] <- 1.
        next
    }    
#
# Prune to the CP value in this row, except so as to avoid rounding error,
# prune to something a little bigger -- say, halfway between the CP in
# this row and the one above.
#
    mytree <- prune.rpart (mytree, 
              cp = (cptbl[best.row, "CP"] + cptbl[best.row - 1,"CP"])/2)
#
# Save the leaf membership values. If this response had NAs, though, we
# will need to generate the full set, first.
#
    if (response.had.NAs) {
        mytree$where.orig <- mytree$where # for debugging
        mytree$where <- rpart.predict.leaves (mytree, dfx, type = "where")
# "where" without missings
        leaf.matrix[,i] <- factor (mytree$where)
    }
    else 
        leaf.matrix[, i] <- factor (mytree$where)
#
# Save the tree if asked -- or if we're using d4.
# 
    if (control$return.trees || d.num == 4)
        big.list.of.trees[[i]] <- mytree
#
# The thing named "dev" really is the deviance for a regression tree,
# but not for a classification tree. Those we have to compute ourselves.
#
    if (is.factor (dfx[,i]))
        devs <- rp.deviance (mytree)
    else
        devs <- mytree$frame$dev
    orig.dev <- devs[1]
    new.dev <- sum (devs[mytree$frame$var == "<leaf>"])

    results[i, "DevRat"] <- (orig.dev - new.dev)/orig.dev
    results[i, "Size"] <- sum (mytree$frame[,"var"] == "<leaf>")
#
# If d.num is 1 or 2, we have everything we need. 
#
    if (d.num <= 2)
        next
#
# For d.num = 3, compute and accumulate distances. For 4, we have to 
# wait until all the trees are assembled.
#
    if (d.num == 3)
        dists <- dists + d3.dist (mytree)

} # end "for"
# -----------------------------
# End of for loop over columns
# ------------------------------

if(!any(results[, "Size"] > 1))
    stop("No tree produced anything! Panic!")
leaf.matrix <- leaf.matrix[, results[, "Size"] > 1., drop=F]
if (control$return.trees || d.num == 4)
    big.list.of.trees <- big.list.of.trees[results[, "Size"] > 1.]
results <- results[results[,"Size"] > 1,, drop=F]
#
# If there's a final algorithm, use it now. Agnes and pam use
# the inter-point distances, plus any additional arguments, so we'll
# compute those. Also compute them if we've asked for them explicitly.
#
if (control$return.dists == TRUE ||
    (!missing (final.algorithm) && is.element (final.algorithm, c("pam", "agnes"))))
{
# If d.num is 1 or 2 we don't need the big.list of trees. For d3, we've 
# already been computing the dists, so we don't need to do anything. For
# d4, we know we have the trees.
#
    if (d.num == 1 || d.num == 2)
        dists <- tcdist (tbl = results, mat = leaf.matrix)
    else if (d.num == 4) {
            dists <- tcdist (tbl = results, mat = leaf.matrix, 
                                 trees = big.list.of.trees, d.num = d.num)
    }
}
if (missing (final.algorithm)) {
    final.algorithm <- "None"
    final.clust <- NULL
} else {
#
# We call "agnes" or "pam" by "do.call", which saves a copy of the dists
# in the call element. That thing is huge and unnecessary, so we remove it.
#
    if (final.algorithm == "agnes") {
        final.clust <- do.call (final.algorithm, list (x = dists, ...))
        final.clust$call$x <- "deleted"
        if (control$cluster.only == TRUE)
            final.clust <- cutree (final.clust, k = k)
    }
    if (final.algorithm == "pam") {
        final.clust <- do.call (final.algorithm, list (x = dists, k = k, ...))  
        final.clust$call$x <- "deleted"
        if (control$cluster.only == TRUE)
            final.clust <- final.clust$clustering
    }
#
# k-means are clara are tricker. We construct the "new" data and run the
# algorithm on that. The new data has p columns for each
# tree with p leaves -- we could make do with one less.
#
    if (is.element (final.algorithm, c("clara", "kmeans"))) {
        newdata <- tcnewdata (tbl = results, mat = leaf.matrix, 
                                 trees = big.list.of.trees, d.num = d.num)

        if (final.algorithm == "kmeans")
            final.clust <- kmeans (x = newdata, centers = k)
        else 
            final.clust <- clara (x = newdata, k = k)
        if (control$cluster.only == TRUE)
                final.clust <- final.clust$cluster
    }
} # end of "final algorithm" stuff

#
# If we were only asked for the clustering, return that.
#
if (control$cluster.only)
    return (final.clust)
#
# Set up return value; add requested stuff.
#
return.val <- list(call = match.call(), d.num = d.num, tbl = results, final.algorithm = final.algorithm, 
                   final.clust = final.clust, additional.args = additional.args)
if (control$return.trees)
    return.val$trees <- big.list.of.trees
if (control$return.dists)
    return.val$dists <- dists
if (control$return.mat)
    if (!missing (final.algorithm) && is.element (final.algorithm, c("clara", "kmeans")))
        return.val$mat <- newdata
    else
        return.val$mat <- leaf.matrix
class(return.val) <- "treeClust"

return(return.val)
}

