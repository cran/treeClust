plot.treeClust <- function (x, ...) 
{
#
# Plot method for objects of class treeClust.
#

dev.rat <- x$tbl[,"DevRat"] / max (x$tbl[,"DevRat"])
plot (dev.rat, type = "n", ylim = c(0, 1.05),
 xlim = c(1,length (dev.rat)), xlab = "Variable Number", ylab = "Scaled Deviance Ratio",
  yaxs = "i", ...)

text (1:length (dev.rat), dev.rat, x$tbl[,"Size"])

}

