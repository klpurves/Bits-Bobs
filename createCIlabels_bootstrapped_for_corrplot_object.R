## Take a correlation matrix, bootstrap CIs and turn those into text matrixes to overlay on a corrplot object.
#KL Purves

library(corrplot)
library(psych)


M <- cor(mtcars[1:5])

res1 <- cor.ci(M)

lowci <- cor.plot.upperLowerCi(res1)
lowci[upper.tri(lowci)] = t(lowci)[upper.tri(lowci)]

upci <- cor.plot.upperLowerCi(res1)
upci[lower.tri(upci)] = t(upci)[lower.tri(lowci)]


conf <- paste0("[", format(lowci, digits=2), ":", format(upci, digits=2), "]")

xs <- row(lowci)
ys <- (ncol(lowci)+1) - col(lowci)


xs[lower.tri(xs,diag=TRUE)] <- 0 
ys[lower.tri(ys,diag=TRUE)] <- 0 


corrplot(cor(mtcars[1:5]), method="number")

text(xs, ys, conf, pos=1, cex=1)


## Turn a vestor into a upper half matrix (row ise)

length(conf)
mat <- matrix(nrow=length(conf),ncol=length(conf),0)
mat[upper.tri(mat,diag=TRUE)] <- conf





