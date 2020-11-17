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








