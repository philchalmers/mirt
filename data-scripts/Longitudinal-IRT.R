library("mirt")

numberItems <- 20
numberItemLevels <- 2

a <- matrix(rlnorm(numberItems, .2, .2))
d <- matrix(rnorm(numberItems*numberItemLevels), numberItems)
d <- t(apply(d, 1, sort, decreasing=TRUE))

#Time 1
library(mvtnorm)
set.seed(1)
Theta <- mvtnorm::rmvnorm(n=1000, 0, matrix(1))
t1 <- simdata(a, d, N=1000, itemtype="graded", Theta=Theta)
t2 <- simdata(a, d, N=1000, itemtype="graded", Theta=Theta +.5 + rnorm(nrow(Theta), 0, .3))
t3 <- simdata(a, d, N=1000, itemtype="graded", Theta=Theta + 1 + rnorm(nrow(Theta), 0, .3))

dat <- data.frame(t1, t2, t3)
save(dat, file='longitudinal-IRT.Rdata')
