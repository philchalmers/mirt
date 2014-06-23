####################################################
# Exercise 1 - dataset with an incorrect scoring key
####################################################

library(mirt)
require(mvtnorm)
set.seed(1234)

N <- 2500
Theta <- matrix(rnorm(N))
a <- matrix(rlnorm(50, .2, .3))
d <- matrix(rnorm(50, 0, .5))
dcat <- matrix(rnorm(50*4, 0, .5), ncol=4)
data <- simdata(a, cbind(d, dcat), N, itemtype = 'nestlogit', Theta=Theta)

key <- sample(1:5, 50, TRUE)

for(i in 1:length(key)){
    s1 <- data[,i] == key[i]
    s2 <- data[,i] == 1
    data[s1,i] <- 1
    data[s2,i] <- key[i]
}

realkey <- key2binary(data, key)
# mess up key 43
key[42] <- key[42] - 1

save(data, key, file='Exercise_01.Rdata')
brokekey <- key2binary(data, key)

##################################################
# Exercise 2 - multidimensional rating scale model
##################################################

#graded rating scale example

#make some data
set.seed(2)
a <- matrix(c(rep(1, 11), numeric(19), rep(1, 10)), ncol=2)
a <- a[sample(1:20),]
d <- matrix(c(1,0.5,-.5,-1), 20, 4, byrow = TRUE)
c <- seq(-2, 2, length.out=20)
data <- simdata(a, d + c, 2000, itemtype = rep('graded',20), sigma=matrix(c(1,.7,.7,1), 2))
names <- ifelse(as.logical(a[, 1]), 'self', 'others')
colnames(data) <- paste(names, 1:20, sep='_')

save(data, file='Exercise_02.Rdata')

##################################
# Exercise 3 - DIF in 3 parameters 
##################################

library(mirt)
require(mvtnorm)
set.seed(1234)

N <- 1000
Theta <- matrix(rnorm(N))
a <- matrix(rlnorm(20, .2, .3))
d <- matrix(rnorm(20, 0, .5))
dat <- simdata(a, d, N, itemtype = 'dich', Theta=Theta)

Theta <- matrix(rnorm(N, 1, 1.5))
a[2, ] <- .5
a[3,] <- 1.5
d[19,] <- .5
dat2 <- simdata(a, d, N, itemtype = 'dich', Theta=Theta)

group <- rep(c('without_training', 'with_training'), each=N)
data <- rbind(dat, dat2)
save(group, data, file='Exercise_03.Rdata')

#######################################
# Exercise 4 - mixed effects regression
#######################################

set.seed(1234)
N <- 1500
Theta <- rmvnorm(N, c(0,0), matrix(c(1,.6,.6,1), 2))
a <- matrix(c(rep(1, 7), numeric(15), rep(1, 8)), 15, 2)
d <- matrix(rnorm(15), 15)
data <- simdata(a, d, N, itemtype = 'dich', Theta=Theta)

midterm_Math <- round(scale(Theta[,1] * 5 + 100  + rnorm(N, 0 , 5))*10 + 75)
midterm_Math[midterm_Math > 100] <- 100
midterm_Science <- round(scale(Theta[,2] * 5 + 100  + rnorm(N, 0 , 5))*15 + 65)
midterm_Science[midterm_Science > 100] <- 100
tmp <- rowSums(Theta)
Gender <- ifelse(tmp < -1, 'Male', 'Female')
ind <- tmp > -1 & tmp < 1
Gender[ind] <- sample(c('Male', 'Female'), size=sum(ind), replace=TRUE, prob=c(.4,.6))
school <- sort(sample(paste0('school_', 1:50), N, replace = TRUE))
covdata <- data.frame(Gender=Gender, midterm_Math=midterm_Math, midterm_Science=midterm_Science, 
                      school=school)

save(covdata, data, file='Exercise_04.Rdata')