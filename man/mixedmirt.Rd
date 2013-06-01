\name{mixedmirt}
\alias{anova,MixedClass-method}
\alias{coef,MixedClass-method}
\alias{mixedmirt}
\alias{summary,MixedClass-method}
\title{Mixed effects modeling for MIRT models}
\usage{
  mixedmirt(data, covdata = NULL, model, fixed = ~1,
    random = NULL, itemtype = "Rasch", itemdesign = NULL,
    constrain = NULL, pars = NULL, ...)
}
\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that
  consists of numerically ordered data, with missing data
  coded as \code{NA}}

  \item{covdata}{a \code{data.frame} that consists of the
  \code{nrow(data)} by \code{K} 'person level' fixed and
  random predictors}

  \item{model}{an object returned from
  \code{confmirt.model()} declaring how the factor model is
  to be estimated. See \code{\link{confmirt.model}} for
  more details}

  \item{fixed}{a standard R formula for specifying the
  fixed effect predictors from \code{covdata} and
  \code{itemdesign}. To estimate the intercepts for each
  item the keyword \code{items} is reserved and
  automatically added to the \code{itemdesign} input}

  \item{random}{a formula similar to the \code{nlme} random
  variable specifications for declaring the random slope
  and intercept predictors. Not currently available, but
  will be available some time in the future}

  \item{itemtype}{same as itemtype in \code{\link{mirt}}}

  \item{itemdesign}{a \code{data.frame} object used to
  create a design matrix for the items, where each
  \code{nrow(itemdesign) == nitems} and the number of
  columns is equal to the number of fixed effect predictors
  (i.e., item intercepts)}

  \item{constrain}{a list indicating parameter equality
  constrains. See \code{\link{mirt}} for more detail}

  \item{pars}{used for parameter starting values. See
  \code{\link{mirt}} for more detail}

  \item{...}{additinonal arguments to be passed to the
  MH-RM estimation engine. See \code{\link{confmirt}} for
  more detail}
}
\description{
  \code{mixedmirt} fits MIRT models using FIML estimation
  to dichotomous and polytomous IRT models conditional on
  fixed and random effect of person and item level
  covariates. The method uses the MH-RM algorithm
  exclusively. The D scaling parameter is automatically
  fixed to 1 so that all coefficients can be interpreted on
  the exponential metric.
}
\examples{
\dontrun{

#make some data
set.seed(1234)
N <- 750
a <- matrix(rlnorm(10,.2,.5),10,1)
d <- matrix(rnorm(10), 10)
Theta <- matrix(sort(rnorm(N)))
pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
group <- factor(rep(c('G1','G2','G3'), each = N/3))
data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
covdata <- data.frame(group, pseudoIQ)
#use cl for parallel computing
library(parallel)
cl <- makeCluster(detectCores())

#specify IRT model
model <- confmirt.model('Theta = 1-10')

#model with no person predictors
mod0 <- mirt(data, model, itemtype = 'Rasch', D = 1)

#group as a fixed effect predictor (aka, uniform dif)
mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, cl=cl)
anova(mod0, mod1)
summary(mod1)
coef(mod1)

#same model as above in lme4
wide <- data.frame(id=1:nrow(data),data,covdata)
long <- reshape2::melt(wide, id.vars = c('id', 'group', 'pseudoIQ'))
library(lme4)
lmod0 <- lmer(value ~ 0 + variable + (1|id), long, family = binomial)
lmod1 <- lmer(value ~ 0 + group + variable + (1|id), long, family = binomial)
anova(lmod0, lmod1)

#model using 2PL items instead of Rasch
mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, itemtype = '2PL', cl=cl)
anova(mod1, mod1b) #much better with 2PL models using all criteria (as expected, given simdata pars)

#continuous predictor and interaction model with group in Rasch model
mod2 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group + pseudoIQ, cl=cl)
mod3 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group * pseudoIQ, cl=cl)
summary(mod2)
anova(mod1b, mod2)
anova(mod2, mod3)

###################################################3
##LLTM, and 2PL version of LLTM
data(SAT12)
data <- key2binary(SAT12,
                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
model <- confmirt.model('Theta = 1-32')

itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))

LLTM <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemdesign = itemdesign, cl=cl)
coef(LLTM)
wald(LLTM)
L <- matrix(c(-1, 1), 1)
wald(LLTM, L) #first half different from second

#compare to items with estimated slopes (2PL)
twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL', itemdesign = itemdesign, cl=cl)
coef(twoPL)
wald(twoPL)
L <- matrix(0, 1, 34)
L[1, 1] <- 1
L[1, 2] <- -1
wald(twoPL, L) #n.s.
anova(twoPL, LLTM)
#twoPL model better than LLTM, and don't draw the incorrect conclusion that the first
#    half of the test is any easier/harder than the last

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}

