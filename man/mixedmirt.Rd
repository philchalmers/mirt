\name{mixedmirt}
\alias{anova,MixedClass-method}
\alias{coef,MixedClass-method}
\alias{mixedmirt}
\alias{summary,MixedClass-method}
\title{Mixed effects modeling for MIRT models}
\usage{
  mixedmirt(data, covdata, model, fixed = ~1,
    random = NULL, itemtype = NULL, fixed.constrain = TRUE,
    ...)
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
  fixed effect predictors from \code{covdata}. By default
  constraints are created so that the fixed effects are
  equal accross items, but this can be disabled using
  \code{fixed.constrain = FALSE}}

  \item{random}{a formula similar to the \code{nlme} random
  variable sepcifications for declaring the random slope
  and intercept predictors. Not currently available, but
  will be soon (hopefully)....}

  \item{itemtype}{same as itemtype in \code{\link{mirt}}}

  \item{fixed.constrain}{logical; constrain the fixed
  effects to be equal accross items?}

  \item{...}{additinonal arguments to be passed to the
  MH-RM estimation engine. See \code{\link{confmirt}} for
  more detail}
}
\description{
  \code{mixedmirt} fits MIRT models using FIML estimation
  to dichotomous and polytomous IRT models conditional on
  fixed and random effect covariates. The method uses the
  MH-RM algorithm exclusively. The D scaling parameter is
  automatically fixed to 1 so that all coefficients can be
  interpreted on the exponential metric.
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

#specify IRT model
model <- confmirt.model()
     Theta = 1-10


#model with no person predictors
mod0 <- mirt(data, model, itemtype = 'Rasch')
#group as a fixed effect predictor (aka, uniform dif)
mod1 <- mixedmirt(data, covdata, model, fixed = ~ group, itemtype = 'Rasch')
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

#model using 2PL items instead of only Rasch
mod1b <- mixedmirt(data, covdata, model, fixed = ~ group)
anova(mod1, mod1b) #much better with 2PL models using all criteria (as expected, given simdata pars)

#global nonuniform dif (group interaction with latent variable)
(dif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta))
anova(mod1b, dif)
#free the interaction terms for each item to detect dif
sv <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, pars = 'values')
constrain <- list(sv$parnum[sv$name == 'groupG2'], sv$parnum[sv$name == 'groupG3']) # main effects
itemdif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, fixed.constrain = FALSE,
     constrain=constrain)
anova(dif, itemdif)

#continuous predictor and interaction model with group
mod2 <- mixedmirt(data, covdata, model, fixed = ~ group + pseudoIQ)
mod3 <- mixedmirt(data, covdata, model, fixed = ~ group * pseudoIQ)
summary(mod2)
anova(mod1b, mod2)
anova(mod2, mod3)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}

