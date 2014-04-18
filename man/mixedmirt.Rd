\name{mixedmirt}
\alias{mixedmirt}
\title{Mixed effects modeling for MIRT models}
\usage{
mixedmirt(data, covdata = NULL, model, fixed = ~1, random = NULL,
  itemtype = "Rasch", itemdesign = NULL, constrain = NULL, pars = NULL,
  return.design = FALSE, SE = TRUE, draws = 5000, TOL = NULL, ...)
}
\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that
  consists of numerically ordered data, with missing data
  coded as \code{NA}}

  \item{covdata}{a \code{data.frame} that consists of the
  \code{nrow(data)} by \code{K} 'person level' fixed and
  random predictors}

  \item{model}{an object returned from \code{mirt.model()}
  declaring how the factor model is to be estimated. See
  \code{\link{mirt.model}} for more details}

  \item{fixed}{a right sided R formula for specifying the
  fixed effect (aka 'explanatory') predictors from
  \code{covdata} and \code{itemdesign}. To estimate the
  intercepts for each item the keyword \code{items} is
  reserved and automatically added to the \code{itemdesign}
  input.  If any polytomous items are being model the
  \code{items} are argument is not valid since all
  intercept parameters are freely estimated and identified
  with the parameterizations found in \code{\link{mirt}},
  and the first column in the fixed design matrix (commonly
  the intercept or a reference group) is omitted}

  \item{random}{a right sided formula or list of formulas
  containing crossed random effects of the form \code{v1 +
  ... v_n | G}, where \code{G} is the grouping variable and
  \code{v_n} are random numeric predictors within each
  group. If no intercept value is specified then by default
  the correlations between the \code{v}'s and \code{G} are
  estimated, but can be suppressed by including the \code{~
  -1 + ...} constant}

  \item{itemtype}{same as itemtype in \code{\link{mirt}},
  expect does not support the following item types:
  \code{c('PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM',
  '4PLNRM')}}

  \item{itemdesign}{a \code{data.frame} object used to
  create a design matrix for the items, where each
  \code{nrow(itemdesign) == nitems} and the number of
  columns is equal to the number of fixed effect predictors
  (i.e., item intercepts). By default an \code{items}
  variable is reserved for modeling the item intercept
  parameters}

  \item{constrain}{a list indicating parameter equality
  constrains. See \code{\link{mirt}} for more detail}

  \item{pars}{used for parameter starting values. See
  \code{\link{mirt}} for more detail}

  \item{return.design}{logical; return the design matrices
  before they have (potentially) been reassigned?}

  \item{SE}{logical; while standard errors are always
  computed for the MHRM algorithm this option forces the
  number of iterations to be no less than 400 so that the
  information matrix is computed more accurately.  An
  alternative approach would be to drop the \code{TOL}
  criteria}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood for the MH-RM algorithm. Default is
  5000}

  \item{TOL}{convergence threshold for MH-RM; default is
  .001}

  \item{...}{additional arguments to be passed to the MH-RM
  estimation engine. See \code{\link{mirt}} for more
  detail}
}
\description{
\code{mixedmirt} fits MIRT models using FIML estimation to
dichotomous and polytomous IRT models conditional on fixed
and random effect of person and item level covariates. This
can also be understood as 'explanatory IRT' if only fixed
effects are modeled, or multilevel/mixed IRT if random and
fixed effects are included. The method uses the MH-RM
algorithm exclusively. Additionally, computation of the
log-likelihood can be sped up by using parallel estimation
via \code{\link{mirtCluster}}.
}
\details{
For dichotomous response models, \code{mixedmirt} follows
the general form

\deqn{P(x = 1|\theta, \psi) = g + \frac{(u - g)}{1 + exp(-1
* [\mathbf{\theta a} + \mathbf{X \beta} + \mathbf{Z
\gamma}])}}

where X is a design matrix with associated \eqn{\beta}
fixed effect coefficients, and Z is a design matrix with
associated \eqn{\gamma} random effects. For simplicity and
easier interpretation, the unique item intercept values
typically found in \eqn{\mathbf{X \beta}} are extracted and
reassigned within mirt's 'intercept' parameters (e.g.,
\code{'d'}).  To observe how the design matrices are
structured prior to reassignment and estimation pass the
argument \code{return.design = TRUE}.

Polytomous IRT models follow a similar format except the
item intercepts are automatically estimated internally,
rendering the \code{items} argument in the fixed formula
redundant and therefore must be omitted from the
specification. If there are a mixture of dichotomous and
polytomous items the intercepts for the dichotomous models
are also estimated for consistency.

To simulate maximum a posteriori estimates for the random
effects use the \code{\link{randef}} function.
}
\examples{
\dontrun{

#make some data
set.seed(1234)
N <- 750
a <- matrix(rlnorm(10,.3,1),10,1)
d <- matrix(rnorm(10), 10)
Theta <- matrix(sort(rnorm(N)))
pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
pseudoIQ <- (pseudoIQ - mean(pseudoIQ))/10  #rescale variable for numerical stability
group <- factor(rep(c('G1','G2','G3'), each = N/3))
data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
covdata <- data.frame(group, pseudoIQ)
#use parallel computing
mirtCluster()

#specify IRT model
model <- mirt.model('Theta = 1-10')

#model with no person predictors
mod0 <- mirt(data, model, itemtype = 'Rasch')

#group as a fixed effect predictor (aka, uniform dif)
mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items)
anova(mod0, mod1)
summary(mod1)
coef(mod1)

#same model as above in lme4
wide <- data.frame(id=1:nrow(data),data,covdata)
long <- reshape2::melt(wide, id.vars = c('id', 'group', 'pseudoIQ'))
library(lme4)
lmod0 <- glmer(value ~ 0 + variable + (1|id), long, family = binomial)
lmod1 <- glmer(value ~ 0 + group + variable + (1|id), long, family = binomial)
anova(lmod0, lmod1)

#model using 2PL items instead of Rasch
mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items, itemtype = '2PL')
anova(mod1, mod1b) #better with 2PL models using all criteria (as expected, given simdata pars)

#continuous predictor with group
mod2 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items + pseudoIQ)
summary(mod2)
anova(mod1b, mod2)

#view fixed design matrix with and without unique item level intercepts
withint <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, return.design = TRUE)
withoutint <- mixedmirt(data, covdata, model, fixed = ~ 0 + group, return.design = TRUE)

#notice that in result above, the intercepts 'items1 to items 10' were reassigned to 'd'
head(withint$X)
tail(withint$X)
head(withoutint$X) #no intercepts design here to be reassigned into item intercepts
tail(withoutint$X)

###################################################
### random effects
#make the number of groups much larger
covdata$group <- factor(rep(paste0('G',1:50), each = N/50))

#random groups
rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group)
summary(rmod1)
coef(rmod1)

#random groups and random items
rmod2 <- mixedmirt(data, covdata, 1, random = list(~ 1|group, ~ 1|items))
summary(rmod2)
eff <- randef(rmod2) #estimate random effects

#random slopes with fixed intercepts (suppressed correlation)
rmod3 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ -1 + pseudoIQ|group)
summary(rmod3)
(eff <- randef(rmod3))

###################################################
##LLTM, and 2PL version of LLTM
data(SAT12)
data <- key2binary(SAT12,
                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
model <- mirt.model('Theta = 1-32')

itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))

#notice that the 'fixed = ~ ... + items' argument is omitted
LLTM <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemdesign = itemdesign,
   SE = TRUE) # SE argument ensures that the information matrix is computed accurately
summary(LLTM)
coef(LLTM)
wald(LLTM)
L <- c(-1, 1, 0)
wald(LLTM, L) #first half not different from second

#compare to items with estimated slopes (2PL)
twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL',
                   itemdesign = itemdesign)
#twoPL not mixing too well (AR should be between .2 and .5), decrease MHcand
twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL',
                  itemdesign = itemdesign, technical = list(MHcand = 0.8))
anova(twoPL, LLTM) #much better fit
summary(twoPL)
coef(twoPL)

wald(twoPL)
L <- matrix(0, 1, 34)
L[1, 1] <- 1
L[1, 2] <- -1
wald(twoPL, L) #n.s.

##LLTM with item error term
LLTMwithError <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, random = ~ 1|items,
    itemdesign = itemdesign)
summary(LLTMwithError)
#large item level variance after itemorder is regressed; not a great predictor of item difficulty
coef(LLTMwithError)

###################################################
### Polytomous example

#make an arbitrary group difference
covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))

#partial credit model
mod <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group)
coef(mod)

#gpcm to estimate slopes
mod2 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
                 itemtype = 'gpcm')
summary(mod2)
anova(mod, mod2)

#graded model
mod3 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
                 itemtype = 'graded')
coef(mod3)

}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
\code{\link{anova-method}}, \code{\link{coef-method}},
\code{\link{summary-method}} \code{\link{randef}},
\code{\link{calcLogLik}}, \code{\link{mirtCluster}},
\code{\link{mirt}}
}

