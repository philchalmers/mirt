\name{simdata}
\alias{simdata}
\title{Simulate response patterns}
\usage{
  simdata(a, d, N, itemtype, sigma = NULL, mu = NULL,
    guess = 0, upper = 1, nominal = NULL, Theta = NULL,
    D = 1.702)
}
\arguments{
  \item{a}{a matrix of slope parameters. If slopes are to
  be constrained to zero then use \code{NA}. \code{a} may
  also be a similar matrix specifying factor loadings if
  \code{factor.loads = TRUE}}

  \item{d}{a matrix of intercepts. The matrix should have
  as many columns as the item with the largest number of
  categories, and filled empty locations with \code{NA}}

  \item{itemtype}{a character vector of length
  \code{nrow(a)} (or 1, if all the item types are the same)
  specifying the type of items to simulate. Can be
  \code{'dich', 'graded', 'gpcm','nominal', 'nestlogit'},
  or \code{'partcomp'}, for dichotomous, graded,
  generalized partial credit, nominal, nested logit, and
  partially compensatory models. Note that for the gpcm,
  nominal, and nested logit models there should be as many
  parameters as desired categories, however to
  parameterized them for meaningful interpretation the
  first category intercept should equal 0 for these models
  (second column for \code{'nestlogit'}, since first column
  is for the correct item traceline). For nested logit
  models the 'correct' category is always the lowest
  category (i.e., == 1)}

  \item{nominal}{a matrix of specific item category slopes
  for nominal models. Should be the dimensions as the
  intercept specification with one less column, with
  \code{NA} in locations where not applicable. Note that
  during estimation the first slope will be constrained to
  0 and the last will be constrained to the number of
  categories minus 1, so it is best to set these as the
  values for the first and last categories as well}

  \item{N}{sample size}

  \item{guess}{a vector of guessing parameters for each
  item; only applicable for dichotomous items. Must be
  either a scalar value that will affect all of the
  dichotomous items, or a vector with as many values as to
  be simulated items}

  \item{upper}{same as \code{guess}, but for upper bound
  parameters}

  \item{D}{logistic scaling parameter, default 1.702}

  \item{sigma}{a covariance matrix of the underlying
  distribution. Default is the identity matrix}

  \item{mu}{a mean vector of the underlying distribution.
  Default is a vector of zeros}

  \item{Theta}{a user specified matrix of the underlying
  ability parameters, where \code{nrow(Theta) == N} and
  \code{ncol(Theta) == ncol(a)}}
}
\description{
  Simulates response patterns for compensatory and
  noncompensatory MIRT models from multivariate normally
  distributed factor (\eqn{\theta}) scores, or from a user
  input matrix of \eqn{\theta}'s.
}
\details{
  Returns a data matrix simulated from the parameters.
}
\examples{
\dontrun{
###Parameters from Reckase (2009), p. 153

set.seed(1234)

a <- matrix(c(
 .7471, .0250, .1428,
 .4595, .0097, .0692,
 .8613, .0067, .4040,
1.0141, .0080, .0470,
 .5521, .0204, .1482,
1.3547, .0064, .5362,
1.3761, .0861, .4676,
 .8525, .0383, .2574,
1.0113, .0055, .2024,
 .9212, .0119, .3044,
 .0026, .0119, .8036,
 .0008, .1905,1.1945,
 .0575, .0853, .7077,
 .0182, .3307,2.1414,
 .0256, .0478, .8551,
 .0246, .1496, .9348,
 .0262, .2872,1.3561,
 .0038, .2229, .8993,
 .0039, .4720, .7318,
 .0068, .0949, .6416,
 .3073, .9704, .0031,
 .1819, .4980, .0020,
 .4115,1.1136, .2008,
 .1536,1.7251, .0345,
 .1530, .6688, .0020,
 .2890,1.2419, .0220,
 .1341,1.4882, .0050,
 .0524, .4754, .0012,
 .2139, .4612, .0063,
 .1761,1.1200, .0870),30,3,byrow=TRUE)

d <- matrix(c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
  .6431,.0122,.0912,.8082,-.1867,.4533,-1.8398,.4139,
  -.3004,-.1824,.5125,1.1342,.0230,.6172,-.1955,-.3668,
  -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329),ncol=1)

mu <- c(-.4, -.7, .1)
sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)

dataset1 <- simdata(a, d, 2000, itemtype = 'dich')
dataset2 <- simdata(a, d, 2000, itemtype = 'dich', mu = mu, sigma = sigma)

#mod <- confmirt(dataset1, 3)
#coef(mod)

###An example of a mixed item, bifactor loadings pattern with correlated specific factors

a <- matrix(c(
.8,.4,NA,
.4,.4,NA,
.7,.4,NA,
.8,NA,.4,
.4,NA,.4,
.7,NA,.4),ncol=3,byrow=TRUE)

d <- matrix(c(
-1.0,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
0.0,-1.0,1.5,  #the first 0 here is the recommended constraint for nominal
0.0,1.0,-1, #the first 0 here is the recommended constraint for gpcm
2.0,0.0,NA),ncol=3,byrow=TRUE)

nominal <- matrix(NA, nrow(d), ncol(d))
nominal[4, ] <- c(0,1.2,2) #the first 0 and last (ncat - 1) = 2 values are the recommended constraints

sigma <- diag(3)
sigma[2,3] <- sigma[3,2] <- .25
items <- c('dich','dich','dich','nominal','gpcm','graded')

dataset <- simdata(a,d,2000,items,sigma=sigma,nominal=nominal)

#mod <- bfactor(dataset, c(1,1,1,2,2,2), itemtype = c(rep('2PL', 3), 'nominal', 'gpcm','graded'))
#coef(mod)

####Unidimensional nonlinear factor pattern

theta <- rnorm(2000)
Theta <- cbind(theta,theta^2)

a <- matrix(c(
.8,.4,
.4,.4,
.7,.4,
.8,NA,
.4,NA,
.7,NA),ncol=2,byrow=TRUE)
d <- matrix(rnorm(6))
itemtype <- rep('dich',6)

nonlindata <- simdata(a,d,2000,itemtype,Theta=Theta)

#model <- confmirt.model('
#F1 = 1-6
#(F1 * F1) = 1-3')
#mod <- confmirt(nonlindata, model)
#coef(mod)

####2PLNRM model for item 4 (with 4 categories), 2PL otherwise

a <- matrix(rlnorm(4,0,.2))

#first column of item 4 is the intercept for the correct category of 2PL model,
#    otherwise nominal model configuration
d <- matrix(c(
-1.0,NA,NA,NA,
 1.5,NA,NA,NA,
 0.0,NA,NA,NA,
 1, 0.0,-0.5,0.5),ncol=4,byrow=TRUE)

nominal <- matrix(NA, nrow(d), ncol(d))
nominal[4, ] <- c(NA,0,.5,.6)

items <- c(rep('dich',3),'nestlogit')

dataset <- simdata(a,d,2000,items,nominal=nominal)

#mod <- mirt(dataset, 1, itemtype = c('2PL', '2PL', '2PL', '2PLNRM'), key=c(NA,NA,NA,1))
#coef(mod)
#itemplot(mod,4)


   }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Reckase, M. D. (2009). \emph{Multidimensional Item
  Response Theory}. New York: Springer.
}
\keyword{data}

