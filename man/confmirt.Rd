\name{confmirt}
\alias{anova,confmirt-method}
\alias{coef,confmirt-method}
\alias{confmirt}
\alias{fitted,confmirt-method}
\alias{residuals,confmirt-method}
\alias{summary,confmirt-method}
\title{Confirmatory Full-Information Item Factor Analysis}
\usage{
  confmirt(data, model, itemtype = NULL, guess = 0, upper =
    1, startvalues = NULL, constrain = NULL, freepars =
    NULL, parprior = NULL, verbose = TRUE, calcLL = TRUE,
    draws = 2000, debug = FALSE, rotate = 'varimax', Target
    = NULL, technical = list(), ...)

  \S4method{coef}{confmirt}(object, rotate = '', Target =
    NULL, allpars = FALSE, digits = 3, ...)

  \S4method{summary}{confmirt}(object, digits = 3, ...)

  \S4method{residuals}{confmirt}(object, restype = 'LD',
    digits = 3, printvalue = NULL, ...)

  \S4method{anova}{confmirt}(object, object2, ...)

  \S4method{fitted}{confmirt}(object, digits = 3, ...)
}
\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that
  consists of numerically ordered data}

  \item{model}{an object returned from
  \code{confmirt.model()} declarating how the factor model
  is to be estimated, or a single numeric value indicating
  the number of exploratory factors to estimate. See
  \code{\link{confmirt.model}} for more details}

  \item{guess}{initial (or fixed) values for the
  pseudo-guessing parameter. Can be entered as a single
  value to assign a global guessing parameter or may be
  entered as a numeric vector for each item}

  \item{upper}{initial (or fixed) upper bound parameters
  for 4-PL model. Can be entered as a single value to
  assign a global upper bound parameter or may be entered
  as a numeric vector corresponding to each item}

  \item{printvalue}{a numeric value to be specified when
  using the \code{res='exp'} option. Only prints patterns
  that have standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{verbose}{logical; display iteration history during
  estimation?}

  \item{calcLL}{logical; calculate the log-likelihood via
  Monte Carlo integration?}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood}

  \item{allpars}{logical; print all the item parameters
  instead of just the slopes?}

  \item{restype}{type of residuals to be displayed. Can be
  either \code{'LD'} for a local dependence matrix (Chen &
  Thissen, 1997) or \code{'exp'} for the expected values
  for the frequencies of every response pattern}

  \item{itemtype}{type of items to be modeled, decalred as
  a vector for each item or a single value which will be
  repeated globally. The NULL default assumes that the
  items are ordinal or 2PL, however they may be changed to
  the following: '1PL', '2PL', '3PL', '3PLu', '4PL',
  'graded', 'gpcm', 'nominal', 'mcm', 'PC2PL' and 'PC3PL'
  for the 1 and 2 parameter logistic, 3 parameter logistic
  (lower asymptote and upper), 4 parameter logistic, graded
  response model, generalized partial credit model, nominal
  model, multiple choice model, and 2 and 3PL
  partially-compensatory models, respectively. Note that
  specifying a '1PL' model should be of length 1 (since
  there is only 1 slope parameter estimated). If
  \code{NULL} the defaul assumes that the data follow a
  '2PL' or 'graded' format}

  \item{constrain}{a list of user declared equallity
  constraints. To see how to define the parameters
  correctly use \code{constrain = 'index'} initially to see
  how the parameters are labeled. To constrain parameters
  to be equal create a list with seperate concatenated
  vectors signifying which parameters to contrain. For
  example, to set parameters 1 and 5 equal, and also set
  parameters 2, 6, and 10 equal use \code{constrain =
  list(c(1,5), c(2,6,10))}}

  \item{parprior}{a list of user declared prior item
  probabilities. To see how to define the parameters
  correctly use \code{parprior = 'index'} initially to see
  how the parameters are labeled. Can define either normal
  (normally for slopes and intercepts) or beta (for
  guessing and upper bounds) prior probabilities. Note that
  for upper bounds the value used in the prior is 1 - u so
  that the lower and upper bounds can function the same. To
  specify a prior the form is c('priortype', ...), where
  normal priors are \code{parprior = list(c(parnumber,
  'norm', mean, sd))} and betas are \code{parprior =
  list(c(parnumber, 'beta', alpha, beta))}.}

  \item{freepars}{a list of user declared logical values
  indicating which parameters to estimate. To see how to
  define the parameters correctly use \code{freepars =
  'index'} initially to see how the parameters are labeled.
  These values may be modified and input back into the
  function by using \code{freepars=newfreepars}. Note that
  user input values must match what the default structure
  would have been}

  \item{startvalues}{a list of user declared start values
  for parameters. To see how to define the parameters
  correctly use \code{startvalues = 'index'} initially to
  see what the defaults would noramlly be. These values may
  be modified and input back into the function by using
  \code{startavlues=newstartvalues}. Note that user input
  values must match what the default structure would have
  been}

  \item{debug}{logical; turn on debugging features?}

  \item{object}{an object of class \code{confmirtClass}}

  \item{object2}{an object of class \code{confmirtClass}}

  \item{digits}{the number of significant digits to be
  rounded}

  \item{rotate}{if \code{model} is numeric (indicating an
  exploratory item FA) then this rotation is used. Default
  is \code{'varimax'}}

  \item{Target}{a dummy variable matrix indicing a target
  rotation pattern}

  \item{technical}{list specifying subtle parameters that
  can be adjusted. These values are \describe{
  \item{NCYCLES}{max number of MH-RM cycles; default 2000}
  \item{BURNIN}{number of burn in cycles (stage 1); default
  150} \item{SEMCYCLES}{number of SEM cycles (stage 2);
  default 50} \item{KDRAWS}{number of paralell MH sets to
  be drawn; default 1} \item{TOL}{minimum threshold
  tolerance for convergence of MH-RM, must occur on three
  consecutive occations; default .001} \item{set.seed}{seed
  number used during estimation. Default is 12345}
  \item{guess.prior.n}{a scalar or vector for the weighting
  of the beta priors for guessing parameters (default is
  50, typical ranges are from 2 to 500). If a scalar is
  specified this is used globally, otherwise a numeric
  vector of size \code{ncol(data)} can be used to
  correspond to particualr items (NA values use the
  default)} \item{gain}{a vector of three values specifying
  the numerator, exponent, and subtracted values for the RM
  gain value. Default is \code{c(0.05,0.5,0.004)}} }}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{confmirt} fits a conditional (i.e., confirmatory)
  full-information maximum-likelihood factor analysis model
  to dichotomous and polychotomous data under the item
  response theory paradigm using Cai's (2010)
  Metropolis-Hastings Robbins-Monro algorithm. Fits
  univariate and multivariate Rasch, 1-4PL, graded,
  (generalized) partial credit, nominal, multiple choice,
  and partially-compensatory models, potentially with
  polynomial and product constructed latent traits.
}
\details{
  \code{confmirt} follows a confirmatory and exploratory
  item factor analysis strategy that uses a stochastic
  version of maximum likelihood estimation described by Cai
  (2010a, 2010b). The general equation used for
  multidimensional item response theory in this function is
  in the logistic form with a scaling correction of 1.702.
  This correction is applied to allow comparison to
  mainstream programs such as TESTFACT (2003) and POLYFACT.
  Missing data are treated as 'missing at random' so that
  each response vector is included in the estimation (i.e.,
  full-information). Residuals are computed using the LD
  statistic (Chen & Thissen, 1997) in the lower diagonal of
  the matrix returned by \code{residuals}, and Cramer's V
  above the diagonal. For computing the log-likelihood more
  accurately see \code{\link{logLik}}.

  \code{coef} displays the item parameters with their
  associated standard errors, while use of \code{summary}
  transforms the slopes into a factor loadings metric and
  if the model is exploratory allows for rotating the
  parameters. Also, nested models may be compared by using
  the \code{anova} function, where a Chi-squared difference
  test and AIC/BIC difference values are displayed.
}
\section{Confirmatory IRT}{
  Specification of the confirmatory item factor analysis
  model follows many of the rules in the SEM framework for
  confirmatory factor analysis. The variances of the latent
  factors are automatically fixed to 1 to help facilitate
  model identification. All parameters may be fixed to
  constant values or set equal to other parameters using
  the appropriate declarations.
}

\section{Exploratory IRT}{
  Specifying a number as the second input to confmirt an
  exploratory IRT model is estimatated and can be viewed as
  a stocastic analogue of \code{mirt}, with much of the
  same behaviour and specifications. Rotatation and target
  matrix options will be used in this subroutine and will
  be passed to the returned object for use in generic
  functions such as \code{summary()} and \code{fscores}.
  Again, factor means and variances are fixed to ensure
  proper identification. See \code{\link{mirt}} for more
  details.
}
\examples{
\dontrun{
#Exploratory model estimation, similar to mirt()
data(LSAT7)
fulldata <- expand.table(LSAT7)
(mod1 <- confmirt(fulldata, 1))

#Confirmatory models

#simulate data
a <- matrix(c(
1.5,NA,
0.5,NA,
1.0,NA,
1.0,0.5,
 NA,1.5,
 NA,0.5,
 NA,1.0,
 NA,1.0),ncol=2,byrow=TRUE)

d <- matrix(c(
-1.0,NA,NA,
-1.5,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
3.0,2.0,-0.5,
2.5,1.0,-1,
2.0,0.0,NA,
1.0,NA,NA),ncol=3,byrow=TRUE)

sigma <- diag(2)
sigma[1,2] <- sigma[2,1] <- .4
items <- c(rep('dich',4), rep('graded',3), 'dich')
dataset <- simdata(a,d,2000,items,sigma)

#analyses
#CIFA for 2 factor crossed structure

model.1 <- confmirt.model()
  F1 = 1-4
  F2 = 4-8
  COV = F1*F2


mod1 <- confmirt(dataset,model.1)
coef(mod1)
summary(mod1)
residuals(mod1)

#####
#bifactor
model.3 <- confmirt.model()
  G = 1-8
  F1 = 1-4
  F2 = 5-8


mod3 <- confmirt(dataset,model.3)
coef(mod3)
summary(mod3)
residuals(mod3)
anova(mod1,mod3)

#####
#polynomial/combinations
data(SAT12)
data <- key2binary(SAT12,
                  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))

model.quad <- confmirt.model()
       F1 = 1-32
  (F1*F1) = 1-32


model.combo <- confmirt.model()
       F1 = 1-16
       F2 = 17-32
  (F1*F2) = 1-8


(mod.quad <- confmirt(data, model.quad))
(mod.combo <- confmirt(data, model.combo))
anova(mod.quad, mod.combo)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Cai, L. (2010a). High-Dimensional exploratory item factor
  analysis by a Metropolis-Hastings Robbins-Monro
  algorithm. \emph{Psychometrika, 75}, 33-57.

  Cai, L. (2010b). Metropolis-Hastings Robbins-Monro
  algorithm for confirmatory item factor analysis.
  \emph{Journal of Educational and Behavioral Statistics,
  35}, 307-335.

  Chalmers, R., P. (2012). mirt: A Multidimensional Item
  Response Theory Package for the R Environment.
  \emph{Journal of Statistical Software, 48}(6), 1-29.

  Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S.
  G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for
  Windows: Test Scoring, Item Statistics, and
  Full-information Item Factor Analysis [Computer
  software]. Lincolnwood, IL: Scientific Software
  International.
}
\seealso{
  \code{\link{expand.table}}, \code{\link{key2binary}},
  \code{\link{simdata}}, \code{\link{fscores}},
  \code{\link{confmirt.model}}
}
\keyword{models}

