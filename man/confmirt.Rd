\name{confmirt}
\alias{anova,confmirt-method}
\alias{coef,confmirt-method}
\alias{confmirt}
\alias{fitted,confmirt-method}
\alias{residuals,confmirt-method}
\alias{summary,confmirt-method}
\title{Confirmatory Full-Information Item Factor Analysis for Mixed Data Formats}
\usage{
  confmirt(data, model, guess = 0, upper = 1, estGuess =
    NULL, estUpper = NULL, ncycles = 2000, burnin = 150,
    SEM.cycles = 50, kdraws = 1, tol = .001, printcycles =
    TRUE, calcLL = TRUE, draws = 2000, returnindex = FALSE,
    debug = FALSE, technical = list(), rotate = 'varimax',
    ...)

  \S4method{coef}{confmirt}(object, SE = TRUE, print.gmeans
    = FALSE, digits = 3, ...)

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

  \item{estGuess}{a logical vector indicating which
  lower-asymptote parameters to be estimated (default is
  null, and therefore is contingent on the values in
  \code{guess}). By default, if any value in \code{guess}
  is greater than 0 then its respective \code{estGuess}
  value is set to \code{TRUE}. Additionally, beta priors
  are automatically imposed for estimated parameters which
  correspond to the input guessing values.}

  \item{estUpper}{same function as \code{estGuess}, but for
  upper bound parameters}

  \item{ncycles}{the maximum number of MH-RM iterations to
  be performed. Default is 2000}

  \item{printvalue}{a numeric value to be specified when
  using the \code{res='exp'} option. Only prints patterns
  that have standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{burnin}{number of burn-in cycles to perform before
  beginning the SEM stage. Default is 150}

  \item{SEM.cycles}{number of stochastic EM cycles to
  perform and average over before beginning the MH-RM
  algorithm. Default is 50}

  \item{kdraws}{number of Metropolis-Hastings imputations
  of the factor scores at each iteration. Default is 1}

  \item{tol}{tolerance that can be reached to terminate the
  model estimation; must be achieved on 3 consecutive
  iterations}

  \item{printcycles}{logical; display iteration history
  during estimation?}

  \item{calcLL}{logical; calculate the log-likelihood via
  Monte Carlo integration?}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood}

  \item{restype}{type of residuals to be displayed. Can be
  either \code{'LD'} for a local dependence matrix (Chen &
  Thissen, 1997) or \code{'exp'} for the expected values
  for the frequencies of every response pattern}

  \item{returnindex}{logical; return the list containing
  the item paramter locations? To be used when specifying
  prior parameter distributions}

  \item{debug}{logical; turn on debugging features?}

  \item{object}{an object of class \code{confmirtClass}}

  \item{object2}{an object of class \code{confmirtClass}}

  \item{SE}{logical; print standard errors?}

  \item{print.gmeans}{logical; print latent factor means?}

  \item{digits}{the number of significant digits to be
  rounded}

  \item{rotate}{if \code{model} is numeric (indicating an
  exploratory item FA) then this rotation is used. Default
  is \code{'varimax'}}

  \item{technical}{list specifying subtle parameters that
  can be adjusted. These values are \describe{
  \item{set.seed}{seed number used during estimation.
  Default is 12345} \item{guess.prior.n}{a scalar or vector
  for the weighting of the beta priors for guessing
  parameters (default is 50, typical ranges are from 2 to
  500). If a scalar is specified this is used globally,
  otherwise a numeric vector of size \code{ncol(data)} can
  be used to correspond to particualr items (NA values use
  the default)} \item{gain}{a vector of three values
  specifying the numerator, exponent, and subtracted values
  for the RM gain value. Default is
  \code{c(0.05,0.5,0.004)}} }}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{confmirt} fits a conditional (i.e., confirmatory)
  full-information maximum-likelihood factor analysis model
  to dichotomous and polychotomous data under the item
  response theory paradigm using Cai's (2010)
  Metropolis-Hastings Robbins-Monro algorithm. If
  requested, lower asymptote parameters are estimated with
  a beta prior included automatically.
}
\details{
  \code{confmirt} follows a confirmatory item factor
  analysis strategy that uses a stochastic version of
  maximum likelihood estimation described by Cai (2010).
  The general equation used for multidimensional item
  response theory in this function is in the logistic form
  with a scaling correction of 1.702. This correction is
  applied to allow comparison to mainstream programs such
  as TESTFACT (2003) and POLYFACT. Missing data are treated
  as 'missing at random' so that each response vector is
  included in the estimation (i.e., full-information).
  Residuals are computed using the LD statistic (Chen &
  Thissen, 1997) in the lower diagonal of the matrix
  returned by \code{residuals}, and Cramer's V above the
  diagonal. For computing the log-likelihood more
  accurately see \code{\link{logLik}}.

  Specification of the confirmatory item factor analysis
  model follows many of the rules in the SEM framework for
  confirmatory factor analysis. The variances of the latent
  factors are automatically fixed to 1 to help facilitate
  model identification. All parameters may be fixed to
  constant values or set equal to other parameters using
  the appropriate declarations. Guessing parameters may be
  specified for dichotomous items and are estimated with
  beta priors automatically, and if a guessing parameter is
  declared for a polychotomous item it is ignored.

  \code{coef} displays the item parameters with their
  associated standard errors, while use of \code{summary}
  transforms the slopes into a factor loadings metric.
  Also, nested models may be compared by using the
  \code{anova} function, where a Chi-squared difference
  test and AIC/BIC difference values are displayed.
}
\examples{
\dontrun{
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
dataset <- simdata(a,d,2000,sigma)

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

#fix first slope at 1.5, and set slopes 7 & 8 to be equal
model.2 <- confmirt.model()
  F1 = 1-4
  F2 = 4-8
  COV = F1*F2
  SLOPE = F1@1 eq 1.5, F2@7 eq F2@8


mod2 <- confmirt(dataset,model.2)
anova(mod2,mod1)

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
#polynomial and combinations
model.linear <- confmirt.model()
      F = 1-8


model.quad <- confmirt.model()
      F = 1-8
  (F*F) = 1-8


model.cube <- confmirt.model()
        F = 1-8
    (F*F) = 1-8
  (F*F*F) = 1-8


model.combo <- confmirt.model()
       F1 = 1-4
       F2 = 5-8
  (F1*F2) = 1-8


mod.linear <- confmirt(dataset, model.linear)
mod.quad <- confmirt(dataset, model.quad)
mod.cube <- confmirt(dataset, model.cube)
mod.combo <- confmirt(dataset, model.combo)

anova(mod.linear,mod.quad)
anova(mod.linear,mod.cube)
anova(mod.linear,mod.combo)
anova(mod.cube,mod.combo)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Cai, L. (2010). Metropolis-Hastings Robbins-Monro
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

