\name{polymirt}
\alias{anova,polymirt-method}
\alias{coef,polymirt-method}
\alias{fitted,polymirt-method}
\alias{plot,polymirt-method}
\alias{polymirt}
\alias{residuals,polymirt-method}
\alias{summary,polymirt-method}
\title{Full-Information Item Factor Analysis for Mixed Data Formats}
\usage{
  polymirt(data, nfact, guess = 0, upper = 1, estGuess =
    NULL, estUpper = NULL, prev.cor = NULL, rotate =
    'varimax', verbose = TRUE, calcLL = TRUE, draws = 2000,
    debug = FALSE, technical = list(), ...)

  \S4method{summary}{polymirt}(object, rotate='', suppress
    = 0, digits = 3, print = FALSE, ...)

  \S4method{coef}{polymirt}(object, rotate = '', SE = TRUE,
    digits = 3, ...)

  \S4method{plot}{polymirt}(x, npts = 50, type = 'info',
    rot = list(x = -70, y = 30, z = 10), ...)

  \S4method{residuals}{polymirt}(object, restype = 'LD',
    digits = 3, printvalue = NULL, ...)

  \S4method{anova}{polymirt}(object, object2, ...)

  \S4method{fitted}{polymirt}(object, digits = 3, ...)
}
\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that
  consists of numerically ordered data}

  \item{nfact}{number of factors to be extracted}

  \item{guess}{starting (or fixed) values for the
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
  are automatically imposed for estimated parameters that
  correspond to the input guessing value.}

  \item{estUpper}{same function as \code{estGuess}, but for
  upper bound parameters}

  \item{prev.cor}{use a previously computed correlation
  matrix to be used to estimate starting values the
  estimation. The input could be any correlation matrix,
  but it is advised to use a matrix of polychoric
  correlations.}

  \item{rotate}{type of rotation to perform after the
  initial orthogonal parameters have been extracted by
  using \code{summary}; default is \code{'varimax'}. See
  \code{\link{mirt}} for list of possible rotations. If
  \code{rotate != ''} in the \code{summary} input then the
  default from the object is ignored and the new rotation
  from the list is used instead}

  \item{SE}{logical; display the standard errors?}

  \item{printvalue}{a numeric value to be specified when
  using the \code{res='exp'} option. Only prints patterns
  that have standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{print}{logical; print output to console?}

  \item{x}{an object of class \code{polymirtClass} to be
  plotted or printed}

  \item{object}{a model estimated from \code{polymirtClass}
  of class \code{polymirt}}

  \item{object2}{a model estimated from \code{polymirt} of
  class \code{polymirtClass}}

  \item{suppress}{a numeric value indicating which
  (possibly rotated) factor loadings should be suppressed.
  Typical values are around .3 in most statistical
  software}

  \item{digits}{the number of significant digits to be
  rounded}

  \item{npts}{number of quadrature points to be used for
  plotting features. Larger values make plots look
  smoother}

  \item{rot}{allows rotation of the 3D graphics}

  \item{verbose}{logical; display iteration history during
  estimation?}

  \item{calcLL}{logical; calculate the log-likelihood?}

  \item{restype}{type of residuals to be displayed. Can be
  either \code{'LD'} for a local dependence matrix (Chen &
  Thissen, 1997) or \code{'exp'} for the expected values
  for the frequencies of every response pattern}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood}

  \item{type}{either \code{'info'} or \code{'infocontour'}
  to plot test information plots}

  \item{Target}{a dummy variable matrix indicing a target
  rotation pattern}

  \item{debug}{logical; turn on debugging features?}

  \item{technical}{list specifying subtle parameters that
  can be adjusted. Can be \describe{ \item{NCYCLES}{max
  number of MH-RM cycles; default 2000}
  \item{BURNIN}{number of burn in cycles (stage 1); default
  150} \item{SEMCYCLES}{number of SEM cycles (stage 2);
  default 50} \item{KDRAWS}{number of paralell MH sets to
  be drawn; default 1} \item{TOL}{minimum threshold
  tolerance for convergence of MH-RM, must occur on three
  consecutive occations; default .001} }}

  \item{...}{additional arguments to be passed to the
  \code{\link{confmirt}} estimation engine}
}
\description{
  \code{polymirt} fits an unconditional (exploratory)
  full-information maximum-likelihood factor analysis model
  to dichotomous and polychotomous data under the item
  response theory paradigm using Cai's (2010)
  Metropolis-Hastings Robbins-Monro algorithm. If
  requested, lower asymptote parameters are estimated with
  a beta prior included automatically.
}
\details{
  \code{polymirt} follows the item factor analysis strategy
  by a stochastic version of maximum likelihood estimation
  described by Cai (2010). The general equation used for
  multidimensional item response theory in this package is
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

  Use of \code{plot} will display the test information
  function for 1 and 2 dimensional solutions. To examine
  individuals item plots use \code{\link{itemplot}}
  (although the \code{\link[plink]{plink}} package is much
  more suitable for IRT graphics) which will also plot
  information and surface functions.

  \code{coef} displays the item parameters with their
  associated standard errors, while use of \code{summary}
  transforms the slopes into a factor loadings metric.
  Also, factor loading values below a specified constant
  can be also be suppressed in \code{summary} to allow
  better visual clarity. Models may be compared by using
  the \code{anova} function, where a Chi-squared difference
  test and AIC difference values are displayed.
}
\examples{
\dontrun{
#load LSAT section 7 data and compute 1 and 2 factor models
data(LSAT7)
fulldata <- expand.table(LSAT7)

(mod1 <- polymirt(fulldata, 1))
summary(mod1)
residuals(mod1)

(mod2 <- polymirt(fulldata, 2))
summary(mod2)
coef(mod2)
anova(mod1,mod2)

###########
#data from the 'ltm' package in numeric format
data(Science)
(mod1 <- polymirt(Science, 1))
summary(mod1)
residuals(mod1)
coef(mod1)

(mod2 <- polymirt(Science, 2, calcLL = FALSE)) #don't calculate log-likelihood
mod2 <- logLik(mod2,5000) #calc log-likelihood here with more draws
summary(mod2, 'promax', suppress = .3)
coef(mod2)
anova(mod1,mod2)


     }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Cai, L. (2010). High-Dimensional exploratory item factor
  analysis by a Metropolis-Hastings Robbins-Monro
  algorithm. \emph{Psychometrika, 75}, 33-57.

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
  \code{\link{polymirt}}, \code{\link{itemplot}}
}
\keyword{models}

