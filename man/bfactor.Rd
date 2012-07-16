\name{bfactor}
\alias{bfactor}
\alias{coef,bfactor-method}
\alias{fitted,bfactor-method}
\alias{residuals,bfactor-method}
\alias{summary,bfactor-method}
\title{Full-Information Item Bifactor Analysis}
\usage{
  bfactor(data, specific, guess = 0, upper = 1, SE = FALSE,
    prev.cor = NULL, par.prior = FALSE, startvalues = NULL,
    quadpts = 15, verbose = FALSE, debug = FALSE, technical
    = list(), ...)

  \S4method{summary}{bfactor}(object, digits = 3, ...)

  \S4method{coef}{bfactor}(object, digits = 3, ...)

  \S4method{fitted}{bfactor}(object, digits = 3, ...)

  \S4method{residuals}{bfactor}(object, restype = 'LD',
    digits = 3, printvalue = NULL, ...)
}
\arguments{
  \item{data}{a complete \code{matrix} or \code{data.frame}
  of item responses that consists of only 0, 1, and
  \code{NA} values to be factor analyzed. If scores have
  been recorded by the response pattern then they can be
  recoded to dichotomous format using the
  \code{\link{key2binary}} function.}

  \item{specific}{a numeric vector specifying which factor
  loads on which item. For example, if for a 4 item test
  with two specific factors, the first specific factor
  loads on the first two items and the second specific
  factor on the last two, then the vector is
  \code{c(1,1,2,2)}.}

  \item{printvalue}{a numeric value to be specified when
  using the \code{res='exp'} option. Only prints patterns
  that have standardized residuals greater than
  \code{abs(printvalue)}. The default (NULL) prints all
  response patterns}

  \item{guess}{fixed pseudo-guessing parameter. Can be
  entered as a single value to assign a global value or may
  be entered as a numeric vector for each item of length
  \code{ncol(data)}.}

  \item{upper}{fixed upper bound parameters for 4-PL model.
  Can be entered as a single value to assign a global
  guessing parameter or may be entered as a numeric vector
  corresponding to each item}

  \item{SE}{logical; estimate parameter standard errors?}

  \item{prev.cor}{uses a previously computed correlation
  matrix to be used to estimate starting values for the EM
  estimation}

  \item{par.prior}{a list declaring which items should have
  assumed priors distributions, and what these prior
  weights are. Elements are \code{slope} and \code{int} to
  specify the coefficients beta prior for the slopes and
  normal prior for the intercepts, and \code{slope.items}
  and \code{int.items} to specify which items to constrain.
  The value in \code{slope} is the \emph{p} meta-parameter
  for the beta distribution (where \emph{p} > 1 constrains
  the slopes), and the two values in \code{int} are the
  normal distribution intercept and variance. Larger values
  of the variance have less impact on the solution. For
  example, if items 2 and 3 were Heywood cases with no
  extreme item facilities, and item 4 had a very large item
  facility (say, greater than .95) then a possible
  constraint might be \code{par.prior = list(int = c(0,2),
  slope = 1.2, int.items = 4, slope.items = c(2,3))}}

  \item{startvalues}{user declared start values for
  parameters}

  \item{quadpts}{number of quadrature points per
  dimension.}

  \item{object}{a model estimated from \code{bfactor} of
  class \code{bfactorClass}}

  \item{restype}{type of residuals to be displayed. Can be
  either \code{'LD'} for a local dependence matrix (Chen &
  Thissen, 1997) or \code{'exp'} for the expected values
  for the frequencies of every response pattern}

  \item{digits}{number of significant digits to be rounded}

  \item{verbose}{logical; print observed log-likelihood
  value at each iteration?}

  \item{debug}{logical; turn on debugging features?}

  \item{technical}{a list containing lower level technical
  parameters for estimation \describe{
  \item{MAXQUAD}{maximum number of quadratures; default
  10000} \item{MSTEPMAXIT}{number of M-step iterations;
  default 25} \item{TOL}{EM convergence threshold; default
  .001} \item{NCYCLES}{maximum number of EM cycles; default
  300} \item{NOWARN}{a logical indicating whether dependent
  packages warnings should be printed; default \code{TRUE}}
  }}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{bfactor} fits a confirmatory maximum likelihood
  bifactor model to dichotomous and polychotomous data
  under the item response theory paradigm. Pseudo-guessing
  parameters may be included but must be declared as
  constant.
}
\details{


  \code{bfactor} follows the item factor analysis strategy
  explicated by Gibbons and Hedeker (1992) and Gibbons et
  al. (2007). Nested models may be compared via an
  approximate chi-squared difference test or by a reduction
  in AIC or BIC (accessible via \code{\link{anova}}); note
  that this only makes sense when comparing class
  \code{bfactorClass} models to class \code{mirtClass} or
  \code{polymirtClass}. The general equation used for item
  bifactor analysis in this package is in the logistic form
  with a scaling correction of 1.702. This correction is
  applied to allow comparison to mainstream programs such
  as TESTFACT 4 (2003) and POLYFACT.

  Unlike TESTFACT 4 (2003) initial start values are
  computed by using information from a quasi-tetrachoric
  correlation matrix, potentially with Carroll's (1945)
  adjustment for chance responses. To begin, a MINRES
  factor analysis with one factor is extracted, and the
  transformed loadings and intercepts (see \link{mirt} for
  more details) are used as starting values for the general
  factor loadings and item intercepts. Values for the
  specific factor loadings are taken to be half the
  magnitude of the extracted general factor loadings. Note
  that while the sign of the loading may be incorrect for
  specific factors (and possibly for some of the general
  factor loadings) the intercepts and general factor
  loadings will be relatively close to the final solution.
  These initial values should be an improvement over the
  TESTFACT initial starting values of 1.414 for all the
  general factor slopes, 1 for all the specific factor
  slopes, and 0 for all the intercepts.

  Factor scores are estimated assuming a normal prior
  distribution and can be appended to the input data matrix
  (\code{full.scores = TRUE}) or displayed in a summary
  table for all the unique response patterns. Fitted and
  residual values can be observed by using the
  \code{fitted} and \code{residuals} functions. To examine
  individuals item plots use \code{\link{itemplot}} which
  will also plot information and surface functions.
  Residuals are computed using the LD statistic (Chen &
  Thissen, 1997) in the lower diagonal of the matrix
  returned by \code{residuals}, and Cramer's V above the
  diagonal.
}
\examples{
\dontrun{

###load SAT12 and compute bifactor model with 3 specific factors
data(SAT12)
data <- key2binary(SAT12,
  key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
mod1 <- bfactor(data, specific)
coef(mod1)

###Try with guessing parameters added
guess <- rep(.1,32)
mod2 <- bfactor(data, specific, guess = guess)
coef(mod2)

#########
#simulate data
a <- matrix(c(
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,0.5,NA,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5,
1,NA,0.5),ncol=3,byrow=TRUE)

d <- matrix(c(
-1.0,NA,NA,
-1.5,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
2.5,1.0,-1,
3.0,2.0,-0.5,
3.0,2.0,-0.5,
3.0,2.0,-0.5,
2.5,1.0,-1,
2.0,0.0,NA,
-1.0,NA,NA,
-1.5,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
1.0,NA,NA),ncol=3,byrow=TRUE)

sigma <- diag(3)
dataset <- simdata(a,d,2000,sigma)

specific <- c(rep(1,7),rep(2,7))
simmod <- bfactor(dataset, specific)
coef(simmod)

    }
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
  Chalmers, R., P. (2012). mirt: A Multidimensional Item
  Response Theory Package for the R Environment.
  \emph{Journal of Statistical Software, 48}(6), 1-29.

  Gibbons, R. D., & Hedeker, D. R. (1992). Full-information
  Item Bi-Factor Analysis. \emph{Psychometrika, 57},
  423-436.

  Gibbons, R. D., Darrell, R. B., Hedeker, D., Weiss, D.
  J., Segawa, E., Bhaumik, D. K., Kupfer, D. J., Frank, E.,
  Grochocinski, V. J., & Stover, A. (2007).
  Full-Information item bifactor analysis of graded
  response data. \emph{Applied Psychological Measurement,
  31}, 4-19

  Carroll, J. B. (1945). The effect of difficulty and
  chance success on correlations between items and between
  tests. \emph{Psychometrika, 26}, 347-372.

  Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S.
  G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for
  Windows: Test Scoring, Item Statistics, and
  Full-information Item Factor Analysis [Computer
  software]. Lincolnwood, IL: Scientific Software
  International.
}
\seealso{
  \code{\link{expand.table}}, \code{\link{key2binary}},
  \code{\link{confmirt}}, \code{\link{fscores}}
}
\keyword{models}

