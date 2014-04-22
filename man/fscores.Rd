\name{fscores}
\alias{fscores}
\title{Methods for Function fscores}
\usage{
fscores(object, rotate = "", full.scores = FALSE, method = "EAP",
  quadpts = NULL, response.pattern = NULL, returnER = FALSE,
  return.acov = FALSE, mean = NULL, cov = NULL, verbose = TRUE,
  scores.only = TRUE, full.scores.SE = FALSE, theta_lim = c(-4, 4),
  MI = 0)
}
\arguments{
  \item{object}{a computed model object of class
  \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
  \code{MultipleGroupClass}}

  \item{full.scores}{if \code{FALSE} (default) then a
  summary table with factor scores for each unique pattern
  is displayed. Otherwise the original data matrix is
  returned with the computed factor scores appended to the
  rightmost column}

  \item{rotate}{rotation declaration to be used when
  estimating the factor scores. If \code{""} then the
  \code{object@rotate} default value is used (only
  applicable to \code{ExploratoryClass} objects)}

  \item{method}{type of factor score estimation method. Can
  be expected a-posteriori (\code{"EAP"}), Bayes modal
  (\code{"MAP"}), weighted likelihood estimation
  (\code{"WLE"}), maximum likelihood (\code{"ML"}), or
  expected a-posteriori for sum scores (\code{"EAPsum"})}

  \item{quadpts}{number of quadratures to use per
  dimension. If not specified, the value will be determined
  from the input model object, otherwise a suitable one
  will be created which decreases as the number of
  dimensions increases (and therefore for estimates such as
  EAP, will be less accurate)}

  \item{theta_lim}{lower and upper range to evaluate latent
  trait integral for each dimension}

  \item{mean}{a vector for custom latent variable means. If
  NULL, the default for 'group' values from the computed
  mirt object will be used}

  \item{cov}{a custom matrix of the latent variable
  covariance matrix. If NULL, the default for 'group'
  values from the computed mirt object will be used}

  \item{MI}{a number indicating how many multiple
  imputation draws to perform. Default is 0, indicating
  that no MI draws will be performed}

  \item{response.pattern}{an optional argument used to
  calculate the factor scores and standard errors for a
  given response vector or matrix}

  \item{returnER}{logical; return empirical reliability
  (also known as marginal reliability) estimates as a
  numeric values?}

  \item{return.acov}{logical; return a list containing
  covariance matricies instead of factors scores?
  \code{impute = TRUE} not supported with this option}

  \item{full.scores.SE}{logical; when \code{full.scores ==
  TRUE}, also return the standard errors associated with
  each respondent? Default is \code{FALSE}}

  \item{verbose}{logical; print verbose output messages?}

  \item{scores.only}{logical; return only the factor scores
  (only applicable when \code{full.scores = TRUE})}
}
\description{
Computes MAP, EAP, ML (Embretson & Reise, 2000), EAP for
sum-scores (Thissen et al., 1995), or WLE (Warm, 1989)
factor scores with a multivariate normal prior distribution
using equally spaced quadrature. Additionally, multiple
imputation variants are possible for each estimator if a
parameter information matrix was computed, which are useful
if the sample size/number of items were small.
}
\details{
The function will return either a table with the computed
scores and standard errors, the original data matrix with
scores appended to the rightmost column, or the scores
only. By default the latent means and covariances are
determined from the estimated object, though these can be
overwritten. Iterative estimation methods can be estimated
in parallel to decrease estimation times if a
\code{\link{mirtCluster}} object is available.
}
\examples{
\dontrun{

mod <- mirt(Science, 1)
tabscores <- fscores(mod)
head(tabscores)
fullscores <- fscores(mod, full.scores = TRUE)
fullscores_with_SE <- fscores(mod, full.scores = TRUE, full.scores.SE=TRUE)
head(fullscores)
head(fullscores_with_SE)

#chage method argument to use MAP estimates
fullscores <- fscores(mod, full.scores = TRUE, method='MAP')
head(fullscores)

#calculate MAP for a given response vector
fscores(mod, method='MAP', response.pattern = c(1,2,3,4))
#or matrix
fscores(mod, method='MAP', response.pattern = rbind(c(1,2,3,4), c(2,2,1,3)))

#use custom latent variable properties (diffuse prior for MAP is very close to ML)
fscores(mod, method='MAP', cov = matrix(1000))
fscores(mod, method='ML')

#WLE estimation, run in parallel using available cores
mirtCluster()
fscores(mod, method='WLE')

#multiple imputation using 30 draws for EAP scores. Requires information matrix
mod <- mirt(Science, 1, SE=TRUE)
fscores(mod, MI = 30)

}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\references{
Embretson, S. E. & Reise, S. P. (2000). Item Response
Theory for Psychologists. Erlbaum.

Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S.
L. (1995). Item Response Theory for Scores on Tests
Including Polytomous Items with Ordered Responses.
\emph{Applied Psychological Measurement, 19}, 39-49

Warm, T. A. (1989). Weighted likelihood estimation of
ability in item response theory. \emph{Psychometrika, 54},
427-450.
}
\keyword{factor.scores}

