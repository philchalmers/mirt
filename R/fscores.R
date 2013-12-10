#' Methods for Function fscores
#'
#' Computes MAP, EAP, ML (Embretson & Reise, 2000), EAP for sum-scores (Thissen et al., 1995),
#' or WLE (Warm, 1989) factor scores with a multivariate normal
#' prior distribution. Will return either a table with the computed scores and standard errors,
#' the original data matrix with scores appended to the rightmost column, or the scores only. By
#' default the latent means are set to be 0 for each factor, and the covariance matrix is set to the
#' identity matrix, though these can be overwritten. Iterative estimation methods can be estimated
#' in parallel to decrease estimation times if a \code{\link{mirtCluster}} object is available.
#'
#'
#' @aliases fscores
#' @param object a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
#' \code{MultipleGroupClass}
#' @param full.scores if \code{FALSE} (default) then a summary table with
#' factor scores for each unique pattern is displayed. Otherwise the original
#' data matrix is returned with the computed factor scores appended to the
#' rightmost column
#' @param rotate rotation declaration to be used when estimating the factor scores. If \code{""} then the
#' \code{object@@rotate} default value is used (only applicable to \code{ExploratoryClass} objects)
#' @param method type of factor score estimation method. Can be expected
#' a-posteriori (\code{"EAP"}), Bayes modal (\code{"MAP"}), weighted likelihood estimation
#' (\code{"WLE"}), maximum likelihood (\code{"ML"}), or expected a-posteriori for sum scores (\code{"EAPsum"})
#' @param quadpts number of quadratures to use per dimension
#' @param mean a vector for custom latent variable means. If NULL, the default for 'group' values from the
#' computed mirt object will be used
#' @param cov a custom matrix of the latent variable covariance matrix. If NULL, the default for 'group' values
#' from the computed mirt object will be used
#' @param degrees the degrees argument to be passed to \code{\link{iteminfo}}, only necessary for
#' multidimensional models when \code{method = 'WLE'}
#' @param response.pattern an optional argument used to calculate the factor scores and standard errors
#' for a given response vector or matrix
#' @param returnER logical; return empirical reliability (also known as marginal reliability)
#' estimates as a numeric values?
#' @param verbose logical; print verbose output messages?
#' @param scores.only logical; return only the factor scores (only applicable when \code{full.scores = TRUE})
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords factor.scores
#' @export fscores
#' @references
#'
#' Embretson, S. E. & Reise, S. P. (2000). Item Response Theory for Psychologists. Erlbaum.
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995). Item Response Theory for Scores
#' on Tests Including Polytomous Items with Ordered Responses. \emph{Applied Psychological
#' Measurement, 19}, 39-49
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
#' \emph{Psychometrika, 54}, 427-450.
#'
#' @examples
#'
#' \dontrun{
#'
#' mod <- mirt(Science, 1)
#' tabscores <- fscores(mod)
#' fullscores <- fscores(mod, full.scores = TRUE)
#' fullscores <- fscores(mod, full.scores = TRUE, method='MAP')
#' scores <- fscores(mod, full.scores = TRUE, scores.only = TRUE)
#' head(scores)
#'
#' #calculate MAP for a given response vector
#' fscores(mod, method='MAP', response.pattern = c(1,2,3,4))
#' #or matrix
#' fscores(mod, method='MAP', response.pattern = rbind(c(1,2,3,4), c(2,2,1,3)))
#'
#' #use custom latent variable properties (diffuse prior for MAP is very close to ML)
#' fscores(mod, method='MAP', cov = matrix(1000))
#' fscores(mod, method='ML')
#'
#' #WLE estimation, run in parallel using available cores
#' mirtCluster()
#' fscores(mod, method='WLE')
#'
#'   }
fscores <- function(object, rotate = '', full.scores = FALSE, method = "EAP",
                    quadpts = NULL, response.pattern = NULL, degrees = NULL,
                    returnER = FALSE, mean = NULL, cov = NULL, verbose = TRUE,
                    scores.only = FALSE)
{
    ret <- fscores.internal(object=object, rotate=rotate, full.scores=full.scores, method=method,
                            quadpts=quadpts, response.pattern=response.pattern, degrees=degrees,
                            verbose=verbose, returnER=returnER, gmean=mean, gcov=cov,
                            scores.only=scores.only)
    ret
}
