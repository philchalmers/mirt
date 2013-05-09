#' Methods for Function fscores
#'
#' Computes MAP, EAP, EAP for sum-scores, WLE, or ML factor scores with a multivariate normal
#' prior distribution. Will return either a table with the computed scores and standard errors, or
#' the original data matrix with scores appended to the rightmost column.
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
#' @param response.vector an optional argument used to calculate the factor scores and standard errors
#' for a given response vector that may or may not have been in the original dataset
#' @param returnER logical; return empirical reliability estimate as a numeric value?
#' @param verbose logical; print verbose output messages?
#' @param scores.only logical; return only the factor scores (only applicable when \code{full.scores = TRUE})
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords factor.scores
#' @export fscores
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
#' fscores(mod, method='MAP', response.vector = c(1,2,3,4))
#' 
#' #use custom latent variable properties (diffuse prior for MAP is very close to ML)
#' fscores(mod, method='MAP', cov = matrix(1000))
#' fscores(mod, method='ML')
#'   }
#'
fscores <- function(object, rotate = '', full.scores = FALSE, method = "EAP",
                    quadpts = NULL, response.vector = NULL, degrees = NULL,
                    returnER = FALSE, mean = NULL, cov = NULL, verbose = TRUE,
                    scores.only = FALSE)
{
    ret <- fscores.internal(object=object, rotate=rotate, full.scores=full.scores, method=method,
                            quadpts=quadpts, response.vector=response.vector, degrees=degrees,
                            verbose=verbose, returnER=returnER, gmean=mean, gcov=cov, 
                            scores.only=scores.only)
    ret
}
