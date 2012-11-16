#' Methods for Function fscores
#' 
#' Computes MAP, EAP, WLE, or ML factor scores with a multivariate normal prior distribution.
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
#' (\code{"WLE"}), or maximum likelihood (\code{"ML"}) 
#' @param quadpts number of quadratures to use per dimension
#' @param degrees the degrees argument to be passed to \code{\link{iteminfo}}, only necessary for 
#' multidimensional models when \code{method = 'WLE'}
#' @param response.vector an optional argument used to calculate the factor scores and standard errors
#' for a given response vector that may or may not have been in the original dataset 
#' @param verbose logical; print verbose output messages?
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
#' 
#' #calculate MAP for a given response vector
#' fscores(mod, method='MAP', response.vector = c(1,2,3,4))
#'   }
#'
fscores <- function(object, rotate = '', full.scores = FALSE, method = "EAP", 
                    quadpts = NULL, response.vector = NULL, degrees = NULL, verbose = TRUE)
{
    ret <- fscores.internal(object=object, rotate=rotate, full.scores=full.scores, method=method, 
                            quadpts=quadpts, response.vector=response.vector, degrees=degrees,
                            verbose=verbose)
    ret    
}
