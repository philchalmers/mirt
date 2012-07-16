#' Methods for Function fscores
#' 
#' Computes MAP, EAP, or ML factor scores for \code{mirt} and \code{bfactor} models,
#' or a stochastic approximation with a multivariate normal prior for \code{polymirt} and 
#' \code{confmirt}. Note that only the general factor scores are computed for bifactor 
#' models.
#'
#' 
#' @usage 
#' fscores(object, ...)
#' 
#' @aliases fscores-method fscores,bfactorClass-method
#' fscores,mirtClass-method fscores,polymirtClass-method
#' fscores,confmirtClass-method
#' @docType methods
#' @section Methods: \describe{ \item{fscores}{\code{signature(object =
#' "bfactorClass")}} \item{fscores}{\code{signature(object = "mirtClass")}}
#' \item{fscores}{\code{signature(object = "polymirtClass")}}
#' \item{fscores}{\code{signature(object = "confmirtClass")}} }
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @rdname fscores-methods   
#' @exportMethod fscores
#' @keywords methods
setGeneric("fscores", 
           def = function(object, ...) standardGeneric("fscores")
)


#------------------------------------------------------------------------------
#' Displays item surface and information plots
#' 
#' \code{itemplot} displays various item based IRT plots.
#' 
#'
#' @usage 
#' itemplot(object, ...)
#' 
#' \S4method{itemplot}{mirtClass}(object, ...)
#'
#' \S4method{itemplot}{bfactorClass}(object, ...) 
#'
#' \S4method{itemplot}{polymirtClass}(object, ...)
#'
#' @aliases itemplot-method itemplot,mirtClass-method 
#' itemplot,polymirtClass-method itemplot,bfactorClass-method
#' @param object a computed model of class \code{bfactorClass},
#' \code{mirtClass}, or \code{polymirtClass}
#' @param ... additional arguments to be passed on to \code{\link[plink]{plink}} generic
#' \code{plot()}. See the \code{\link[plink]{plink}} package for further details.  
#' @section Methods: \describe{ \item{itemplot}{\code{signature(object =
#' "bfactorClass")}} \item{itemplot}{\code{signature(object =
#' "mirtClass")}} \item{itemplot}{\code{signature(object =
#' "polymirtClass")}} }
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{plot}}, \code{\link{mirt}}, \code{\link{bfactor}},
#' \code{\link{polymirt}}
#' @rdname itemplot-methods  
#' @keywords plot
#' @docType methods
#' @exportMethod itemplot
#' @examples
#' 
#' \dontrun{
#' 
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7)
#' mod1 <- mirt(fulldata,1)
#' mod2 <- mirt(fulldata,2)
#' 
#' itemplot(mod1)
#' itemplot(mod1, combine = 5, auto.key=list(space="right"))
#'
#' itemplot(mod2, drape = TRUE)
#' itemplot(mod2, type = "vectorplot1")
#' 
#'     }
#' 
setGeneric("itemplot", 
           def = function(object, ...) standardGeneric("itemplot")
)

