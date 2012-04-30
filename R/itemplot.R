#' Displays item surface and information plots
#' 
#' \code{itemplot} displays various item based IRT plots.
#' 
#'
#' @usage 
#' itemplot(object, ...)
#' 
#' \S4method{itemplot}{mirtClass}(object,
#'   items = NULL, ...)
#'
#' \S4method{itemplot}{bfactorClass}(object,
#'   items = NULL, ...) 
#'
#' \S4method{itemplot}{polymirtClass}(object,
#'   items = NULL, ...)
#'
#' @aliases itemplot-method itemplot,mirtClass-method 
#' itemplot,polymirtClass-method itemplot,bfactorClass-method
#' @param object a computed model of class \code{bfactorClass},
#' \code{mirtClass}, or \code{polymirtClass}
#' @param items an integer vector which item(s) to plot. Default is to plot all items
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

# Methods for Function itemplot
# 
# Plot individual items for fitted \code{mirt}, \code{bfactor}, or
# \code{polymirt} models.
# 
# 
# @name itemplot
# @docType methods
# @rdname itemplot-methods  
#' @export itemplot
# @keywords methods
setMethod(
	f = "itemplot",
	signature = signature(object = 'mirtClass'),
	definition = function(object, items = NULL, ...)
	{  			
		if(is.null(items)) items <- 1:length(object@K)
		x <- itemplot.main(object, items)
		ret <- plot(x, ...)
		invisible(ret)		
	}
)

# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, items = NULL, ...)
	{
		if(is.null(items)) items <- 1:length(object@K)
		x <- itemplot.main(object, items)
		ret <- plot(x, ...)
		invisible(ret)
	}
)

# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'polymirtClass'),
	definition = function(object, items = NULL, ...)
	{
		if(is.null(items)) items <- 1:length(object@K)
		x <- itemplot.main(object, items)
		ret <- plot(x, ...)
		invisible(ret)
	}
)

itemplot.main <- function(object, items)
{
	K <- object@K
	nitems <- length(items)
	guess <- object@guess
	if(class(object) == 'polymirtClass'){
		lambdas <- object@parlist$lambdas[items, , drop=FALSE]
		zetas <- object@parlist$zetas
	} else {	
		lambdas <- object@pars$lambdas[items, , drop=FALSE]
		zetas <- object@pars$zetas
	}
	if(class(object) == 'bfactorClass')
		lambdas <- cbind(lambdas[ ,1], rowSums(lambdas[,2:ncol(lambdas)])) 		
	nfact <- ncol(lambdas)
	pars <- matrix(NA, length(items), nfact + max(K))
	pars[ ,1:nfact] <- lambdas
	for(i in items){
		len <- K[i] - 1
		pars[i, (nfact+1):(nfact+len)] <- zetas[[i]]
		if(len == 1) pars[i, nfact + 2] <- guess[i]	
	}	
	if(all(is.na(pars[,ncol(pars)]))) pars <- pars[, -ncol(pars)]	
	index <- 1:nitems
	drm <- index[K == 2]
	grm <- index[K != 2]	
	if(nitems == 1){
		if(K[items] == 2) x <- drm(pars, dimensions = nfact) 
		else x <- plink::grm(pars, dimensions = nfact)
		return(x)
	}
	if(length(grm) > 0 && length(drm) > 0)
		pm <- plink::as.poly.mod(nitems, c('drm','grm'), list(drm, grm))
	else if(length(drm) > 0)
		pm <- plink::as.poly.mod(nitems, 'drm')
	else
		pm <- plink::as.poly.mod(nitems, 'grm')
	x <- plink::mixed(pars, K, pm, dimensions=nfact)
	x
}


