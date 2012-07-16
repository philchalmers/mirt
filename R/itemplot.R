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
	definition = function(object, ...)
	{  			
		x <- read.mirt(object)
		ret <- plot(x, ...)
		invisible(ret)		
	}
)

#------------------------------------------------------------------------------
# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'bfactorClass'),
	definition = function(object, ...)
	{
		x <- read.mirt(object)
		ret <- plot(x, ...)
		invisible(ret)
	}
)

#------------------------------------------------------------------------------
# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'polymirtClass'),
	definition = function(object, ...)
	{
		x <- read.mirt(object)
		ret <- plot(x, ...)
		invisible(ret)
	}
)

