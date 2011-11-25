#Turned off until roxygen2 has support for @slots and methods

# Class "mirtClass"
# 
# Defines the object returned from \code{\link{mirt}}.
# 
# 
# @name mirtClass-class
# @aliases mirtClass-class anova,mirtClass-method coef,mirtClass-method
# fitted,mirtClass-method plot,mirtClass,missing-method print,mirtClass-method
# residuals,mirtClass-method show,mirtClass-method summary,mirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("mirtClass", ...).}.
# @method Emiter number of EM iterations
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass mirtClass
# @keywords classes
setClass(
	Class = 'mirtClass',
	representation = representation(EMiter = 'numeric', pars = 'matrix', guess = 'numeric', 
		X2 = 'numeric', df = 'numeric', p = 'numeric', AIC = 'numeric', log.lik = 'numeric',
		F = 'matrix', h2 = 'numeric', tabdata = 'matrix', Theta = 'matrix', Pl = 'numeric',
		fulldata = 'matrix', cormat = 'matrix', facility = 'numeric', converge = 'numeric', 
		quadpts = 'numeric', BIC = 'numeric', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

# Class "bfactorClass"
# 
# Defines the object returned from \code{\link{bfactor}}.
# 
# 
# @name bfactorClass-class
# @aliases bfactorClass-class coef,bfactorClass-method
# fitted,bfactorClass-method print,bfactorClass-method
# residuals,bfactorClass-method show,bfactorClass-method
# summary,bfactorClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("bfactorClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass bfactorClass
# @keywords classes
setClass(
	Class = 'bfactorClass',
	representation = representation(EMiter = 'numeric', pars = 'matrix', guess = 'numeric', 
		AIC = 'numeric', X2 = 'numeric', df = 'numeric', log.lik = 'numeric', p = 'numeric', 
		F = 'matrix', h2 = 'numeric', itemnames = 'character', tabdata = 'matrix', 
		N = 'numeric', Pl = 'numeric', Theta = 'matrix', fulldata = 'matrix', 
		logicalfact = 'matrix', facility = 'numeric', specific = 'numeric', BIC = 'numeric',
		cormat = 'matrix', converge = 'numeric', par.prior = 'matrix', quadpts = 'numeric', 
		Call = 'call'),	
	validity = function(object) return(TRUE)
)	

# Class "polymirtClass"
# 
# Defines the object returned from \code{\link{polymirt}}.
# 
# 
# @name polymirtClass-class
# @aliases polymirtClass-class coef,polymirtClass-method
# plot,polymirtClass,missing-method print,polymirtClass-method
# residuals,polymirtClass-method show,polymirtClass-method
# summary,polymirtClass-method anova,polymirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("polymirtClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass polymirtClass
# @keywords classes
setClass(
	Class = 'polymirtClass',
	representation = representation(pars = 'matrix', guess = 'numeric', SEpars = 'matrix', 
		cycles = 'numeric', Theta = 'matrix', fulldata = 'matrix', data = 'matrix', 
		K = 'numeric', F = 'matrix', h2 = 'numeric', itemloc = 'numeric', AIC = 'numeric',
		converge = 'numeric', logLik = 'numeric', SElogLik = 'numeric', df = 'integer', 
		G2 = 'numeric', p = 'numeric', tabdata = 'matrix', BIC = 'numeric', estGuess = 'logical', 
		Call = 'call'),	
	validity = function(object) return(TRUE)
)	

# Class "confmirtClass"
# 
# Defines the object returned from \code{\link{confmirt}}.
# 
# 
# @name confmirtClass-class
# @aliases confmirtClass-class coef,confmirtClass-method
# print,confmirtClass-method residuals,confmirtClass-method
# show,confmirtClass-method summary,confmirtClass-method
# anova,confmirtClass-method
# @docType class
# @section Objects from the Class: Objects can be created by calls of the form
# \code{new("confmirtClass", ...)}.
# @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportClass confmirtClass
# @keywords classes
setClass(
	Class = 'confmirtClass',
	representation = representation(pars = 'matrix', guess = 'numeric', SEpars = 'matrix', 
		SEg = 'numeric', gpars = 'list', SEgpars = 'list', estpars = 'list',cycles = 'numeric', 
		Theta = 'matrix', fulldata = 'matrix', data = 'matrix', K = 'numeric', itemloc = 'numeric',
		h2 = 'numeric',F = 'matrix', converge = 'numeric', logLik = 'numeric',SElogLik = 'numeric',
		df = 'integer', AIC = 'numeric', nconstvalues = 'integer', G2 = 'numeric', p = 'numeric',
		tabdata = 'matrix', BIC = 'numeric', estComp = 'logical', Call = 'call'),	
	validity = function(object) return(TRUE)
)	

#' Methods for Function fscores
#' 
#' Save tabulated or full data factor scores for \code{mirt} or \code{bfactor}
#' using EAP or MAP scoring.
#' 
#' 
#' @name fscores-methods
#' @aliases fscores-methods fscores,bfactorClass-method
#' fscores,mirtClass-method fscores,polymirtClass-method
#' fscores,confmirtClass-method
#' @docType methods
#' @section Methods: \describe{ \item{fscores}{\code{signature(object =
#' "bfactorClass")}} \item{fscores}{\code{signature(object = "mirtClass")}}
#' \item{fscores}{\code{signature(object = "polymirtClass")}}
#' \item{fscores}{\code{signature(object = "confmirtClass")}} }
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportMethod fscores
#' @keywords methods
setGeneric("fscores", 
	def = function(object, ...) standardGeneric("fscores")
)

#' Displays item surface and information plots
#' 
#' \code{itemplot} displays 3D surface plots if the number of factors is 2, or
#' standard 2D plots if the number of factors is equal to one.
#' 
#' 
#' @aliases itemplot itemplot,mirt-method itemplot,bfactor-method
#' itemplot,polymirt-method
#' @param object a computed model of class \code{bfactorClass},
#' \code{mirtClass}, or \code{polymirtClass}
#' @param item a single numerical value indicating which item to plot
#' @param type type of graphic to plot, may be either \code{'curve'} for the
#' response surface curve, \code{'info'} for the item information,
#' \code{'contour'} for item probability contours, or \code{'infocontour'} for
#' information contour plot. Note that \code{polymirtClass} objects can only
#' display information plots
#' @param npts number of points to use per dimension. A higher value will make
#' the graphs look smoother
#' @param rot allows rotation of the 3D graphics
#' @param ... additional arguments to be passed
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{plot}}
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
#' #first item unidimensional
#' itemplot(mod1, 1)
#' itemplot(mod1, 1, type = 'info')
#' 
#' #first item multidimensional
#' itemplot(mod2, 1)
#' itemplot(mod2, 1, type = 'info')
#' #turn it a little
#' itemplot(mod2, 1, type = 'info', rot = list(x = -70, y = -15, z = -10))
#' 
#'     }
#' 
setGeneric("itemplot", 
	def = function(object, item, ...) standardGeneric("itemplot")
)

