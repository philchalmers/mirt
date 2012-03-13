#' Displays item surface and information plots
#' 
#' \code{itemplot} displays 3D surface plots if the number of factors is two, or
#' standard 2D plots if the number of factors is equal to one.
#' 
#'
#' @usage 
#' itemplot(object, item, ...)
#' 
#' \S4method{itemplot}{mirtClass,numeric}(object,
#'   item, type = "info", npts = 50, rot = list(), ...)
#'
#' \S4method{itemplot}{bfactorClass,numeric}(object,
#'   item, type = "info", npts = 50, rot = list(), ...)
#'
#' \S4method{itemplot}{polymirtClass,numeric}(object,
#'   item, type = "info", npts = 50, rot = list(), ...)
#' @aliases itemplot-method itemplot,mirtClass,numeric-method 
#' itemplot,polymirtClass,numeric-method itemplot,bfactorClass,numeric-method
#' @param object a computed model of class \code{bfactorClass},
#' \code{mirtClass}, or \code{polymirtClass}
#' @param item a single numerical value indicating which item to plot
#' @param type type of graphic to plot, may be either \code{'curve'} for the
#' response surface curve, \code{'info'} for the item information,
#' \code{'contour'} for item probability contours, or \code{'infocontour'} for
#' information contour plot. Note that \code{polymirtClass} objects can only
#' display information plots
#' @param npts number of points to use per dimension. A higher value will make
#' the graphs look smoother. Default is 50
#' @param rot allows rotation of the 3D graphics. Default is \code{list(x=-70,y=30,z=10)}
#' @param ... additional arguments to be passed
#' @section Methods: \describe{ \item{itemplot}{\code{signature(object =
#' "bfactorClass", item = "numeric")}} \item{itemplot}{\code{signature(object =
#' "mirtClass", item = "numeric")}} \item{itemplot}{\code{signature(object =
#' "polymirtClass", item = "numeric")}} }
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
	signature = signature(object = 'mirtClass', item = 'numeric'),
	definition = function(object, item, type = 'info', npts = 50,
		rot = list(), ...)
	{  
		if (!type %in% c('curve','info','contour','infocontour')) stop(type, " is not a valid plot type.")
		nfact <- ncol(object@Theta)
		a <- as.matrix(object@pars[ ,1:nfact])
		d <- object@pars[ ,ncol(object@pars)]
		g <- object@guess
		g[is.na(g)] <- 0
		A <- as.matrix(sqrt(apply(a^2,1,sum)))
		B <- -d/A
		if(nfact > 2 ) stop("Can't plot high dimentional solutions.\n")
		theta <- seq(-4,4,length.out=npts)
		Theta <- thetaComb(theta, nfact) 
		P <- P.mirt(a[item, ],d[item],as.matrix(Theta),g[item])  
		Pstar <- P.mirt(a[item, ],d[item],as.matrix(Theta),0)  		  
		if(nfact == 2){			
			I <- (P * (1 - P)) * Pstar/P * A[item,]^2
			if(type == 'info'){				 
				plt <- data.frame(cbind(I,Theta))
				colnames(plt) <- c('I','Theta1','Theta2')
				return(wireframe(I ~ Theta1 + Theta2, data = plt, main = paste("Item", item,"Information"),
					zlab = "I", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
					screen = rot))
			} 
			if(type == 'curve'){	
				plt <- data.frame(cbind(P,Theta))
				colnames(plt) <- c('P','Theta1','Theta2')
				return(wireframe(P ~ Theta1 + Theta2, data = plt, main = paste("Item", item,
					"Characteristic Surface"), zlab = "P", xlab = "Theta 1", ylab = "Theta 2", 
					scales = list(arrows = FALSE), screen = rot))
			}		
			if(type == 'contour'){
				plt <- data.frame(cbind(P,Theta))
				colnames(plt) <- c('P','Theta1','Theta2')
				contour(theta, theta, matrix(P,length(theta),length(theta)), 
					main = paste("Item", item,"Probability Contour"), xlab = "Theta 1", ylab = "Theta 2")
			}
			if(type == 'infocontour'){
				plt <- data.frame(cbind(I,Theta))
				colnames(plt) <- c('I','Theta1','Theta2')
				contour(theta, theta, matrix(I,length(theta),length(theta)), 
					main = paste("Item", item,"Information Contour"), xlab = "Theta 1", ylab = "Theta 2")
			}
		} else {
			if(type == 'curve')  
				plot(Theta, P, type='l',main = paste("Item", item, "Characteristic Curve"), 
					xlab = 'Theta', ylab='Probability')
			if(type == 'info'){ 
				I <- (P * (1 - P)) * Pstar/P * a[item,]^2 
				plot(Theta, I, type='l',main = paste('Item', item, 'Information'), xlab = 'Theta', 
					ylab='Information')
			}
			if(type=='contour' || type=='infocontour') cat('No \'contour\' plots for 1-dimensional models\n')
		}  
	}
)

# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'bfactorClass', item = 'numeric'),
	definition = function(object, item, type = 'info', npts = 50, 
		rot = list(), ...)
	{
		if (!type %in% c('curve','info','contour','infocontour')) stop(type, " is not a valid plot type.")
		a <- as.matrix(object@pars[ ,1:(ncol(object@pars) - 1)])
		d <- object@pars[ ,ncol(object@pars)]
		g <- object@guess
		logicalfact <- object@logicalfact
		A <- as.matrix(sqrt(apply(a^2,1,sum)))[item]
		B <- -d/A  
		theta <- seq(-4,4,length.out=npts)
		Theta <- thetaComb(theta, 2)		
		P <- P.bfactor(a[item,],d[item],Theta,g[item],logicalfact[item, ])
		Pstar <- P.bfactor(a[item,],d[item],Theta,0,logicalfact[item, ])		    
		I <- (P * (1 - P)) * Pstar/P * A^2		
		if(type == 'info'){			 
			plt <- data.frame(cbind(I,Theta))
			colnames(plt) <- c('I','Theta1','Theta2')		
			return(wireframe(I ~ Theta1 + Theta2, data = plt, main = paste("Item",item,"Information"), 
				zlab = "I", xlab = "General", ylab = "Specific", scales = list(arrows = FALSE)))
		}		
		if(type == 'curve'){
			plt <- data.frame(cbind(P,Theta))	
			colnames(plt) <- c('P','Theta1','Theta2')		
			return(wireframe(P ~ Theta1 + Theta2, data=plt, main = paste("Item",item, "Probability Surface"), 
				zlab = "P", xlab = "General", ylab = "Specific", scales = list(arrows = FALSE)))
		}	  
		if(type == 'contour'){
			plt <- data.frame(cbind(P,Theta))
			colnames(plt) <- c('P','Theta1','Theta2')
			contour(theta, theta, matrix(P,length(theta),length(theta)), 
				main = paste("Item", item,"Probability Contour"), xlab = "General", ylab = "Specific")	
		}
		if(type == 'infocontour'){
			plt <- data.frame(cbind(I,Theta))
			colnames(plt) <- c('I','Theta1','Theta2')
			contour(theta, theta, matrix(I,length(theta),length(theta)), 
				main = paste("Item", item,"Information Contour"), xlab = "General", ylab = "Specific")
		}
	}
)





