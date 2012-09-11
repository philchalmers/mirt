#' Displays item surface and information plots
#' 
#' \code{itemplot} displays various item based IRT plots.
#' 
#'
#' @usage 
#' itemplot(object, ...)
#' 
#' \S4method{itemplot}{ExploratoryClass}(object, item, type = 'trace', degrees = NULL, ...)
#'
#' \S4method{itemplot}{ConfirmatoryClass}(object, item, type = 'trace' , degrees = NULL, ...) 
#'
#' @aliases itemplot-method itemplot,ExploratoryClass-method 
#' itemplot,ConfirmatoryClass-method
#' @param object a computed model object
#' @param item a single numeric value indicating which item to plot
#' @param type plot type to use, information (\code{'info'}), information contours \code{('infocontour')},
#'  or item trace lines (\code{'trace'})
#' @param degrees the degrees argument to be used if there are exactly two factors. See \code{\link{iteminfo}}
#' for more detail
#' @param ... additional arguments to be passed to lattice 
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
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
#' 
#' itemplot(mod1, 2)
#'     }
#' 
setGeneric("itemplot", 
           def = function(object, ...) standardGeneric("itemplot")
)

#------------------------------------------------------------------------------
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
	signature = signature(object = 'ExploratoryClass'),
	definition = function(object, item, type = 'trace', degrees = 45, ...)
	{  			
		x <- itemplot.main(object, item, type, degrees, ...)		        
		return(invisible(x))
	}
)

#------------------------------------------------------------------------------
# @rdname itemplot-methods  
setMethod(
	f = "itemplot",
	signature = signature(object = 'ConfirmatoryClass'),
	definition = function(object, item, type = 'trace', degrees = 45, ...)
	{
	    x <- itemplot.main(object, item, type, degrees, ...)    	
	    return(invisible(x))
	}
)

itemplot.main <- function(x, item, type, degrees = 45, ...){        
    nfact <- ncol(x@F)
    if(nfact > 2 && !x[[1]]@bfactor) stop('Can not plot high dimensional models')
    if(nfact == 2 && is.null(degrees)) stop('Please specify a vector of angles that sum to 90')    
    theta <- seq(-4,4, length.out=40)
    Theta <- thetaComb(theta, nfact)   
    P <- ProbTrace(x=x@pars[[item]], Theta=Theta)         
    info <- 0 
    if(nfact == 2){
        for(i in 1:length(degrees))
            info <- info + iteminfo(x=x@pars[[item]], Theta=Theta, degrees=c(degrees[i], 
                                                                             90 - degrees[i]))
    }
    if(nfact == 1){
        if(type == 'trace'){            
            plot(Theta, P[,1], col = 1, type='l', main = paste('Item', item), 
                 ylab = expression(P(theta)), xlab = expression(theta), ylim = c(0,1), las = 1, 
                 ...)
            for(i in 2:ncol(P))
                lines(Theta, P[,i], col = i)                 
        }
        if(type == 'info'){            
            plot(Theta, info, col = 1, type='l', main = paste('Information for item', item), 
                 ylab = expression(I(theta)), xlab = expression(theta), las = 1)
        }
        if(type == 'infocontour') stop('Cannot draw contours for 1 factor models')        
    } else {        
        plt <- data.frame(info = info, Theta1 = Theta[,1], Theta2 = Theta[,2])
        plt2 <- data.frame(P = P, Theta1 = Theta[,1], Theta2 = Theta[,2])
        colnames(plt2) <- c(paste("P", 1:ncol(P), sep=''), "Theta1", "Theta2")
        plt2 <- reshape(plt2, direction='long', varying = paste("P", 1:ncol(P), sep=''), v.names = 'P', 
                times = paste("P", 1:ncol(P), sep=''))
        colnames(plt) <- c("info", "Theta1", "Theta2")            
        if(type == 'infocontour')												
            return(contourplot(info ~ Theta1 * Theta2, data = plt, 
                               main = paste("Item", item, "Information Contour"), xlab = expression(theta[1]), 
                               ylab = expression(theta[2]), ...))
        if(type == 'info')
            return(lattice::wireframe(info ~ Theta1 + Theta2, data = plt, main = paste("Item", item, "Information"), 
                             zlab=expression(I(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                             scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...))
        if(type == 'trace'){
            return(lattice::wireframe(P ~ Theta1 + Theta2|time, data = plt2, main = paste("Item", item, "Trace"), 
                             zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                             scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...))            
        }        
    }    
}
