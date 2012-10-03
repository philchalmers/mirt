#' Displays item surface and information plots
#' 
#' \code{itemplot} displays various item based IRT plots.
#' 
#'
#' @aliases itemplot
#' @param object a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or 
#' \code{MultipleGroupClass}
#' @param item a single numeric value indicating which item to plot
#' @param type plot type to use, information (\code{'info'}), or item trace lines (\code{'trace'}),
#' or information contours \code{('infocontour')} (not for \code{MultipleGroupClass} objects)
#' @param degrees the degrees argument to be used if there are exactly two factors. See \code{\link{iteminfo}}
#' for more detail
#' @param ... additional arguments to be passed to lattice 
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords plot
#' @export itemplot
#' @examples
#' \dontrun{
#' 
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7)
#' mod1 <- mirt(fulldata,1)
#' 
#' itemplot(mod1, 2)
#' itemplot(mod1, 2, type = 'info')
#'     }
#' 
itemplot <- function(object, item, type = 'trace', degrees = 45, ...){
    ret <- itemplot.internal(object=object, item=item, type=type, degrees=degrees, ...)
    if(is.null(ret)) return(invisible(ret))
    else return(ret)
}


#------------------------------------------------------------------------------
setGeneric("itemplot.internal", 
           def = function(object, ...) standardGeneric("itemplot.internal")
)

setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'ExploratoryClass'),
	definition = function(object, item, type = 'trace', degrees = 45, ...)
	{  			
		x <- itemplot.main(object, item, type, degrees, ...)		        
		return(invisible(x))
	}
)

#------------------------------------------------------------------------------
setMethod(
	f = "itemplot.internal",
	signature = signature(object = 'ConfirmatoryClass'),
	definition = function(object, item, type = 'trace', degrees = 45, ...)
	{
	    x <- itemplot.main(object, item, type, degrees, ...)    	
	    return(invisible(x))
	}
)

#------------------------------------------------------------------------------
setMethod(
    f = "itemplot.internal",
    signature = signature(object = 'MultipleGroupClass'),
    definition = function(object, item, type = 'trace', degrees = 45, ...)
    {           
        Pinfo <- list()        
        gnames <- object@groupNames
        nfact <- object@nfact        
        K <- object@cmods[[1]]@pars[[item]]@ncat
        for(g in 1:length(gnames)){
            Pinfo[[g]] <- itemplot.main(object@cmods[[g]], item=item, type='RETURN', 
                                        degrees=degrees, ...)
            Pinfo[[g]]$group <- rep(gnames[g], nrow(Pinfo[[g]]))
        }
        dat <- Pinfo[[1]]        
        for(i in 2:length(gnames))
            dat <- rbind(dat, Pinfo[[g]])           
        Plist <- unclass(dat[, 1:K])
        P <- c()
        dat2 <- dat[, (K+1):ncol(dat)]
        for(i in 1:length(Plist))
            P <- c(P, Plist[[i]])
        for(i in 2:length(Plist))
            dat2 <- rbind(dat2, dat[, (K+1):ncol(dat)])
        dat2$P <- P
        dat2$cat <- rep(as.character(0:(length(Plist)-1)), each = nrow(dat))
        if(nfact == 1){
            if(type == 'info')            
                return(lattice::xyplot(info ~ Theta, dat, group=group, type = 'l', 
                                       auto.key = TRUE, main = paste('Information for item', item), 
                                       ylab = expression(I(theta)), xlab = expression(theta), ...))            
            if(type == 'trace')
                return(lattice::xyplot(P ~ Theta | cat, dat2, group=group, type = 'l', 
                                auto.key = TRUE, main = paste("Item", item, "Trace"), 
                                ylab = expression(P(theta)), xlab = expression(theta), ...))
        }
        if(nfact == 2){
            Names <- colnames(dat)
            Names[c(length(Names) - 2,length(Names) - 1)] <- c('Theta1', 'Theta2')
            Names2 <- colnames(dat2)
            Names2[2:3] <- c('Theta2', 'Theta1')
            colnames(dat) <- Names
            colnames(dat2) <- Names2            
            if(type == 'info')            
                return(lattice::wireframe(info ~ Theta1 + Theta2, data = dat, group=group, 
                                          main=paste("Item", item, "Information"), 
                                          zlab=expression(I(theta)), xlab=expression(theta[1]), 
                                          ylab=expression(theta[2]), 
                                          scales = list(arrows = FALSE), 
                                          auto.key = TRUE, ...))            
            if(type == 'trace')
                return(lattice::wireframe(P ~ Theta1 + Theta2|cat, data = dat2, group = group, 
                                          main = paste("Item", item, "Trace"), 
                                          zlab=expression(P(theta)), 
                                          xlab=expression(theta[1]), 
                                          ylab=expression(theta[2]), 
                                          scales = list(arrows = FALSE), 
                                          auto.key = TRUE, ...))           
        }
    }
)


itemplot.main <- function(x, item, type, degrees = 45, ...){        
    nfact <- ncol(x@F)
    if(nfact > 2) stop('Can not plot high dimensional models')
    if(nfact == 2 && is.null(degrees)) stop('Please specify a vector of angles that sum to 90')    
    theta <- seq(-4,4, length.out=40)
    Theta <- ThetaFull <- thetaComb(theta, nfact)   
    prodlist <- attr(x@pars, 'prodlist')
    if(length(prodlist) > 0)        
        ThetaFull <- prodterms(Theta,prodlist)
    P <- ProbTrace(x=x@pars[[item]], Theta=Theta)         
    info <- 0 
    if(nfact == 2){
        for(i in 1:length(degrees))
            info <- info + iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=c(degrees[i], 
                                                                             90 - degrees[i]))
    } else {
        info <- iteminfo(x=x@pars[[item]], Theta=ThetaFull, degrees=0)
    }
    if(type == 'RETURN') return(data.frame(P=P, info=info, Theta=Theta))
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
            return(lattice::wireframe(P ~ Theta1 + Theta2, data = plt2, group = time, main = paste("Item", item, "Trace"), 
                             zlab=expression(P(theta)), xlab=expression(theta[1]), ylab=expression(theta[2]), 
                             scales = list(arrows = FALSE), colorkey = TRUE, drape = TRUE, ...))            
        }        
    }    
}
