#' Methods for Function itemplot
#' 
#' Plot individual items for fitted \code{mirt}, \code{bfactor}, or
#' \code{polymirt} models.
#' 
#' 
#' @name itemplot-methods
#' @aliases itemplot-methods itemplot,bfactorClass,numeric-method
#' itemplot,mirtClass,numeric-method itemplot,polymirtClass,numeric-method
#' @docType methods
#' @section Methods: \describe{ \item{itemplot}{\code{signature(object =
#' "bfactorClass", item = "numeric")}} \item{itemplot}{\code{signature(object =
#' "mirtClass", item = "numeric")}} \item{itemplot}{\code{signature(object =
#' "polymirtClass", item = "numeric")}} }
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @exportMethod itemplot
#' @export itemplot
#' @keywords methods
setMethod(
	f = "itemplot",
	signature = signature(object = 'mirtClass', item = 'numeric'),
	definition = function(object, item, type = 'info', npts = 50,
		rot = list(x = -70, y = 30, z = 10), ...)
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
			if(type == 'contour' || type == 'infocontour') cat('No \'contour\' plots for 1-dimensional models\n')
		}  
	}
)

setMethod(
	f = "itemplot",
	signature = signature(object = 'bfactorClass', item = 'numeric'),
	definition = function(object, item, type = 'info', npts = 50, 
		rot = list(x = -70, y = 30, z = 10), ...)
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
			return(wireframe(P ~ Theta1 + Theta2, data = plt, main = paste("Item",item, "Probability Surface"), 
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

setMethod(
	f = "itemplot",
	signature = signature(object = 'polymirtClass', item = 'numeric'),
	definition = function(object, item, type = 'info', npts = 50,
		rot = list(x = -70, y = 30, z = 10), ...)
	{		
		if (!type %in% c('info','infocontour')) stop(type, " is not a valid plot type.")
		if(object@K[item] > 2){
			K <- object@K		
			nfact <- ncol(object@Theta)
			a <- as.matrix(object@pars[ ,1:nfact])
			d <- as.matrix(object@pars[ ,(nfact+1):ncol(object@pars)])			
			A <- as.matrix(sqrt(apply(a^2,1,sum)))[item,]
			nzeta <- K[item] - 1
			theta <- seq(-4,4,length.out=npts)
			Theta <- thetaComb(theta, nfact)		
			P <- P.poly(a[item,], d[item,], Theta, itemexp = FALSE)
			info <- rep(0,nrow(P))
			for(i in 1:K[item]){
				w1 <- P[,i]*(1-P[,i])*A
				w2 <- P[,i+1]*(1-P[,i+1])*A
				I <- ((w1 - w2)^2) / (P[,i] - P[,i+1]) * P[,i]
				info <- info + I
			}	
			plt <- data.frame(cbind(info,Theta))		
			if(nfact == 1)	
				plot(Theta, info, type='l',main = paste('Item', item,'Information'), 
					xlab = 'Theta', ylab='Information')
			else {					
				colnames(plt) <- c('info','Theta1','Theta2')
				if(type == 'info')
					return(wireframe(info ~ Theta1 + Theta2, data = plt, main = paste("Item",item,"Information"), 
						zlab = "I", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
						screen = rot))				
				if(type == 'infocontour'){										
					contour(theta, theta, matrix(info,length(theta),length(theta)), 
						main = paste("Item", item,"Information Contour"), xlab = "Theta 1", ylab = "Theta 2")					
				}
			}	
		} else {
			class(object) <- 'mirtClass'
			itemplot(object,item,type,npts,rot)		 
		}	
	}
)

