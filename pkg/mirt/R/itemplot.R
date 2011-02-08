itemplot <- function(object, ...)
{
  UseMethod("itemplot")
}

itemplot.mirt <- function(object, item, type = 'curve', npts = 30,
  rot = list(x = -70, y = 30, z = 10), ...)
{  
  if (!type %in% c('curve','info')) stop(type, " is not a valid plot type.")
  a <- as.matrix(object$pars[ ,1:(ncol(object$pars) - 1)])
  d <- object$pars[ ,ncol(object$pars)]
  g <- object$guess
  A <- as.matrix(sqrt(apply(a^2,1,sum)))
  B <- -d/A
  if(ncol(a) > 2 ) stop("Can't plot high dimentional solutions.\n")
  theta <- seq(-4,4,length.out=npts)
  Theta <- thetaComb(theta, ncol(a))
  P <- matrix(0, ncol=length(g), nrow = nrow(as.matrix(Theta)))
  for(i in 1:nrow(a)) P[ ,i] <- P.mirt(a[i, ],d[i],as.matrix(Theta),g[i])  
  P <- P[ ,item] 
  
  if(ncol(a) == 2){
    require(lattice)
    if(type == 'info'){
      I <- (P * (1 - P)) * A[item,]^2 
	  plt <- cbind(I,Theta)
	  wireframe(I ~ Theta[ ,1] + Theta[ ,2], data = plt, main = "Item Information", 
	    zlab = "I", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
		screen = rot)
    } else {  
	  plt <- cbind(P,Theta)			
	  wireframe(P ~ Theta[ ,1] + Theta[ ,2], data = plt, main = "Item Characteristic Surface", 
	    zlab = "P", xlab = "Theta 1", ylab = "Theta 2", scales = list(arrows = FALSE),
		screen = rot)
	}		
  } else {
    if(type == 'curve'){  
	  plot(Theta, P, type='l',main = 'Item Characteristic Curve', xlab = 'Theta', ylab='Probability')
	} else {
      I <- (P * (1 - P)) * a[item,]^2 
	  plot(Theta, I, type='l',main = 'Item Information', xlab = 'Theta', ylab='Information')
    } 	
  }  
}

itemplot.bfactor <- function(object, item, type = 'curve', npts = 30, 
  rot = list(x = -70, y = 30, z = 10), ...)
{
  if (!type %in% c('curve','info')) stop(type, " is not a valid plot type.")
  a <- as.matrix(object$pars[ ,1:(ncol(object$pars) - 1)])
  d <- object$pars[ ,ncol(object$pars)]
  g <- object$guess
  logicalfact <- object$logicalfact
  A <- as.matrix(sqrt(apply(a^2,1,sum)))[item]
  B <- -d/A  
  theta <- seq(-4,4,length.out=npts)
  Theta <- thetaComb(theta, 2)
  P <- matrix(0, ncol=length(g), nrow = nrow(Theta))
  for(i in 1:nrow(a)) P[ ,i] <- Pbfactor(a[i,],d[i],Theta,g[i],logicalfact[i, ])
  P <- P[ ,item]     
  require(lattice)
  
  if(type == 'info'){
    I <- (P * (1 - P)) * A^2 
	plt <- cbind(I,Theta)
	wireframe(I ~ Theta[,1] + Theta[,2], data = plt, main = "Item Information", 
	  zlab = "I", xlab = "General", ylab = "Specific", scales = list(arrows = FALSE))
  } else { 	    
	plt <- cbind(P,Theta)			
	wireframe(P ~ Theta[ ,1] + Theta[ ,2], data = plt, main = "Item Probability Surface", 
	  zlab = "P", xlab = "General", ylab = "Specific", scales = list(arrows = FALSE))
  }	  
  

}

