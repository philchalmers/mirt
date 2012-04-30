#' Translate mirt parameters for plink package
#' 
#' A plotting function for displaying the individuals trajectories and their 
#' modelled functional form. Useful for detecting aberrant individual trajectories.
#' 
#' 
#' @aliases read.mirt
#' @param x an object returned from \code{mirt, bfactor, polymirt}, or \code{confmirt}
#' @param loc.out if \code{TRUE}, the step/threshold parameters will be reformated to be deviations 
#' from a location parameter
#' @param as.irt.pars if \code{TRUE}, the parameters will be output as an \code{irt.pars} object 
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords plink
#' @export read.mirt
#' @examples 
#' 
#' \dontrun{
#' data(LSAT7)
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, 1))
#' plinkpars <- read.mirt(mod1)
#' 
#' }
read.mirt <- function (x, loc.out = FALSE, as.irt.pars = TRUE) 
{
    cls <- class(x)    
	if(any(cls == c('mirtClass', 'bfactorClass', 'polymirtClass'))){
		cat <- x@K
		nitems <- length(cat)
		guess <- x@guess
		lambdas <- x@pars$lambdas
		zetas <- x@pars$zetas
		dimensions <- ncol(lambdas)
		pars <- matrix(NA, nitems, dimensions + max(cat-1,2))
		pars[ ,1:dimensions] <- lambdas
		for(i in 1:nitems){
			len <- cat[i] - 1
			pars[i, (dimensions+1):(dimensions+len)] <- zetas[[i]]
			if(len == 1) pars[i, dimensions + 2] <- guess[i]	
		}		
		if (all(cat == 2)) class(cls) <- 'drm'
		else if (all(cat > 2)) class(cls) <- 'grm'
		else class(cls) <- c('drm','grm')				
	} else if (cls == 'confmirtClass'){
		stop('confmirtClass objects not yet supported')
	}	
	if(setequal(attr(cls, "class"), c('drm','grm'))){
		index <- 1:nitems
		drm <- index[cat == 2]
		grm <- index[cat != 2]
		pm <- as.poly.mod(nrow(pars), c('drm','grm'), list(drm, grm))		
	} else pm <- as.poly.mod(nrow(pars), attr(cls, "class"))
    pars <- sep.pars(pars, cat = cat, poly.mod = pm, dimensions = dimensions, 
        loc.out = loc.out)
    pars <- as.irt.pars(pars)
    if (as.irt.pars == FALSE) {
        pars <- pars@pars
    }
    return(pars)
}
