#' Translate mirt parameters for plink package
#' 
#' A plotting function for displaying the individuals trajectories and their 
#' modelled functional form. Useful for detecting aberrant individual trajectories.
#' 
#' 
#' @aliases read.mirt
#' @param x an object returned from \code{mirt, bfactor, confmirt}, or \code{multipleGroup}
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
read.mirt <- function (x, as.irt.pars = TRUE) 
{    
    if(!require(plink)) stop('You must load the plink package.')
    cls <- class(x) 
    if(class(x) == 'MultipleGroupClass'){ 
        pars <- vector(length(x@cmods))
        for(i in 1:length(pars))
            pars[[i]] <- read.mirt(x@cmods, as.irt.pars=as.irt.pars)
        return(pars)
    }
    if(length(x@prodlist)) 
        stop('Polynomial factor models not supported in plink')
    listpars <- coef(x)    
    nitems <- length(listpars) - 1
    mirt.items <- as.character(lapply(x@pars, class))    
    mirt.items <- plink.items <- mirt.items[-(nitems+1)]
    cat <- numeric(nitems)
    D <- x@pars[[1]]@D               
    nfact <- x@pars[[1]]@nfact
    pars <- matrix(NA, nitems, 40)
    for(i in 1:nitems){
        
        if(mirt.items[i] == 'dich'){            
            plink.items[i] <- 'drm' 
            cat[i] <- 2
            abc <- listpars[[i]][1:(nfact+2)]
            if(nfact == 1)
                abc[2] <- -abc[2]/(abc[1] * D)
            pars[i, 1:length(abc)] <- abc
        }
        
        if(mirt.items[i] == 'graded'){
            plink.items[i] <- 'grm'
            ab <- listpars[[i]]
            cat[i] <- x@pars[[i]]@ncat 
            if(nfact == 1)
                for(j in 1:(cat[i] - 1))
                    ab[j + 1] <- -ab[j+1]/(ab[1] * D)
            pars[i, 1:length(ab)] <- ab            
        }
        
        if(mirt.items[i] == 'rating'){
            plink.items[i] <- 'grm'            
            ab <- listpars[[i]]
            adj <- ab[length(ab)]
            ab <- ab[-length(ab)]
            cat[i] <- x@pars[[i]]@ncat 
            for(j in 1:(cat[i] - 1))
                ab[j + nfact] <- ab[j+1] + adj
            if(nfact == 1)
                for(j in 1:(cat[i] - 1))
                    ab[j + 1] <- -ab[j+1]/(ab[1] * D)
            pars[i, 1:length(ab)] <- ab 
        }
        
        if(mirt.items[i] == 'gpcm'){
            ab <- listpars[[i]]
            a <- ab[1:nfact]
            ab <- ab[-(1:nfact)]
            cat[i] <- x@pars[[i]]@ncat 
            newab <- c()
            for(j in 1:nfact)              
                newab <- c(newab, a[j] * ab[1:cat[i]])
            newab <- c(newab, ab[(length(ab)-(cat[i]-1)):length(ab)])
            pars[i, 1:length(newab)] <- newab 
        }
        
        if(mirt.items[i] == 'nominal'){            
            plink.items[i] <- 'nrm'
            ab <- listpars[[i]]
            a <- ab[1:nfact]
            ab <- ab[-(1:nfact)]
            cat[i] <- x@pars[[i]]@ncat             
            newab <- c()
            for(j in 1:nfact)              
                newab <- c(newab, a[j] * ab[1:cat[i]])
            newab <- c(newab, ab[(length(ab)-(cat[i]-1)):length(ab)])
            pars[i, 1:length(newab)] <- newab
        }
        
        if(mirt.items[i] == 'mcm'){                                     
            acd <- listpars[[i]]
            a <- acd[1:nfact]
            acd <- acd[-(1:nfact)]
            cat[i] <- x@pars[[i]]@ncat 
            newacd <- c()
            for(j in 1:nfact)              
                newacd <- c(newacd, a[j] * acd[1:(cat[i]+1)])
            newacd <- c(newacd, acd[(length(acd) - 2*cat[i]):length(acd)])                                              
            pars[i, 1:length(newacd)] <- newacd
        }
        
        if(mirt.items[i] == 'partcomp'){
            stop('Partially compensatory models not supported in plink')
        }
        
    }    
    model <- unique(plink.items)
    items <- vector('list', length(model))
    index <- 1:nitems
    for(i in 1:length(model))
        items[[i]] <- index[model[i] == plink.items]
    pm <- as.poly.mod(nitems, model=model, items=items)
    pars <- pars[ , colSums(is.na(pars)) != nitems]
    if(as.irt.pars) 
        pars <- as.irt.pars(pars, cat=cat, poly.mod=pm, dimensions=nfact)
    return(pars)
}

