#' Convert an esimated mirt model to special data.frame
#' 
#' Given an estimated model from any of mirt's model fitting functions this function will convert 
#' the model parameters into the design data frame of starting values and other parameter characteristics 
#' (similar to using the \code{pars = 'values'} for obtaining starting values).
#' 
#'
#' @aliases mod2values
#' @param x an estimated model x from the mirt package
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords convert model
#' @export mod2values
#' @examples 
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#' values <- mod2values(mod)
#' 
#' #use the converted values as starting values in a new model
#' mod2 <- mirt(dat, 1, pars = values)
#' 
#' }
mod2values <- function(x){   
    if(is(x, 'MultipleGroupClass')){
        PrepList <- x@cmods
        names(PrepList) <- x@groupNames        
        MG <- TRUE
    } else {
        PrepList <- list(pars=x@pars)
        names(PrepList) <- 'full'
        MG <- FALSE        
    }
    itemnames <- colnames(x@data)
    parnum <- par <- est <- item <- parname <- gnames <- itemtype <- 
        lbound <- ubound <- c()                                        
    for(g in 1:length(PrepList)){
        if(MG) tmpgroup <- PrepList[[g]]@pars                                
        else tmpgroup <- PrepList[[g]]
        for(i in 1:length(tmpgroup)){
            if(i <= length(itemnames))
                item <- c(item, rep(itemnames[i], length(tmpgroup[[i]]@parnum)))
            parname <- c(parname, names(tmpgroup[[i]]@parnum))
            parnum <- c(parnum, tmpgroup[[i]]@parnum) 
            par <- c(par, tmpgroup[[i]]@par)
            est <- c(est, tmpgroup[[i]]@est)
            lbound <- c(lbound, tmpgroup[[i]]@lbound)
            ubound <- c(ubound, tmpgroup[[i]]@ubound)
        }
        item <- c(item, rep('GROUP', length(tmpgroup[[i]]@parnum)))                                
    }
    gnames <- rep(names(PrepList), each = length(est)/length(PrepList))
    ret <- data.frame(group=gnames, item = item, name=parname, parnum=parnum, value=par, 
                      lbound=lbound, ubound=ubound, est=est)
    ret
}