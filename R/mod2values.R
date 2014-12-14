#' Convert an estimated mirt model to special data.frame
#'
#' Given an estimated model from any of mirt's model fitting functions this function will convert
#' the model parameters into the design data frame of starting values and other parameter
#' characteristics (similar to using the \code{pars = 'values'} for obtaining starting values).
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
#' values
#'
#' #use the converted values as starting values in a new model, and reduce TOL
#' mod2 <- mirt(dat, 1, pars = values, TOL=1e-5)
#'
#' }
mod2values <- function(x){
    if(is(x, 'MultipleGroupClass') || is(x, 'DiscreteClass')){
        PrepList <- x@pars
        names(PrepList) <- x@Data$groupNames
        MG <- TRUE
    } else {
        PrepList <- list(pars=x@pars)
        names(PrepList) <- 'all'
        MG <- FALSE
    }
    itemnames <- colnames(x@Data$data)
    parnum <- par <- est <- item <- parname <- gnames <- class <-
        lbound <- ubound <- prior.type <- prior_1 <- prior_2 <- c()
    for(g in 1L:length(PrepList)){
        if(MG) tmpgroup <- PrepList[[g]]@pars
        else tmpgroup <- PrepList[[g]]
        for(i in 1L:length(tmpgroup)){
            if(i <= length(itemnames))
                item <- c(item, rep(itemnames[i], length(tmpgroup[[i]]@parnum)))
            class <- c(class, rep(class(tmpgroup[[i]]), length(tmpgroup[[i]]@parnum)))
            parname <- c(parname, names(tmpgroup[[i]]@est))
            parnum <- c(parnum, tmpgroup[[i]]@parnum)
            par <- c(par, tmpgroup[[i]]@par)
            est <- c(est, tmpgroup[[i]]@est)
            lbound <- c(lbound, tmpgroup[[i]]@lbound)
            ubound <- c(ubound, tmpgroup[[i]]@ubound)
            prior.type <- c(prior.type, tmpgroup[[i]]@prior.type)
            prior_1 <- c(prior_1, tmpgroup[[i]]@prior_1)
            prior_2 <- c(prior_2, tmpgroup[[i]]@prior_2)
        }
        item <- c(item, rep('GROUP', length(tmpgroup[[i]]@parnum)))
    }
    if(is(x, 'MixedClass')){
        tmpgroup <- x@random
        if(length(tmpgroup)){
            for(i in 1L:length(tmpgroup)){
                parname <- c(parname, names(tmpgroup[[i]]@est))
                parnum <- c(parnum, tmpgroup[[i]]@parnum)
                par <- c(par, tmpgroup[[i]]@par)
                est <- c(est, tmpgroup[[i]]@est)
                lbound <- c(lbound, tmpgroup[[i]]@lbound)
                ubound <- c(ubound, tmpgroup[[i]]@ubound)
                prior.type <- c(prior.type, tmpgroup[[i]]@prior.type)
                prior_1 <- c(prior_1, tmpgroup[[i]]@prior_1)
                prior_2 <- c(prior_2, tmpgroup[[i]]@prior_2)
                item <- c(item, rep('RANDOM', length(tmpgroup[[i]]@est)))
                class <- c(class, rep('RandomPars', length(tmpgroup[[i]]@parnum)))
            }
        }
    }
    if(length(x@lrPars) > 0L){
        lrPars <- x@lrPars
        parname <- c(parname, names(lrPars@est))
        parnum <- c(parnum, lrPars@parnum)
        par <- c(par, lrPars@par)
        est <- c(est, lrPars@est)
        lbound <- c(lbound, lrPars@lbound)
        ubound <- c(ubound, lrPars@ubound)
        tmp <- sapply(as.character(lrPars@prior.type),
                      function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', 'none'))
        prior.type <- c(prior.type, tmp)
        prior_1 <- c(prior_1, lrPars@prior_1)
        prior_2 <- c(prior_2, lrPars@prior_2)
        class <- c(class, rep('lrPars', length(lrPars@parnum)))
        item <- c(item, rep('BETA', length(lrPars@parnum)))
    }
    gnames <- rep(names(PrepList), each = length(est)/length(PrepList))
    par[parname %in% c('g', 'u')] <- antilogit(par[parname %in% c('g', 'u')])
    prior.type <- sapply(as.character(prior.type),
                         function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', 'none'))
    ret <- data.frame(group=gnames, item=item, class=class, name=parname, parnum=parnum, value=par,
                      lbound=lbound, ubound=ubound, est=est, prior.type=prior.type,
                      prior_1=prior_1, prior_2=prior_2)
    ret
}
