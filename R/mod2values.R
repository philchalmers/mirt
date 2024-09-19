#' Convert an estimated mirt model to a data.frame
#'
#' Given an estimated model from any of mirt's model fitting functions this function will convert
#' the model parameters into the design data frame of starting values and other parameter
#' characteristics (similar to using the \code{pars = 'values'} for obtaining starting values).
#'
#'
#' @aliases mod2values
#' @param x an estimated model x from the mirt package
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords convert model
#' @export mod2values
#' @seealso \code{\link{extract.mirt}}
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, "F=1-5
#'                   CONSTRAIN=(1-5, a1)")
#' values <- mod2values(mod)
#' values
#'
#' # use the converted values as starting values in a new model, and reduce TOL
#' mod2 <- mirt(dat, 1, pars = values, TOL=1e-5)
#' coef(mod2, simplify=TRUE)
#'
#' # use parameters on different dataset
#' mod3 <- mirt(expand.table(LSAT6), pars=values)
#' coef(mod3, simplify=TRUE)
#'
#' # supports differing itemtypes on second model
#' sv <- mirt(Science, itemtype=c('graded', rep('gpcm', 3)), pars='values')
#' mod3 <- mirt(Science, pars = sv)  # itemtype omitted
#' coef(mod3, simplify=TRUE)$items
#' extract.mirt(mod3, 'itemtype')
#'
#'
#' }
mod2values <- function(x){
    if(is(x, 'MultipleGroupClass') || is(x, 'DiscreteClass') || is(x, 'MixtureClass')){
        PrepList <- x@ParObjects$pars
        names(PrepList) <- x@Data$groupNames
        MG <- TRUE
    } else {
        PrepList <- list(pars=x@ParObjects$pars)
        names(PrepList) <- 'all'
        MG <- FALSE
    }
    itemnames <- colnames(x@Data$data)
    parnum <- par <- est <- item <- parname <- gnames <- class <-
        lbound <- ubound <- prior.type <- prior_1 <- prior_2 <- c()
    for(g in 1L:length(PrepList)){
        if(MG) tmpgroup <- PrepList[[g]]@ParObjects$pars
        else tmpgroup <- PrepList[[g]]
        for(i in 1L:length(tmpgroup)){
            if(i <= length(itemnames))
                item <- c(item, rep(itemnames[i], length(tmpgroup[[i]]@parnum)))
            class <- c(class, rep(class(tmpgroup[[i]]), length(tmpgroup[[i]]@parnum)))
            parname <- c(parname, tmpgroup[[i]]@parnames)
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
        tmpgroup <- x@ParObjects$random
        if(length(tmpgroup)){
            for(i in 1L:length(tmpgroup)){
                parname <- c(parname, tmpgroup[[i]]@parnames)
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
    if(length(x@Model$lrPars) > 0L){
        lrPars <- x@Model$lrPars
        parname <- c(parname, lrPars@parnames)
        parnum <- c(parnum, lrPars@parnum)
        par <- c(par, lrPars@par)
        est <- c(est, lrPars@est)
        lbound <- c(lbound, lrPars@lbound)
        ubound <- c(ubound, lrPars@ubound)
        tmp <- sapply(as.character(lrPars@prior.type),
                      function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', '4'='expbeta', 'none'))
        prior.type <- c(prior.type, tmp)
        prior_1 <- c(prior_1, lrPars@prior_1)
        prior_2 <- c(prior_2, lrPars@prior_2)
        class <- c(class, rep('lrPars', length(lrPars@parnum)))
        item <- c(item, rep('BETA', length(lrPars@parnum)))
    }
    if(is(x, 'MixedClass')){
        tmpgroup <- x@ParObjects$lr.random
        if(length(tmpgroup)){
            for(i in 1L:length(tmpgroup)){
                parname <- c(parname, tmpgroup[[i]]@parnames)
                parnum <- c(parnum, tmpgroup[[i]]@parnum)
                par <- c(par, tmpgroup[[i]]@par)
                est <- c(est, tmpgroup[[i]]@est)
                lbound <- c(lbound, tmpgroup[[i]]@lbound)
                ubound <- c(ubound, tmpgroup[[i]]@ubound)
                prior.type <- c(prior.type, tmpgroup[[i]]@prior.type)
                prior_1 <- c(prior_1, tmpgroup[[i]]@prior_1)
                prior_2 <- c(prior_2, tmpgroup[[i]]@prior_2)
                item <- c(item, rep('LRRANDOM', length(tmpgroup[[i]]@est)))
                class <- c(class, rep('LRRandomPars', length(tmpgroup[[i]]@parnum)))
            }
        }
    }
    gnames <- rep(names(PrepList), each = length(est)/length(PrepList))
    par[parname %in% c('g', 'u')] <- antilogit(par[parname %in% c('g', 'u')])
    lbound[parname %in% c('g', 'u')] <- antilogit(lbound[parname %in% c('g', 'u')])
    ubound[parname %in% c('g', 'u')] <- antilogit(ubound[parname %in% c('g', 'u')])
    prior.type <- sapply(as.character(prior.type),
                         function(x) switch(x, '1'='norm', '2'='lnorm', '3'='beta', '4'='expbeta', 'none'))
    clist <- extract.mirt(x, 'constrain')
    constrain <- nconstrain <- rep("none", length(gnames))
    if(!is.null(clist)){
        for(i in seq_len(length(clist)))
            constrain[clist[[i]]] <- i
    }
    nclist <- extract.mirt(x, 'nconstrain')
    if(!is.null(nclist)){
        for(i in seq_len(length(nclist)))
            nconstrain[nclist[[i]]] <- i
    }
    ret <- data.frame(group=gnames, item=item, class=class, name=parname, parnum=parnum,
                      value=par, lbound=lbound, ubound=ubound, est=est, const=constrain, nconst=nconstrain,
                      prior.type=prior.type, prior_1=prior_1, prior_2=prior_2, stringsAsFactors = FALSE)
    ret <- as.mirt_df(ret)
    attr(ret, 'itemtype') <- extract.mirt(x, 'itemtype')
    ret
}
