#' Translate mirt parameters into suitable structure for plink package
#'
#' This function exports item parameters from the \code{mirt} package to the
#' \code{plink} package.
#'
#' @aliases read.mirt
#' @param x a single object (or list of objects) returned from \code{mirt, bfactor}, or a
#'   single object returned by \code{multipleGroup}
#' @param as.irt.pars if \code{TRUE}, the parameters will be output as an \code{irt.pars} object
#' @param ... additional arguments to be passed to \code{coef()}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords plink
#' @export read.mirt
#' @examples
#'
#' \donttest{
#'
#' ## unidimensional
#' library(plink)
#'
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, 1))
#' plinkpars <- read.mirt(mod1)
#' plot(plinkpars)
#' plot(mod1, type = 'trace')
#'
#' # graded
#' mod2 <- mirt(Science, 1)
#' plinkpars <- read.mirt(mod2)
#' plot(plinkpars)
#' plot(mod2, type = 'trace')
#'
#' # gpcm
#' mod3 <- mirt(Science, 1, itemtype = 'gpcm')
#' plinkpars <- read.mirt(mod3)
#' plot(plinkpars)
#' plot(mod3, type = 'trace')
#'
#' # nominal
#' mod4 <- mirt(Science, 1, itemtype = 'nominal')
#' plinkpars <- read.mirt(mod4)
#' plot(plinkpars)
#' plot(mod4, type = 'trace')
#'
#' ## multidimensional
#'
#' data <- expand.table(LSAT7)
#' (mod1 <- mirt(data, 2))
#' plinkpars <- read.mirt(mod1)
#' plinkpars
#' plot(plinkpars)
#' plot(mod1, type = 'trace')
#'
#' cmod <- mirt.model('
#'    F1 = 1,4,5
#'    F2 = 2-4')
#' model <- mirt(data, cmod)
#' plot(read.mirt(model))
#' itemplot(model, 1)
#'
#' # graded
#' mod2 <- mirt(Science, 2)
#' plinkpars <- read.mirt(mod2)
#' plinkpars
#' plot(plinkpars)
#' plot(mod2, type = 'trace')
#'
#' ### multiple group equating example
#' set.seed(1234)
#' dat <- expand.table(LSAT7)
#' group <- sample(c('g1', 'g2'), nrow(dat), TRUE)
#' dat1 <- dat[group == 'g1', ]
#' dat2 <- dat[group == 'g2', ]
#' mod1 <- mirt(dat1, 1)
#' mod2 <- mirt(dat2, 1)
#'
#' # convert and combine pars
#' plinkMG <- read.mirt(list(g1=mod1, g2=mod2))
#'
#' # equivalently:
#' # mod <- multipleGroup(dat, 1, group)
#' # plinkMG <- read.mirt(mod)
#'
#' combine <- matrix(1:5, 5, 2)
#' comb <- combine.pars(plinkMG, combine, grp.names=unique(group))
#' out <- plink(comb, rescale="SL")
#' equate(out)
#' equate(out, method = 'OSE')
#'
# # gpcm
# mod3 <- mirt(Science, 2, itemtype = 'gpcm')
# plinkpars <- read.mirt(mod3)
# plot(plinkpars)
# itemplot(mod3, 1)
#
# # nominal
# mod4 <- mirt(Science, 2, itemtype = 'nominal')
# plinkpars <- read.mirt(mod4)
# plot(plinkpars)
# itemplot(mod4, 1)
#' }
read.mirt <- function (x, as.irt.pars = TRUE, ...)
{
    if(requireNamespace("plink", quietly = TRUE)){
        if(is(x, 'MultipleGroupClass')){
            pars <- vector('list', length(extract.mirt(x, 'pars')))
            for(i in 1:length(pars)){
                tmp <- extract.group(x, group=i)
                pars[[i]] <- read.mirt(tmp, as.irt.pars=as.irt.pars, ...)
            }
            names(pars) <- extract.mirt(x, 'groupNames')
            return(pars)
        } else if(is.list(x)){
            pars <- vector('list', length(x))
            names(pars) <- names(x)
            for(i in 1:length(pars))
                pars[[i]] <- read.mirt(x[[i]], as.irt.pars=as.irt.pars, ...)
            return(pars)
        }
        if(is(x, 'MixedClass'))
            stop('Mixed effect models not supported.', call.=FALSE)
        if(length(extract.mirt(x, 'prodlist')))
            stop('Polynomial factor models not supported in plink', call.=FALSE)
        #converts unidimensional parameters to classic IRT (if possible)
        nfact <- extract.mirt(x, 'nfact')
        nitems <- extract.mirt(x, 'nitems')
        listpars <- coef(x, IRTpars=ifelse(nfact == 1, TRUE, FALSE), rotate='none', verbose=FALSE, ...)
        if(!is(listpars[[1]], 'matrix'))
            for(i in 1:nitems)
                listpars[[i]] <- t(matrix(listpars[[i]]))
        mirt.items <- as.character(lapply(extract.mirt(x, 'pars'), class))
        mirt.items <- plink.items <- mirt.items[-(nitems+1)]
        cat <- numeric(nitems)
        pars <- matrix(NA, nitems, 40)
        theta <- -4:4
        Theta <- thetaComb(theta, nfact)
        K <- extract.mirt(x, 'K')
        for(i in 1:nitems){
            if(mirt.items[i] == 'dich'){
                plink.items[i] <- 'drm'
                cat[i] <- 2
                abc <- listpars[[i]][1, 1:(nfact+2)]
                pars[i, 1:length(abc)] <- abc
                next
            }
            if(mirt.items[i] == 'graded'){
                plink.items[i] <- 'grm'
                ab <- listpars[[i]][1, ]
                cat[i] <- K[i]
                pars[i, 1:length(ab)] <- ab
                next
            }
            if(mirt.items[i] == 'rsm'){
                stop('Rasch rating scale models not supported for now', call.=FALSE)
            }
            if(mirt.items[i] == 'nestlogit'){
                stop('nestlogit models not supported in plink', call.=FALSE)
            }
            if(mirt.items[i] == 'rating'){
                stop('rating model not supported for now', call.=FALSE)
                #not converted to classic IRT form for now
                plink.items[i] <- 'grm'
                ab <- listpars[[i]][1, ]
                adj <- ab[length(ab)]
                ab <- ab[-length(ab)]
                cat[i] <- K[i]
                for(j in 1:(cat[i] - 1))
                    ab[j + nfact] <- ab[j+1] + adj
                pars[i, 1:length(ab)] <- ab
                next
            }
            if(mirt.items[i] == 'gpcm'){
                ab <- listpars[[i]][1, ]
                a <- ab[1:nfact]
                ab <- ab[-(1:nfact)]
                cat[i] <- K[i]
                if(nfact == 1L){
                    pars[i, 1:length(ab)] <- ab
                } else {
                    stop('Multidimensional gpcm not yet supported', call.=FALSE)
                }
                next
            }
            if(mirt.items[i] == 'nominal'){
                plink.items[i] <- 'nrm'
                ab <- listpars[[i]][1, ]
                cat[i] <- K[i]
                if(nfact == 1L){
                    pars[i, 1:length(ab)] <- ab
                } else {
                    stop('Multidimensional nrm not yet supported', call.=FALSE)
                }
                next
            }
            if(mirt.items[i] == 'partcomp'){
                stop('Partially compensatory models not supported in plink', call.=FALSE)
            }
            if(mirt.items[i] == 'custom'){
                stop('User defined models not supported in plink', call.=FALSE)
            }
        }
        model <- unique(plink.items)
        items <- vector('list', length(model))
        index <- 1:nitems
        for(i in 1:length(model))
            items[[i]] <- index[model[i] == plink.items]
        pm <- plink::as.poly.mod(nitems, model=model, items=items)
        pars <- pars[ , colSums(is.na(pars)) != nitems]
        if(as.irt.pars)
            pars <- plink::as.irt.pars(pars, cat=cat, poly.mod=pm, dimensions=nfact, ability=Theta)
        return(pars)
    } else {
        stop('plink package is not available. Please install.', call.=FALSE)
    }
}