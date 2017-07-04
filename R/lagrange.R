#' Lagrange test for freeing parameters
#'
#' Lagrange (i.e., score) test to test whether parameters should be freed from a
#' more constrained baseline model.
#'
#' @param mod an estimated model
#' @param parnum a vector, or list of vectors, containing one or more parameter
#'   locations/sets of locations to be tested.
#'   See objects returned from \code{\link{mod2values}} for the locations
#' @param SE.type type of information matrix estimator to use. See \code{\link{mirt}} for
#'   further details
#' @param type type of numerical algorithm passed to \code{\link{numerical_deriv}} to
#'   obtain the gradient terms
#' @param ... additional arguments to pass to \code{\link{mirt}}
#'
#' @seealso \code{\link{wald}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords Lagrange test
#' @export lagrange
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1, 'Rasch')
#' (values <- mod2values(mod))
#'
#' #test all fixed slopes individually
#' parnum <- values$parnum[values$name == 'a1']
#' lagrange(mod, parnum)
#'
#' # compare to LR test for first two slopes
#' mod2 <- mirt(dat, 'F = 1-5
#'                    FREE = (1, a1)', 'Rasch')
#' coef(mod2, simplify=TRUE)$items
#' anova(mod, mod2)
#'
#' mod2 <- mirt(dat, 'F = 1-5
#'                    FREE = (2, a1)', 'Rasch')
#' coef(mod2, simplify=TRUE)$items
#' anova(mod, mod2)
#'
#' mod2 <- mirt(dat, 'F = 1-5
#'                    FREE = (3, a1)', 'Rasch')
#' coef(mod2, simplify=TRUE)$items
#' anova(mod, mod2)
#'
#' # test slopes first two slopes and last three slopes jointly
#' lagrange(mod, list(parnum[1:2], parnum[3:5]))
#'
#' # test all 5 slopes and first + last jointly
#' lagrange(mod, list(parnum[1:5], parnum[c(1, 5)]))
#'
#' }
lagrange <- function(mod, parnum, SE.type = 'Oakes', type = 'central', ...){
    fn <- function(par, mod, pn, dat, model, parprior, PrepList, large, sv, ObJeCtIvE, MG, group, ...){
        sv2 <- sv
        sv2$value[pn] <- par
        if(MG){
            mod2 <- mirt::multipleGroup(dat, model, group=group,
                                        pars = sv2, verbose = FALSE, parprior=parprior, PrepList=PrepList,
                                        large=large, TOL=NaN, calcNull=FALSE,
                                        technical=list(message=FALSE, warn=FALSE,
                                                       parallel=FALSE, Etable=FALSE), ...)

        } else {
            mod2 <- mirt::mirt(dat, model, pars = sv2, verbose = FALSE, parprior=parprior, PrepList=PrepList,
                               large=large, TOL=NaN, calcNull=FALSE, technical=list(message=FALSE, warn=FALSE,
                                                                       parallel=FALSE, Etable=FALSE), ...)
        }
        ret <- extract.mirt(mod2, 'logLik')
        ret
    }
    grads <- function(mod, pn, dat, model, parprior, PrepList, large, sv, ObJeCtIvE, type = 'central', MG, group, ...){
        par <- sv$value[pn]
        g <- numerical_deriv(par, fn, mod=mod, pn=pn, dat=dat, model=model, parprior=parprior,
                             PrepList=PrepList, large=large, sv=sv, MG=MG, group=group,
                             ObJeCtIvE=ObJeCtIvE, type=type, ...)
        g
    }

    if(missing(mod)) missingMsg('mod')
    if(missing(parnum)) missingMsg('parnum')
    if(SE.type %in% c('SEM', 'MHRM', 'FMHRM'))
        stop('SE.type not supported for Lagrange tests')
    ObJeCtIvE <- extract.mirt(mod, 'logLik')
    group <- extract.mirt(mod, 'group')
    parnum <- as.list(parnum)
    df <- sapply(parnum, length)
    MG <- ifelse(class(mod) == 'MultipleGroupClass', TRUE, FALSE)
    dat <- extract.mirt(mod, 'data')
    large <- if(MG) mirt::multipleGroup(dat, 1, group=group, large=TRUE)
       else mirt::mirt(dat, 1, large=TRUE)
    model <- extract.mirt(mod, "model")
    parprior <- extract.mirt(mod, "parprior")
    PrepList <- mirt(mod@Data$data, mod@Model$model, Return_PrepList=TRUE)
    sv <- mod2values(mod)
    sv2 <- sv
    sv2$est[do.call(c, parnum)] <- TRUE
    nms <- vector('list', length(parnum))
    for(i in 1:length(parnum)){
        nms[[i]] <- apply(cbind(as.character(sv2$name[parnum[[i]]]), '.', parnum[[i]]),
                            1, paste0, collapse='')
    }
    modInfo <- if(MG) mirt::multipleGroup(dat, extract.mirt(mod, 'model'), TOL = NaN, sv=sv2,
                                        group=group, SE=TRUE, SE.type = SE.type, large=large,
                                        technical = list(infoAsVcov=TRUE, warn=FALSE), verbose=FALSE)
    else mirt::mirt(dat, extract.mirt(mod, 'model'), SE=TRUE, SE.type = SE.type, sv=sv2,
                    large=large, TOL = NaN, technical = list(infoAsVcov=TRUE, warn=FALSE), verbose=FALSE)
    info <- vcov(modInfo)
    X2 <- mySapply(1:length(parnum), function(i, mod, dat, model, parprior, MG, group,
                                              PrepList, large, sv, ObJeCtIvE, info, nms, type, ...){
        pn <- parnum[[i]]
        g <- grads(mod=mod, pn=pn, dat=dat, model=model, parprior=parprior, MG=MG, group=group,
                   PrepList=PrepList, large=large, sv=sv, ObJeCtIvE=ObJeCtIvE, type=type, ...)
        pick <- colnames(info)
        tmp <- do.call(c, nms[-i])
        tmp <- tmp[!(tmp %in% nms[[i]])]
        pick <- pick[!(pick %in% tmp)]
        vcov <- solve(info[pick, pick])
        h <- vcov[colnames(vcov) %in% nms[[i]], colnames(vcov) %in% nms[[i]], drop=FALSE]
        as.numeric(g %*% h %*% g)
    }, mod=mod, dat=dat, model=model, parprior=parprior, MG=MG, group=group, type=type,
    PrepList=PrepList, large=large, sv=sv, ObJeCtIvE=ObJeCtIvE, info=info, nms=nms, ...)
    ret <- data.frame(X2=as.numeric(X2), df=df, p=sapply(1:length(parnum),
                                             function(ind, X2, df) 1 - pchisq(X2[ind], df[ind]),
                                             X2=X2, df=df))
    rownames(ret) <- sapply(parnum, function(x) paste0(x, collapse = '.'))
    ret
}
