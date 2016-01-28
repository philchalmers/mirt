#' Lagrange test for freeing parameters
#'
#' Lagrange (i.e., score) test to test whether parameters should be freed from a
#' more constrained baseline model. The required derivative terms are evaluated numerically
#' from the \code{\link{numerical_deriv}} function. Defining a \code{\link{mirtCluster}} will allow the
#' independent tests to be run in parallel.
#'
#' @param mod an estimated model
#' @param parnum a vector, or list of vectors, containing one or more parameter
#'   locations/sets of locations to be tested.
#'   See objects returned from \code{\link{mod2values}} for the locations
#' @param type type of numerical derivatives to use to find the associated gradient terms. Default is
#'  'forward', but can also be 'Richardson' and 'central' (see \code{\link{numerical_deriv}})
#' @param ... additional arguments to pass to \code{\link{numerical_deriv}}
#'
#' @seealso \code{\link{Wald}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords Lagrange test
#' @export lagrange
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1, 'Rasch')
#' (values <- mod2values(mod))
#'
#' #test all slopes individually
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
#' # test slopes first two slopes and last three slopes jointly
#' lagrange(mod, list(parnum[1:2], parnum[3:5]))
#'
#' # DIF test
#' set.seed(1234)
#' n <- 30
#' N <- 500
#'
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('Group_1', N), rep('Group_2', N))
#'
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group,
#'                      invariance=c('free_means', 'free_var', colnames(dat)))
#' coef(mod, simplify=TRUE)
#' values <- mod2values(mod)
#'
#' # mirtCluster()
#' lagrange(mod, list(c(123, 124), c(239,240)))
#'
#' }
lagrange <- function(mod, parnum, type = 'forward', ...){
    fn <- function(par, mod, pn, dat, model, parprior, PrepList, large, sv, ObJeCtIvE, MG, group, ...){
        sv2 <- sv
        sv2$value[pn] <- par
        if(MG){
            mod2 <- mirt::multipleGroup(dat, model, group=group,
                                        pars = sv2, verbose = FALSE, parprior=parprior, PrepList=PrepList,
                                        large=large, TOL=NaN, calcNull=FALSE,
                                        technical=list(message=FALSE, warn=FALSE,
                                                       parallel=FALSE), ...)

        } else {
            mod2 <- mirt::mirt(dat, model, pars = sv2, verbose = FALSE, parprior=parprior, PrepList=PrepList,
                               large=large, TOL=NaN, calcNull=FALSE, technical=list(message=FALSE, warn=FALSE,
                                                                       parallel=FALSE), ...)
        }
        ret <- extract.mirt(mod2, 'logLik')
        ret
    }
    grads <- function(mod, pn, dat, model, parprior, PrepList, large, sv, ObJeCtIvE, type, MG, group, ...){
        par <- sv$value[pn]
        g <- numerical_deriv(par, fn, mod=mod, pn=pn, dat=dat, model=model, parprior=parprior,
                             PrepList=PrepList, large=large, sv=sv, MG=MG, group=group,
                             ObJeCtIvE=ObJeCtIvE, type=type, ...)
        g
    }
    hess <- function(mod, pn, dat, model, parprior, PrepList, large, sv, ObJeCtIvE, type, MG, group, ...){
        par <- sv$value[pn]
        h <- numerical_deriv(par, fn, mod=mod, pn=pn, dat=dat, model=model, parprior=parprior, group=group,
                             PrepList=PrepList, large=large, sv=sv, MG=MG, ObJeCtIvE=ObJeCtIvE, type=type,
                             gradient=FALSE, ...)
        h
    }

    if(missing(mod)) missingMsg('mod')
    if(missing(parnum)) missingMsg('parnum')
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
    X2 <- mySapply(1:length(parnum), function(i, mod, dat, model, parprior, MG, group,
                                              PrepList, large, sv, ObJeCtIvE, type, ...){
        pn <- parnum[[i]]
        g <- grads(mod=mod, pn=pn, dat=dat, model=model, parprior=parprior, MG=MG, group=group,
                   PrepList=PrepList, large=large, sv=sv, ObJeCtIvE=ObJeCtIvE, type=type, ...)
        h <- -solve(hess(mod=mod, pn=pn, dat=dat, model=model, parprior=parprior, MG=MG, group=group,
                         PrepList=PrepList, large=large, sv=sv, ObJeCtIvE=ObJeCtIvE, type=type, ...))
        as.numeric(g %*% h %*% g)
    }, mod=mod, dat=dat, model=model, parprior=parprior, MG=MG, group=group,
    PrepList=PrepList, large=large, sv=sv, ObJeCtIvE=ObJeCtIvE, type=type, ...)
    ret <- data.frame(X2=as.numeric(X2), df=df, p=sapply(1:length(parnum),
                                             function(ind, X2, df) 1 - pchisq(X2[ind], df[ind]),
                                             X2=X2, df=df))
    rownames(ret) <- sapply(parnum, function(x) paste0(x, collapse = '.'))
    ret
}
