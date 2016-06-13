#' Compute profiled-likelihood (or posterior) confidence intervals
#'
#' Computes profiled-likelihood based confidence intervals. Supports the inclusion of
#' equality constraints. Object returns the confidence intervals
#' and whether the respective interval could be found.
#'
#' @aliases PLCI.mirt
#' @param mod a converged mirt model
#' @param alpha two-tailed alpha critical level
#' @param parnum a numeric vector indicating which parameters to estimate.
#'   Use \code{\link{mod2values}} to determine parameter numbers. If \code{NULL}, all possible
#'   parameters are used
#' @param fun a function used to define a transformation of 1 or more parameters of the
#'   form \code{fun(par)}, where \code{par} is one or more parameters indicated by
#'   the \code{parnum} input. This function must return a single numeric value
#' @param inv.fun the inverse transformation function of \code{fun} of the form
#'   \code{inv.fun(newpar, par)}, where \code{par} is defined as before, however
#'   \code{newpar} is the the transformed parameter(s) (e.g., returned from \code{fun(par)}).
#'   Note that the length of the returned vector should be the same as the original length
#'   of \code{par}
#' @param fun.lower lower bound used when \code{fun} is supplied. Default is -Inf
#' @param fun.upper upper bound used when \code{fun} is supplied. Default is Inf
#' @param ... additional arguments to pass to the estimation functions
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords profiled likelihood
#' @export PLCI.mirt
#' @seealso
#' \code{\link{boot.mirt}}
#'
#' @examples
#'
#' \dontrun{
#' mirtCluster() #use all available cores to estimate CI's in parallel
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#'
#' result <- PLCI.mirt(mod)
#' result
#'
#' mod2 <- mirt(Science, 1)
#' result2 <- PLCI.mirt(mod2)
#' result2
#'
#' #only estimate CI's slopes
#' sv <- mod2values(mod2)
#' parnum <- sv$parnum[sv$name == 'a1']
#' result3 <- PLCI.mirt(mod2, parnum=parnum)
#' result3
#'
#' #####
#' # transformation example for computing PL CI of difficulty parameter
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#'
#' fun <- function(par) -par[2]/par[1]
#' inv.fun <- function(newpar, par) c(par[1], -newpar * par[1])
#'
#' PLCI.mirt(mod, parnum=c(1,2), fun=fun, inv.fun=inv.fun)
#'
#' }
PLCI.mirt <- function(mod, alpha = .05, parnum = NULL, fun = NULL, inv.fun = NULL,
                      fun.lower = -Inf, fun.upper = Inf, ...){

    #silently accepts print_debug = TRUE for printing the minimization criteria

    compute.LL <- function(dat, model, sv, large, parprior, PrepList, technical, ...){
        if(missing(technical))
            technical <- list(message=FALSE, warn=FALSE, parallel=FALSE)
        else {
            technical$message <- technical$warn <- technical$parallel <- FALSE
        }
        tmpmod <- mirt::mirt(dat, model, pars = sv, verbose = FALSE, parprior=parprior, PrepList=PrepList,
                                        large=large, calcNull=FALSE, technical=technical, ...)
        coef(tmpmod, simplify=TRUE)
        ret <- list(LL=tmpmod@Fit$logLik + tmpmod@Fit$logPrior, vals=mod2values(tmpmod))
        ret
    }

    f.min <- function(value, dat, model, which, sv, get.LL, large, parprior, parnames, asigns,
                      PrepList, fun, inv.fun, print_debug = FALSE, ...){
        if(!is.null(fun)){
            sv$est[parnums] <- FALSE
            sv$value[parnums] <- inv.fun(value, sv$value[parnums])
        } else {
            sv$est[which] <- FALSE
            sv$value[which] <- value
            if(sv$class[which] == 'graded'){
                if(!(sv$name[which] %in% paste0('a', 1L:30L))){
                    itemname <- sv$item[which]
                    itemnum <- sv$parnum[which]
                    itemsv <- sv[sv$item == itemname & !(sv$name %in% paste0('a', 1L:30L)), ]
                    ds <- dsnew <- itemsv$value
                    srtds <- sort(ds, decreasing = TRUE)
                    if(!all(ds == srtds)){
                        est <- itemsv$est
                        if(est[1]) dsnew[1L] <- max(ds) + 1
                        if(est[length(est)]) dsnew[length(est)] <- min(ds) - 1
                        seqd <- seq(from=dsnew[1L], to=dsnew[length(est)], length.out = length(ds))
                        ds[est] <- seqd[est]
                        sv$value[itemsv$parnum] <- ds
                    }
                }
            }
        }
        got.LL <- try(compute.LL(dat=dat, model=model, sv=sv, large=large, parprior=parprior,
                                 PrepList=PrepList, ...), silent=TRUE)
        sv2 <- got.LL$vals
        got.LL <- got.LL$LL
        as <- matrix(sv2$value[sv2$name %in% paste0('a', 1L:30L)], ncol(dat))
        if(sum(asigns * sign(as)) < 0L) return(-1e10)
        ret <- (got.LL - get.LL)
        if(print_debug) cat('parnum = ', which, '; value = ', round(value, 3),
                            '; min = ', round(ret, 3), '\n')
        attr(ret, 'value') <- value
        ret
    }

    LLpar <- function(parnum, parnums, parnames, lbound, ubound, dat, model, large,
                      sv, get.LL, parprior, asigns, PrepList, pars, maxLL,
                      fun, inv.fun, ...){
        TOL <- .001
        if(is.null(fun)){
            lower <- ifelse(lbound[parnum] == -Inf, -15, lbound[parnum])
            upper <- ifelse(ubound[parnum] == Inf, 15, ubound[parnum])
            mid <- pars[parnum]
            if(parnames[parnum] %in% c('g', 'u')){
                lower <- 0
                upper <- 1
            } else if(parnames[parnum] %in% paste0('COV_', 1:30, 1:30)){
                lower <- 1e-4
            }
        } else {
            lower <- ifelse(lbound == -Inf, -30, lbound)
            upper <- ifelse(ubound == Inf, 30, ubound)
            mid <- fun(pars[parnums])
        }

        if(mid > lower){
            opt.lower <- try(uniroot(f.min, c(lower, mid), dat=dat, model=model,
                                 large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                 parprior=parprior, parnames=parnames, asigns=asigns,
                                 PrepList=PrepList, fun=fun, inv.fun=inv.fun,
                                 ..., f.upper=maxLL-get.LL, tol = TOL/10),
                             silent = TRUE)
            if(is(opt.lower, 'try-error')) opt.lower <- list(root = lower, f.root=1e10)
        } else opt.lower <- list(root = lower, f.root=1e10)
        if(mid < upper){
            opt.upper <- try(uniroot(f.min, c(mid, upper), dat=dat, model=model,
                                 large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                 parprior=parprior, parnames=parnames, asigns=asigns,
                                 PrepList=PrepList, fun=fun, inv.fun=inv.fun,
                                 ..., f.lower=maxLL-get.LL, tol = TOL/10),
                             silent = TRUE)
            if(is(opt.upper, 'try-error')) opt.upper <- list(root = upper, f.root=1e10)
        } else opt.upper <- list(root = upper, f.root=1e10)
        conv_upper <- conv_lower <- TRUE
        if(abs(opt.lower$f.root) > TOL)
            conv_lower <- FALSE
        if(abs(opt.upper$f.root) > TOL)
            conv_upper <- FALSE
        c(lower=opt.lower$root, upper=opt.upper$root,
          conv_lower=conv_lower, conv_upper=conv_upper)
    }

    if(.hasSlot(mod@Model$lrPars, 'beta'))
        stop('Latent regression models not yet supported')
    dat <- mod@Data$data
    model <- mod@Model$model
    parprior <- mod@Model$parprior
    if(length(parprior))
        stop('Confidence intervals cannot be computed for models that include priors')
    if(length(parprior) == 0L) parprior <- NULL
    sv <- mod2values(mod)
    PrepList <- mirt(mod@Data$data, mod@Model$model, Return_PrepList=TRUE)
    large <- mirt(mod@Data$data, mod@Model$model, large = TRUE)
    as <- matrix(sv$value[sv$name %in% paste0('a', 1L:30L)], ncol(dat))
    asigns <- sign(as)
    if(!is.null(fun) && is.null(parnum))
        stop('When supplying a tranformation function you must supply parnum')
    if(!is.null(parnum)){
        tmp <- sv$parnum %in% parnum
        pars <- sv$value[tmp]
        parnums <- sv$parnum[tmp]
        itemtypes <- sv$class[tmp]
        if(!is.null(fun)){
            lbound <- fun.lower
            ubound <- fun.upper
            which_parnums <- 1L
        } else {
            parnames <- sv$name[tmp]
            lbound <- sv$lbound[tmp]
            ubound <- sv$ubound[tmp]
            which_parnums <- 1L:length(parnums)
        }
    } else {
        pars <- sv$value
        parnums <- sv$parnum[sv$est]
        itemtypes <- sv$class[sv$est]
        parnames <- sv$name[sv$est]
        lbound <- sv$lbound[sv$est]
        ubound <- sv$ubound[sv$est]
        which_parnums <- 1L:length(parnums)
    }
    LL <- mod@Fit$logLik
    get.LL <- LL - qchisq(1-alpha, 1)/2
    result <- mySapply(X=which_parnums, FUN=LLpar, pars=pars, parnums=parnums, asigns=asigns,
                       parnames=parnames, lbound=lbound, ubound=ubound, dat=dat,
                       model=model, large=large, sv=sv, get.LL=get.LL, parprior=parprior,
                       PrepList=PrepList, maxLL=LL, fun=fun, inv.fun=inv.fun, ...)
    colnames(result) <- c(paste0('lower_', alpha/2*100), paste0('upper_', (1-alpha/2)*100),
                          'lower_conv', 'upper_conv')
    ret <- if(!is.null(fun)){
        data.frame(parnam='fun(par)', value=fun(pars), result, row.names=NULL)
    } else {
        data.frame(Item=sv$item[parnums], class=itemtypes, parnam=sv$name[parnums],
                   parnum=parnums, value=pars, result, row.names=NULL)
    }
    ret$lower_conv <- as.logical(ret$lower_conv)
    ret$upper_conv <- as.logical(ret$upper_conv)
    ret
}
