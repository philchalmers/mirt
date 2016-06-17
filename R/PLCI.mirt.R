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
#' @param inf2val a numeric used to change parameter bounds which are infinity to a finite number.
#'   Decreasing this too much may not allow a suitable bound to be located. Default is 30
#' @param search_bound logical; use a fixed grid of values around the ML estimate to
#'   determine more suitable optimization bounds? Using this has much better behaviour
#'   than setting fixed upper/lower bound values and searching from more extreme ends
#' @param step magnitude of steps used when \code{search_bound} is \code{TRUE}.
#'   Smaller values create more points to search a suitable bound for (up to the
#'   lower bound value visible with \code{\link{mod2values}})
#' @param lower logical; search for the lower CI?
#' @param upper logical; search for the upper CI?
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
#' }
PLCI.mirt <- function(mod, alpha = .05, parnum = NULL,
                      search_bound = TRUE, step = .5,
                      lower = TRUE, upper = TRUE, inf2val = 30, ...){

    #silently accepts print_debug = TRUE for printing the minimization criteria

    compute.LL <- function(dat, model, sv, large, parprior, PrepList, itemtype,
                           technical, ...){
        if(missing(technical))
            technical <- list(message=FALSE, warn=FALSE, parallel=FALSE)
        else {
            technical$message <- technical$warn <- technical$parallel <- FALSE
        }
        tmpmod <- mirt::mirt(dat, model, itemtype=itemtype, pars = sv, verbose = FALSE,
                             parprior=parprior, PrepList=PrepList, large=large, calcNull=FALSE,
                             technical=technical, ...)
        # coef(tmpmod, simplify=TRUE)
        ret <- list(LL=tmpmod@Fit$logLik + tmpmod@Fit$logPrior, vals=mod2values(tmpmod))
        ret
    }

    f.min <- function(value, dat, model, which, sv, get.LL, large, parprior, parnames, asigns,
                      PrepList, itemtype, print_debug = FALSE, ...){
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
        got.LL <- try(compute.LL(dat=dat, model=model, itemtype=itemtype,
                                 sv=sv, large=large, parprior=parprior,
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
                      sv, get.LL, parprior, asigns, PrepList, pars, itemtype, inf2val,
                      maxLL, estlower, estupper, search_bound, step, ...){
        TOL <- .001
        lower <- ifelse(lbound[parnum] == -Inf, -inf2val, lbound[parnum])
        upper <- ifelse(ubound[parnum] == Inf, inf2val, ubound[parnum])
        mid <- pars[parnum]
        if(parnames[parnum] %in% c('g', 'u')){
            lower <- 0
            upper <- 1
        } else if(parnames[parnum] %in% paste0('COV_', 1:30, 1:30)){
            lower <- 1e-4
        }
        if(estlower && mid > lower){
            possible_bound <- TRUE
            if(search_bound){
                grid <- mid - cumsum(rep(step, floor(abs(lower/step))))
                for(g in grid){
                    Xval <- f.min(g, dat=dat, model=model,
                                   large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                   parprior=parprior, parnames=parnames, asigns=asigns,
                                   PrepList=PrepList, itemtype=itemtype, ...)
                    if(abs(Xval) > abs(get.LL - maxLL)){
                        lower <- g
                        break
                    }
                    if(Xval == -1e10){
                        possible_bound <- FALSE
                        break
                    }
                }
                if(g == grid[length(grid)] && abs(Xval) < abs(get.LL - maxLL))
                    possible_bound <- FALSE
            }
            if(possible_bound)
                opt.lower <- try(uniroot(f.min, c(lower, mid), dat=dat, model=model,
                                     large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                     parprior=parprior, parnames=parnames, asigns=asigns,
                                     PrepList=PrepList, itemtype=itemtype,
                                     ..., f.upper=maxLL-get.LL, tol = TOL/10),
                                 silent = TRUE)
            else opt.lower <- try(uniroot(), TRUE)
            if(is(opt.lower, 'try-error')) opt.lower <- list(root = lower, f.root=1e10)
        } else opt.lower <- list(root = lower, f.root=1e10)
        if(estupper && mid < upper){
            possible_bound <- TRUE
            if(search_bound){
                grid <- mid + cumsum(rep(step, floor(abs(upper/step))))
                for(g in grid){
                    Xval <- f.min(g, dat=dat, model=model,
                                   large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                   parprior=parprior, parnames=parnames, asigns=asigns,
                                   PrepList=PrepList, itemtype=itemtype, ...)
                    if(abs(Xval) > abs(get.LL - maxLL)){
                        upper <- g
                        break
                    }
                    if(Xval == -1e10){
                        possible_bound <- FALSE
                        break
                    }
                }
                if(g == grid[length(grid)] && abs(Xval) < abs(get.LL - maxLL))
                    possible_bound <- FALSE
            }
            if(possible_bound)
                opt.upper <- try(uniroot(f.min, c(mid, upper), dat=dat, model=model,
                                     large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                     parprior=parprior, parnames=parnames, asigns=asigns,
                                     PrepList=PrepList, itemtype=itemtype,
                                     ..., f.lower=maxLL-get.LL, tol = TOL/10),
                                 silent = TRUE)
            else opt.upper <- try(uniroot(), TRUE)
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
    stopifnot(lower | upper)
    dat <- mod@Data$data
    model <- mod@Model$model
    parprior <- mod@Model$parprior
    if(length(parprior))
        stop('Confidence intervals cannot be computed for models that include priors')
    if(length(parprior) == 0L) parprior <- NULL
    sv <- mod2values(mod)
    itemtype <- extract.mirt(mod, 'itemtype')
    PrepList <- mirt(mod@Data$data, mod@Model$model, itemtype=itemtype,
                     Return_PrepList=TRUE, ...)
    large <- mirt(mod@Data$data, mod@Model$model, large = TRUE)
    as <- matrix(sv$value[sv$name %in% paste0('a', 1L:30L)], ncol(dat))
    asigns <- sign(as)
    if(!is.null(parnum)){
        tmp <- sv$parnum %in% parnum
        pars <- sv$value[tmp]
        LL.upper.crit <- LL.lower.crit <- pars
        parnums <- sv$parnum[tmp]
        itemtypes <- sv$class[tmp]
        parnames <- sv$name[tmp]
        lbound <- sv$lbound[tmp]
        ubound <- sv$ubound[tmp]
    } else {
        pars <- sv$value
        pars <- LL.upper.crit <- LL.lower.crit <- pars[sv$est]
        parnums <- sv$parnum[sv$est]
        itemtypes <- sv$class[sv$est]
        parnames <- sv$name[sv$est]
        lbound <- sv$lbound[sv$est]
        ubound <- sv$ubound[sv$est]
    }
    LL <- mod@Fit$logLik
    get.LL <- LL - qchisq(1-alpha, 1)/2
    result <- mySapply(X=1L:length(parnums), FUN=LLpar, pars=pars, parnums=parnums, asigns=asigns,
                       parnames=parnames, lbound=lbound, ubound=ubound, dat=dat,
                       model=model, large=large, sv=sv, get.LL=get.LL, parprior=parprior,
                       PrepList=PrepList, itemtype=itemtype, inf2val=inf2val, maxLL=LL,
                       estlower=lower, estupper=upper, search_bound=search_bound,
                       step=step, ...)
    colnames(result) <- c(paste0('lower_', alpha/2*100), paste0('upper_', (1-alpha/2)*100),
                                 'lower_conv', 'upper_conv')
    ret <- data.frame(Item=sv$item[parnums], class=itemtypes, parnam=sv$name[parnums],
                      parnum=parnums, value=pars, result, row.names=NULL)
    ret$lower_conv <- as.logical(ret$lower_conv)
    ret$upper_conv <- as.logical(ret$upper_conv)
    if(!lower) ret <- ret[,!grepl('lower', colnames(ret)), drop=FALSE]
    if(!upper) ret <- ret[,!grepl('upper', colnames(ret)), drop=FALSE]
    ret
}
