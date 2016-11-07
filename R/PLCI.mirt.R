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
#'   lower bound value visible with \code{\link{mod2values}}). When upper/lower bounds are detected
#'   this value will be adjusted accordingly
#' @param lower logical; search for the lower CI?
#' @param upper logical; search for the upper CI?
#' @param NealeMiller logical; use the Neale and Miller 1997 approximation? Default is \code{FALSE}
#' @param ... additional arguments to pass to the estimation functions
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords profiled likelihood
#' @export PLCI.mirt
#' @seealso
#' \code{\link{boot.mirt}}
#'
#' @references
#'
#' Neale, M. C. & Miller, M. B. (1997). The use of likelihood-based confidence intervals in genetic
#' models. Behavior Genetircs, 27, 113-120
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
#' # model with constraints
#' mod <- mirt(dat, 'F = 1-5
#'                   CONSTRAIN = (1-5, a1)')
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
#' result3 <- PLCI.mirt(mod2, parnum)
#' result3
#'
#' }
PLCI.mirt <- function(mod, parnum = NULL, alpha = .05,
                      search_bound = TRUE, step = .5,
                      lower = TRUE, upper = TRUE, inf2val = 30,
                      NealeMiller = FALSE, ...){

    #silently accepts print_debug = TRUE for printing the minimization criteria

    compute.LL <- function(dat, model, sv, large, parprior, PrepList, itemtype, constrain,
                           technical, ...){
        if(missing(technical))
            technical <- list(message=FALSE, warn=FALSE, parallel=FALSE, PLCI=TRUE)
        else {
            technical$message <- technical$warn <- technical$parallel <- FALSE
            technical$PLCI <- TRUE
        }
        tmpmod <- mirt::mirt(dat, model, itemtype=itemtype, pars = sv, verbose = FALSE,
                             parprior=parprior, PrepList=PrepList, large=large, calcNull=FALSE,
                             technical=technical, constrain=constrain, ...)
        # coef(tmpmod, simplify=TRUE)
        ret <- list(LL=tmpmod@Fit$logLik + tmpmod@Fit$logPrior, vals=mod2values(tmpmod))
        ret
    }

    f.min <- function(value, dat, model, which, sv, get.LL, large, parprior, asigns,
                      PrepList, itemtype, constrain, print_debug = FALSE, ...){
        sv$est[which] <- FALSE
        sv$value[which] <- value
        drop_contraint <- logical(length(constrain))
        if(length(constrain)){
            for(i in length(constrain):1L){
                if(which %in% constrain[[i]]){
                    sv$est[constrain[[i]]] <- FALSE
                    sv$value[constrain[[i]]] <- value
                    constrain[[i]] <- NULL
                }
            }
        }
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
                                 PrepList=PrepList, constrain=constrain, ...), silent=TRUE)
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

    f.min2 <- function(value, upperBound, ...){
        if(upperBound) value <- -value
        f.min(value=value, ...)^2 + value
    }

    LLpar <- function(X, parnums, dat, model, large, constrain, direction,
                      sv, get.LL, parprior, asigns, PrepList, inf2val, itemtype,
                      maxLL, estlower, estupper, search_bound, step, NealeMiller, ...){
        parnum <- parnums[X]
        TOL <- .001
        lower <- ifelse(sv$lbound[parnum] == -Inf, -inf2val, sv$lbound[parnum])
        upper <- ifelse(sv$ubound[parnum] == Inf, inf2val, sv$ubound[parnum])
        mid <- sv$value[parnum]
        if(closeEnough(lower, -1e-2, 1e-2) && closeEnough(upper, 1 + -1e-2, 1 + 1e-2))
            step <- step/10
        if(closeEnough(lower, -1e-2, 1e-2) && upper > 10L) step <- step/3
        conv <- TRUE
        if(direction[X] == 'lower'){
            possible_bound <- TRUE
            if(search_bound){
                grid <- mid - cumsum(rep(step, floor(abs(upper/step))))
                for(g in grid){
                    Xval <- f.min(g, dat=dat, model=model, constrain=constrain,
                                   large=large, which=parnum, sv=sv, get.LL=get.LL,
                                   parprior=parprior, asigns=asigns,
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
            if(possible_bound){
                if(NealeMiller){
                    opt <- try(optimize(f.min2, c(lower, mid), upperBound=FALSE, dat=dat, model=model,
                                           large=large, which=parnum, sv=sv, get.LL=get.LL,
                                           parprior=parprior, asigns=asigns,
                                           PrepList=PrepList, itemtype=itemtype, constrain=constrain,
                                           ..., f.upper=maxLL-get.LL, tol = TOL),
                                     silent = TRUE)
                    if(!is(opt, 'try-error'))
                        opt <- list(root=opt$minimum, f.root=TOL/10)
                } else {
                    opt <- try(uniroot(f.min, c(lower, mid), dat=dat, model=model,
                                             large=large, which=parnum, sv=sv, get.LL=get.LL,
                                             parprior=parprior, asigns=asigns,
                                             PrepList=PrepList, itemtype=itemtype, constrain=constrain,
                                             ..., f.upper=maxLL-get.LL, tol = TOL/10),
                                     silent = TRUE)
                }
            } else opt <- try(uniroot(), TRUE)
            if(is(opt, 'try-error')) opt <- list(root = lower, f.root=1e10)
        } else if(direction[X] == 'upper'){
            possible_bound <- TRUE
            if(search_bound){
                grid <- mid + cumsum(rep(step, floor(abs(upper/step))))
                for(g in grid){
                    Xval <- f.min(g, dat=dat, model=model, constrain=constrain,
                                   large=large, which=parnum, sv=sv, get.LL=get.LL,
                                   parprior=parprior, asigns=asigns,
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
            if(possible_bound){
                if(NealeMiller){
                    opt <- try(optimize(f.min2, c(mid, upper), upperBound=TRUE, dat=dat, model=model,
                                              large=large, which=parnum, sv=sv, get.LL=get.LL,
                                              parprior=parprior, asigns=asigns,
                                              PrepList=PrepList, itemtype=itemtype, constrain=constrain,
                                              ..., f.upper=maxLL-get.LL, tol = TOL),
                                     silent = TRUE)
                    if(!is(opt, 'try-error'))
                        opt <- list(root=opt$minimum, f.root=TOL/10)
                } else {
                    opt <- try(uniroot(f.min, c(mid, upper), dat=dat, model=model,
                                             large=large, which=parnum, sv=sv, get.LL=get.LL,
                                             parprior=parprior, asigns=asigns,
                                             PrepList=PrepList, itemtype=itemtype, constrain=constrain,
                                             ..., f.lower=maxLL-get.LL, tol = TOL/10),
                                     silent = TRUE)
                }
            } else opt <- try(uniroot(), TRUE)
            if(is(opt, 'try-error')) opt <- list(root = upper, f.root=1e10)
        } else opt <- list(root = upper, f.root=1e10)
        if(abs(opt$f.root) > TOL) conv <- FALSE
        c(CI=opt$root, conv=conv)
    }

    stopifnot(extract.mirt(mod, 'converged'))
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
        parnums <- sv$parnum[tmp]
    } else {
        parnums <- sv$parnum[sv$est]
    }
    LL <- mod@Fit$logLik
    get.LL <- LL - qchisq(1-alpha, 1)/2
    constraints <- extract.mirt(mod, 'constrain')
    if(length(constraints)){
        for(i in 1L:length(constraints)){
            match <- which(parnums %in% constraints[[i]])
            if(length(match) > 1L){
                match <- match[-1L]
                parnums <- parnums[-match]
            }
        }
    }
    direction <- rep(c('lower', 'upper'), length.out = length(parnums)*2)
    if(!lower) direction <- direction[-seq(1L, length(direction), by=2)]
    if(!upper) direction <- direction[-seq(2L, length(direction), by=2)]
    parnumsold <- parnums
    if(lower && upper)
        parnums <- rep(parnums, each = 2)
    X <- 1L:length(parnums)
    result <- mySapply(X=X, FUN=LLpar, parnums=parnums, asigns=asigns,
                       dat=dat, constrain=constraints, itemtype=itemtype,
                       model=model, large=large, sv=sv, get.LL=get.LL, parprior=parprior,
                       PrepList=PrepList, inf2val=inf2val, maxLL=LL,
                       direction=direction, search_bound=search_bound,
                       step=step, NealeMiller=NealeMiller, ...)
    lowerCIs <- lowerconv <- upperCIs <- upperconv <- NA
    if(lower){
        lowerCIs <- result[direction == 'lower', 'CI']
        lowerconv <- as.logical(result[direction == 'lower', 'conv'])
    }
    if(upper){
        upperCIs <- result[direction == 'upper', 'CI']
        upperconv <- as.logical(result[direction == 'upper', 'conv'])
    }
    result <- data.frame(lower=lowerCIs, upper=upperCIs, lower_conv=lowerconv,
                         upper_conv=upperconv)
    colnames(result) <- c(paste0('lower_', alpha/2*100), paste0('upper_', (1-alpha/2)*100),
                                 'lower_conv', 'upper_conv')
    ret <- data.frame(Item=sv$item[parnumsold], class=sv$class[parnumsold], parnam=sv$name[parnumsold],
                      parnum=parnumsold, value=sv$value[parnumsold], result, row.names=NULL)
    if(!lower) ret <- ret[,!grepl('lower', colnames(ret)), drop=FALSE]
    if(!upper) ret <- ret[,!grepl('upper', colnames(ret)), drop=FALSE]
    ret
}
