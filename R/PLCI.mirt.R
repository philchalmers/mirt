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
#' @param plot logical; plot the parameter relationship in the likelihood space for two parameters?
#' @param npts number of points to evaluate and plot if \code{plot = TRUE}
#' @param ... additional arguments to pass to the estimation functions
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
#' # plot the confidence envelope for parameters 1 and 2
#' PLCI.mirt(mod2, parnum=c(1,2), plot=TRUE)
#'
#' }
PLCI.mirt <- function(mod, alpha = .05, parnum = NULL, plot = FALSE, npts = 24, ...){

    #silently accepts print_debug = TRUE for printing the minimization criteria

    compute.LL <- function(dat, model, sv, large, parprior, ...){
        tmpmod <- mirt::mirt(dat, model, pars = sv, verbose = FALSE, parprior=parprior,
                                        large=large, calcNull=FALSE, technical=list(message=FALSE, warn=FALSE,
                                                                                    parallel=FALSE), ...)
        ret <- list(LL=tmpmod@Fit$logLik + tmpmod@Fit$logPrior, vals=mod2values(tmpmod))
        ret
    }

    f.min <- function(value, dat, model, which, sv, get.LL, large, parprior, parnames, asigns,
                      print_debug = FALSE, ...){
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
        got.LL <- try(compute.LL(dat=dat, model=model, sv=sv, large=large, parprior=parprior, ...),
                      silent=TRUE)
        if(is(got.LL, 'try-error')) return(1e10)
        sv2 <- got.LL$vals
        got.LL <- got.LL$LL
        as <- matrix(sv2$value[sv2$name %in% paste0('a', 1L:30L)], ncol(dat))
        if(sum(asigns * sign(as)) < 0L) return(1e10)
        ret <- (got.LL - get.LL)^2
        ret <- ifelse(is.finite(ret), ret, 1e10)
        if(print_debug) cat('parnum = ', which, '; value = ', round(value, 3),
                            '; min = ', round(ret, 3), '\n')
        attr(ret, 'value') <- value
        ret
    }

    LLpar <- function(parnum, parnums, parnames, lbound, ubound, dat, model, large,
                      sv, get.LL, parprior, asigns, single=FALSE, force = FALSE, ...){
        lower <- ifelse(lbound[parnum] == -Inf, -15, lbound[parnum])
        upper <- ifelse(ubound[parnum] == Inf, 15, ubound[parnum])
        mid <- pars[parnum]
        if(parnames[parnum] %in% c('g', 'u')){
            lower <- 0
            upper <- 1
        } else if(parnames[parnum] %in% paste0('COV_', 1:30, 1:30)){
            lower <- 0
        }
        if(single){
            return(optimize(f.min, lower = lower, upper = upper, dat=dat, model=model,
                            large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                            parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)$minimum)
        }
        if(mid > lower){
            opt.lower <- optimize(f.min, lower = lower, upper = mid, dat=dat, model=model,
                                  large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                  parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)
            if(opt.lower$objective > .01){
                tmp <- optim(mid - abs((mid - lower) * .01), f.min, dat=dat, model=model,
                             large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                             parprior=parprior, parnames=parnames, asigns=asigns, ...,
                             method = 'L-BFGS-B', lower = lbound[parnum], upper = mid,
                             control = list(factr=1e10))
                opt.lower$minimum <- tmp$par; opt.lower$objective <- tmp$value
            }
        } else opt.lower <- list(minimum = lower, objective=0)
        if(mid < upper){
            opt.upper <- optimize(f.min, lower = mid, upper = upper, dat=dat, model=model,
                                  large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                  parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)
            if(opt.upper$objective > .01){
                tmp <- optim(mid + abs((upper - mid) * .01), f.min, dat=dat, model=model,
                             large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                             parprior=parprior, parnames=parnames, asigns=asigns, ...,
                             method = 'L-BFGS-B', lower = mid, upper = ubound[parnum],
                             control = list(factr=1e10))
                opt.upper$minimum <- tmp$par; opt.upper$objective <- tmp$value
            }
        } else opt.upper <- list(minimum = upper, objective=0)
        if(force){ #TODO this is pretty hacky, but it works for the most part
            if(opt.upper$objective > .01){
                opt.upper <- optimize(f.min, lower = (opt.lower$minimum + mid)/2, upper = mid, dat=dat, model=model,
                                      large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                      parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)
                if(opt.upper$objective > .01){
                    opt.upper <- optimize(f.min, upper = opt.lower$minimum, lower = lbound[2], dat=dat, model=model,
                                          large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                          parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)
                }
            } else if(opt.lower$objective > .01){
                opt.lower <- optimize(f.min, lower = mid, upper = (opt.upper$minimum + mid)/2, dat=dat, model=model,
                                      large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                      parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)
                if(opt.lower$objective > .01){
                    opt.lower <- optimize(f.min, lower = opt.upper$minimum, upper = ubound[2L], dat=dat, model=model,
                                          large=large, which=parnums[parnum], sv=sv, get.LL=get.LL,
                                          parprior=parprior, parnames=parnames, asigns=asigns, ..., tol = .01)
                }
            }
        }
        conv_upper <- conv_lower <- TRUE
        if(opt.lower$objective > .01){
            conv_lower <- FALSE
            if(force) opt.lower$minimum <- NA
        }
        if(opt.upper$objective > .01){
            conv_upper <- FALSE
            if(force) opt.upper$minimum <- NA
        }
        c(lower=opt.lower$minimum, upper=opt.upper$minimum,
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
    if(plot){
        if(length(parnum) != 2L)
            stop('parnum input must contain exactly two parameter numbers', call.=FALSE)
    }
    LL <- mod@Fit$logLik
    get.LL <- LL - qchisq(1-alpha, 1 + plot)/2
    result <- mySapply(X=1L:length(parnums), FUN=LLpar, parnums=parnums, asigns=asigns,
                       parnames=parnames, lbound=lbound, ubound=ubound, dat=dat,
                       model=model, large=large, sv=sv, get.LL=get.LL, parprior=parprior, ...)
    colnames(result) <- c(paste0('lower_', alpha/2*100), paste0('upper_', (1-alpha/2)*100),
                          'lower_conv', 'upper_conv')
    ret <- data.frame(Item=sv$item[parnums], class=itemtypes, parnam=sv$name[parnums],
                      parnum=parnums, value=pars, result, row.names=NULL)
    ret$lower_conv <- as.logical(ret$lower_conv)
    ret$upper_conv <- as.logical(ret$upper_conv)
    if(plot){
        ret <- rbind(ret[ret$parnum == parnum[1],], ret[ret$parnum != parnum[1],])
        parnames <- ret$parnam; parnums <- ret$parnum
        xrange <- seq(from=ret[1L, 6L], to=ret[1L, 7L], length.out = floor((npts-2)/2))
        xrange[1L] <- mean(xrange[1:2])
        xrange[length(xrange)] <- mean(xrange[length(xrange):(length(xrange)-1)])
        lbound[2L] <- ret[2L, 6]
        ubound[2L] <- ret[2L, 7]
        sv2 <- sv
        sv2$est[sv2$parnum == parnums[1L]] <- FALSE
        collect <- matrix(NA, length(xrange), 2L)
        for(i in 1L:length(xrange)){
            sv2$value[sv2$parnum == parnums[1L]] <- xrange[i]
            result <- mySapply(X=2L, FUN=LLpar, parnums=parnums, asigns=asigns,
                               parnames=parnames, lbound=lbound, ubound=ubound, dat=dat,
                               model=model, large=large, sv=sv2, get.LL=get.LL, parprior=parprior,
                               force = TRUE, ...)
            collect[i, ] <- result[1:2]
        }
        sv2$value[sv2$parnum == parnums[1L]] <- ret[1L, 6L]
        lp <- mySapply(X=2L, FUN=LLpar, parnums=parnums, asigns=asigns,
                           parnames=parnames, lbound=lbound, ubound=ubound, dat=dat,
                           model=model, large=large, sv=sv2, get.LL=get.LL, parprior=parprior,
                           force = TRUE, single=TRUE, ...)
        sv2$value[sv2$parnum == parnums[1L]] <- ret[1L, 7L]
        up <- mySapply(X=2L, FUN=LLpar, parnums=parnums, asigns=asigns,
                       parnames=parnames, lbound=lbound, ubound=ubound, dat=dat,
                       model=model, large=large, sv=sv2, get.LL=get.LL, parprior=parprior,
                       force = TRUE, single=TRUE, ...)
        dat <- data.frame(x=xrange, y=as.numeric(collect))
        dat <- rbind(dat, c(ret[1L, 6L], lp), c(ret[1L, 7L], up))
        dat <- rbind(dat, ret[,'value'])
        dat$group <- factor(c(rep('pts', nrow(dat)-1), 'est'))
        return(xyplot(y ~ x, dat, type = 'p', groups=dat$group, col=c('black', 'blue'),
               main = 'Likelihood Confidence Envelope',
               xlab = paste0(ret[1,'parnam'], ' (#', ret[1,'parnum'], ')'),
               ylab = paste0(ret[2,'parnam'], ' (#', ret[2,'parnum'], ')')))
    }
    ret
}
