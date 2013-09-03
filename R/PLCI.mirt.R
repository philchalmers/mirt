#' Compute profiled-likelihood confidence intervals
#' 
#' Computes profiled-likelihood based confidence intervals. Supports the inclusion of prior parameter 
#' distributions as well as equality constraints. For multidimensional models, the CI's for the slopes 
#' are not estimated due to the possibility of signs flipping during estimation. In unidimenisonal models, the 
#' slope parameters are assumed to be greater than zero, and a lower bound is imposed to ensure that sign
#' flipping does not occur.
#' 
#' @aliases PLCI.mirt
#' @param mod a converged mirt model
#' @param alpha two-tailed alpha critical level
#' @param parnum a numeric vector indicating which parameters to estimate. Use \code{\link{mod2values}}
#' to determine parameter numbers. If \code{NULL}, all possible parameters are used
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
PLCI.mirt <- function(mod, alpha = .05, parnum = NULL){
    
    compute.LL <- function(dat, model, sv, large, parprior){        
        tmpmod <- suppressMessages(mirt::mirt(dat, model, pars = sv, verbose = FALSE, parprior=parprior,
                                        large=large, calcNull=FALSE))
        tmpmod@logLik
    }
    
    f.min <- function(value, dat, model, which, sv, get.LL, large, parprior, parnames){        
        sv$est[which] <- FALSE
        sv$value[which] <- value
        if(any(parnames[which] %in% paste0('d', 1L:10L))){            
            intnum <- as.numeric(strsplit(as.character(parnames[which]), 'd')[[1L]][2L])
            #test if graded model            
            if(sv$name[which-intnum] != 'd0'){    
                about <- parnums[which]
                if(intnum != 1L) about <- sort(which:(which-intnum+1L))
                len <- do.call(c, lapply(strsplit(as.character(parnames[(which+1L):(which+10L)]), 'd'), length))
                if(which(len == 1L)[1L] > 1L)
                    about <- c(about, max(about) + 1L:(which(len == 1L)[1L]-1L))
                vals <- rep(value, length(about))
                if(intnum == 1L){                    
                    vals[2L:length(vals)] <- vals[2L:length(vals)] - seq(0.1, 0.5, length.out = length(vals)-1L)
                } else if(intnum == length(about)){
                    vals[1L:(length(vals)-1L)] <- vals[1L:(length(vals)-1L)] + 
                        seq(0.5, 0.1, length.out = length(1L:(length(vals)-1L)))
                } else {
                    vals[1L:(intnum-1L)] <- vals[1L:(intnum-1L)] + seq(0.5, 0.1, length.out = intnum-1L)
                    vals[(intnum+1L):length(vals)] <- vals[(intnum+1L):length(vals)] -
                        seq(0.1, 0.5, length.out = length((intnum+1L):length(vals)))                                        
                }
                sv$value[about] <- vals
            }
        }            
        got.LL <- try(compute.LL(dat=dat, model=model, sv=sv, large=large, parprior=parprior), silent=TRUE)
        if(is(got.LL, 'try-error')) return(1e8)
        ret <- (got.LL - get.LL)^2
        attr(ret, 'value') <- value
        ret
    }
    
    LLpar <- function(parnum, parnums, parnames, lbound, ubound, dat, model, large, sv, get.LL, parprior){
        lower <- ifelse(lbound[parnum] == -Inf, -10, lbound[parnum])
        upper <- ifelse(ubound[parnum] == Inf, 10, ubound[parnum])
        opt.lower <- optimize(f.min, lower = lower, upper = pars[parnum], dat=dat, model=model, large=large,
                              which=parnums[parnum], sv=sv, get.LL=get.LL, parprior=parprior, parnames=parnames,
                              tol = .01)        
        opt.upper <- optimize(f.min, lower = pars[parnum], upper = upper, dat=dat, model=model, large=large,
                              which=parnums[parnum], sv=sv, get.LL=get.LL, parprior=parprior, parnames=parnames,
                              tol = .01)
        c(lower=opt.lower$minimum, upper=opt.upper$minimum)
    }
    
    dat <- mod@data
    model <- mod@model[[1L]]    
    parprior <- mod@parprior 
    if(length(parprior) == 0L) parprior <- NULL
    sv <- mod2values(mod)
    large <- mirt(mod@data, mod@model[[1L]], large = TRUE)
    #set lbounds to 0 to avoid sign flipping in slopes
    sv$lbound[sv$name == 'a1'] <- 0
    if(!is.null(parnum)){
        tmp <- sv$parnum %in% parnum
        pars <- sv$value[tmp]    
        LL.upper.crit <- LL.lower.crit <- pars
        parnums <- sv$parnum[tmp]
        parnames <- sv$name[tmp]
        lbound <- sv$lbound[tmp]
        ubound <- sv$ubound[tmp]
    } else {
        pars <- sv$value    
        pars <- LL.upper.crit <- LL.lower.crit <- pars[sv$est]    
        parnums <- sv$parnum[sv$est]
        parnames <- sv$name[sv$est]
        lbound <- sv$lbound[sv$est]
        ubound <- sv$ubound[sv$est]        
    }    
    if(mod@nfact > 1L && is.null(parnum)){
        #don't estimate slopes of multidimensional
        tmp <- !(sv$name %in% paste0('a', 1:20)) & sv$est
        pars <- sv$value[tmp]
        parnames <- sv$name[tmp]
        parnums <- sv$parnum[tmp]
        lbound <- sv$lbound[tmp]
        ubound <- sv$ubound[tmp]
    }
    LL <- mod@logLik
    get.LL <- LL - qchisq(1-alpha, 1)/2
    if(!is.null(globalenv()$MIRTCLUSTER)){
        result <- t(parallel::parSapply(cl=globalenv()$MIRTCLUSTER, 1L:length(parnums), LLpar, 
                                        parnums=parnums, parnames=parnames, lbound=lbound, ubound=ubound, 
                                        dat=dat, model=model,
                                        large=large, sv=sv, get.LL=get.LL, parprior=parprior))
    } else {
        result <- t(sapply(1L:length(parnums), LLpar, parnums=parnums, parnames=parnames, lbound=lbound, 
                           ubound=ubound, dat=dat, model=model, large=large, sv=sv, get.LL=get.LL, parprior=parprior))           
    }
    colnames(result) <- c(paste0('lower_', alpha/2*100), paste0('upper_', (1-alpha/2)*100)) 
    ret <- data.frame(Item=sv$item[parnums], parnam=sv$name[parnums], value=pars, result, row.names=NULL)
    ret
}
