#' Compute profiled-likelihood confidence intervals
#' 
#' Computes profiled-likelihood based confidence intervals. Supports the inclusion of prior parameter 
#' distributions as well as equality constraints. Currently only supports unidimensional 
#' dichotomous response models. 
#' 
#' @aliases PLCI.mirt
#' @param mod a converged mirt model
#' @param alpha two-tailed alpha critical level
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
#' }
PLCI.mirt <- function(mod, alpha = .05){
    
    compute.LL <- function(dat, model, sv, large, parprior){        
        tmpmod <- suppressMessages(mirt::mirt(dat, model, pars = sv, verbose = FALSE, parprior=parprior,
                                        large=large, calcNull=FALSE))
        tmpmod@logLik
    }
    
    f.min <- function(value, dat, model, which, sv, get.LL, large, parprior){
        sv$est[which] <- FALSE
        sv$value[which] <- value
        got.LL <- compute.LL(dat=dat, model=model, sv=sv, large=large, parprior=parprior)
        ret <- (got.LL - get.LL)^2
        attr(ret, 'value') <- value
        ret
    }
    
    LLpar <- function(parnum, parnums, lbound, ubound, dat, model, large, sv, get.LL, parprior){        
        lower <- ifelse(lbound[parnum] == -Inf, -5, lbound[parnum])
        upper <- ifelse(ubound[parnum] == Inf, 10, ubound[parnum])
        opt.lower <- optimize(f.min, lower = lower, upper = pars[parnum], dat=dat, model=model, large=large,
                              which=parnums[parnum], sv=sv, get.LL=get.LL, parprior=parprior, 
                              tol = .01)        
        opt.upper <- optimize(f.min, lower = pars[parnum], upper = upper, dat=dat, model=model, large=large,
                              which=parnums[parnum], sv=sv, get.LL=get.LL, parprior=parprior, 
                              tol = .01)
        c(lower=opt.lower$minimum, upper=opt.upper$minimum)
    }
    
    if(!all(sapply(mod@pars, class) %in% c('dich', 'GroupPars')))
        stop('Likelihood confidence intervals only available for unidimensional dichotomous models.')
    dat <- mod@data
    model <- mod@model[[1L]]    
    parprior <- mod@parprior 
    if(length(parprior) == 0L) parprior <- NULL
    sv <- mod2values(mod)
    large <- mirt(mod@data, mod@model[[1L]], large = TRUE)
    #set lbounds to 0 to avoid sign flipping in slopes
    sv$lbound[sv$name == 'a1'] <- 0
    pars <- sv$value    
    pars <- LL.upper.crit <- LL.lower.crit <- pars[sv$est]    
    parnums <- sv$parnum[sv$est]
    lbound <- sv$lbound[sv$est]
    ubound <- sv$ubound[sv$est]
    LL <- mod@logLik
    get.LL <- LL - qchisq(1-alpha, 1)/2
    if(!is.null(globalenv()$MIRTCLUSTER)){
        result <- t(parallel::parSapply(cl=globalenv()$MIRTCLUSTER, 1L:length(parnums), LLpar, 
                                        parnums=parnums, lbound=lbound, ubound=ubound, dat=dat, model=model,
                                        large=large, sv=sv, get.LL=get.LL, parprior=parprior))
    } else {
        result <- t(sapply(1L:length(parnums), LLpar, parnums=parnums, lbound=lbound, ubound=ubound, 
                           dat=dat, model=model, large=large, sv=sv, get.LL=get.LL, parprior=parprior))           
    }
    colnames(result) <- c(paste0('lower_', alpha/2*100), paste0('upper_', (1-alpha/2)*100)) 
    ret <- data.frame(Item=sv$item[sv$est], parnam=sv$name[sv$est], value=pars, result, row.names=NULL)
    ret
}
