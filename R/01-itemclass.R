###################################################

# valid itemtype inputs

# flag to indicate an experimental item type (requires an S4 initializer in the definitions below)
Experimental_itemtypes <- function() c('experimental', 'grsmIRT', 'crm')

Valid_iteminputs <- function() c('Rasch', '2PL', '3PL', '3PLu', '4PL', 'graded', 'grsm', 'gpcm', 'gpcmIRT',
                                 'rsm', 'nominal', 'PC2PL','PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM',
                                 'ideal', 'lca', 'spline', 'monopoly', 'ggum', 'sequential', 'Tutz', Experimental_itemtypes())

ordinal_itemtypes <- function() c('dich', 'graded', 'gpcm', 'sequential', 'ggum', 'rating', 'spline', 'monopoly',
                                  'partcomp', 'rsm', 'ideal', 'gpcmIRT', 'grsmIRT')

Continuous_itemtypes <- function() c('crm')

# Indicate which functions should use the R function instead of those written in C++
Use_R_ProbTrace <- function() c('custom', 'spline', 'sequential', 'Tutz', Experimental_itemtypes())

Use_R_Deriv <- function() c('custom', 'rating', 'partcomp', 'nestlogit', 'spline', 'sequential', 'Tutz',
                            Experimental_itemtypes())

###################################################
#Generic Item class

setClass("AllItemsClass",
         representation(par='numeric',
                        SEpar='numeric',
                        parnames='character',
                        est='logical',
                        dps='function',
                        dps2='function',
                        constr='logical',
                        itemclass='integer',
                        parnum='numeric',
                        nfact='integer',
                        nfixedeffects='numeric',
                        fixed.design='matrix',
                        dat='matrix',
                        orgdat='matrix',
                        ncat='integer',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric',
                        'VIRTUAL')
)

#--------------------------------------------------------------------------

#Generics

setGeneric('ProbTrace', function(x, Theta, ...) standardGeneric("ProbTrace"))

setGeneric('ExtractLambdas', function(x) standardGeneric("ExtractLambdas"))

setGeneric('ExtractZetas', function(x) standardGeneric("ExtractZetas"))

setGeneric('Deriv', function(x, Theta, ...) standardGeneric("Deriv"))

setGeneric('DerivTheta', function(x, Theta) standardGeneric("DerivTheta"))

setGeneric('dP', function(x, Theta) standardGeneric("dP"))

setGeneric('calcLogLik', function(object, ...) standardGeneric("calcLogLik"))

setGeneric('set_null_model', function(x) standardGeneric('set_null_model'))

setGeneric("itemplot.internal",  function(object, ...) standardGeneric("itemplot.internal"))

setGeneric("fscores.internal", function(object, ...) standardGeneric("fscores.internal"))

setGeneric("summary", function(object, ...) standardGeneric("summary"))

setGeneric("vcov", function(object, ...) standardGeneric("vcov"))

setGeneric("logLik", function(object, ...) standardGeneric("logLik"))

setGeneric('DrawValues', function(x, Theta, ...) standardGeneric("DrawValues"))

setGeneric('GenRandomPars', function(x) standardGeneric("GenRandomPars"))

setGeneric('CheckIntercepts', function(x) standardGeneric("CheckIntercepts"))


# ----------------------------------------------------------------
# helper functions

EML <- function(par, obj, Theta){
    obj@par[obj@est] <- par
    itemtrace <- ProbTrace(x=obj, Theta=Theta)
    LL <- sum(obj@dat * log(itemtrace))
    LL <- LL.Priors(x=obj, LL=LL)
    return(LL)
}

EML2 <- function(x, Theta, pars, tabdata, freq, itemloc, CUSTOM.IND, bfactor_info){
    obj <- pars[[length(pars)]]
    obj@par[obj@est] <- x
    gp <- ExtractGroupPars(obj)
    mu <- gp$gmeans
    sigma <- gp$gcov
    prior <- mirt_dmvnorm(Theta, mean=mu, sigma=sigma)
    prior <- prior/sum(prior)
    if(obj@dentype == 'bfactor'){
        J <- length(itemloc) - 1L
        sitems <- bfactor_info$sitems
        nfact <- bfactor_info$nfact
        theta <- pars[[J+1L]]@theta
        Thetabetween <- pars[[J+1L]]@Thetabetween
        p <- matrix(0, nrow(Theta), ncol(sitems))
        pp <- matrix(0, nrow(theta), ncol(sitems))
        for(i in seq_len(ncol(sitems))){
            sel <- c(seq_len(nfact-ncol(sitems)), i + nfact - ncol(sitems))
            p[,i] <- mirt_dmvnorm(Theta[ ,sel], gp$gmeans[sel], gp$gcov[sel,sel,drop=FALSE])
            pp[,i] <- dnorm(theta, gp$gmeans[sel[length(sel)]],
                            sqrt(gp$gcov[sel[length(sel)],sel[length(sel)],drop=FALSE]))
        }
        pb <- mirt_dmvnorm(Thetabetween, gp$gmeans[seq_len(ncol(Thetabetween))],
                           gp$gcov[seq_len(ncol(Thetabetween)), seq_len(ncol(Thetabetween)), drop=FALSE])
        Priorbetween <- pb / sum(pb)
        prior <- t(t(pp) / colSums(pp))
        rlist <- Estep.bfactor(pars=pars, tabdata=tabdata, freq=freq,
                               Theta=Theta, prior=prior,
                               Priorbetween=Priorbetween, specific=bfactor_info$specific,
                               sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND, omp_threads=1L)
    } else {
        rlist <- Estep.mirt(pars=pars, tabdata=tabdata, freq=freq,
                            Theta=Theta, prior=prior, itemloc=itemloc,
                            CUSTOM.IND=CUSTOM.IND, full=FALSE, omp_threads=1L)
    }
    tmp <- log(rlist$expected)
    pick <- is.finite(tmp)
    LL <- sum(freq[pick]*tmp[pick])
    LL <- LL.Priors(x=obj, LL=LL)
    return(LL)
}

difexp <- function(x) x * (1 - x)

dif2exp <- function(x) 2 * (x * (1 - x)^2)

numDeriv_DerivTheta <- function(item, Theta){
    P <- function(Theta, item, cat) probtrace(item, Theta)[cat]
    grad <- hess <- vector('list', item@ncat)
    tmp <- tmp2 <- matrix(0, nrow(Theta), ncol(Theta))
    for(j in seq_len(item@ncat)){
        for(i in seq_len(nrow(Theta))){
            tmp[i, ] <- numerical_deriv(Theta[i, , drop=FALSE], P, item=item, cat=j)
            tmp2[i, ] <- diag(numerical_deriv(Theta[i, , drop=FALSE], P, item=item, cat=j,
                                              gradient=FALSE))
        }
        grad[[j]] <- tmp
        hess[[j]] <- tmp2
    }
    return(list(grad=grad, hess=hess))
}

numDeriv_dP <- function(item, Theta){
    P <- function(par, Theta, item, cat){
        item@par[item@est] <- par
        sum(ProbTrace(item, Theta)[cat:item@ncat])
    }
    par <- item@par[item@est]
    ret <- matrix(0, nrow(Theta), length(item@par))
    for(i in seq_len(nrow(Theta))){
        tmp <- numeric(length(par))
        for(j in seq_len(item@ncat))
            tmp <- tmp + numerical_deriv(par, P, Theta=Theta[i, , drop=FALSE],
                                         item=item, cat=j)
        ret[i, item@est] <- tmp
    }
    ret
}

numDeriv_dP2 <- function(item, Theta){
    P <- function(par, Theta, item, cat){
        item@par[item@est] <- par
        ProbTrace(item, Theta)[cat]
    }
    par <- item@par[item@est]
    tmpmat <- matrix(0, nrow(Theta), length(item@par))
    ret <- lapply(2L:item@ncat - 1L, function(x) tmpmat)
    for(i in seq_len(nrow(Theta))){
        for(j in 2L:item@ncat)
            ret[[j-1L]][i, item@est] <-
                numerical_deriv(par, P, Theta=Theta[i, , drop=FALSE],
                                item=item, cat=j)
    }
    ret
}

symbolicGrad_par <- function(x, Theta, dp1 = NULL, P = NULL){
    if(is.null(P)) P <- ProbTrace(x, Theta)
    xLength <- length(x@par)
    r_P <- x@dat / P
    if(is.null(dp1))
        dp1 <- array(x@dps(x@par, Theta, x@ncat), c(nrow(Theta),x@ncat,xLength))
    grad <- numeric(length(x@par))
    for (i in 1L:xLength)
        grad[i] <- sum(r_P * dp1[,,i])
    grad
}

symbolicHessian_par <- function(x, Theta, dp1 = NULL, dp2 = NULL, P = NULL){
    if(is.null(P)) P <- ProbTrace(x, Theta)
    xLength <- length(x@par)
    ThetaLength <- length(Theta)
    if(is.null(dp1))
        dp1 <- array(x@dps(x@par, Theta, x@ncat), c(ThetaLength,x@ncat,xLength))
    if(is.null(dp2))
        dp2 <- array(x@dps2(x@par, Theta, x@ncat), c(ThetaLength,x@ncat,xLength,xLength))
    H <- matrix(0,xLength,xLength)
    P2 <- P^2
    for (i in 1L:xLength){
        for (j in i:xLength){
            H[i,j] <- sum(x@dat*dp2[,,i,j]/P + x@dat*dp1[,,i]*(-dp1[,,j]/P2))
            H[j,i] <- H[i,j]
        }
    }
    H
}
