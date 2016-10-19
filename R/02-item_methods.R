# ----------------------------------------------------------------
# helper functions

EML <- function(par, obj, Theta){
    obj@par[obj@est] <- par
    itemtrace <- ProbTrace(x=obj, Theta=Theta)
    LL <- sum(obj@dat * log(itemtrace))
    LL <- LL.Priors(x=obj, LL=LL)
    return(LL)
}

EML2 <- function(x, Theta, pars, tabdata, freq, itemloc, CUSTOM.IND){
    obj <- pars[[length(pars)]]
    obj@par[obj@est] <- x
    gpars <- ExtractGroupPars(obj)
    mu <- gpars$gmeans
    sigma <- gpars$gcov
    prior <- mirt_dmvnorm(Theta, mean=mu, sigma=sigma)
    prior <- prior/sum(prior)
    rlist <- Estep.mirt(pars=pars, tabdata=tabdata, freq=freq,
                        Theta=Theta, prior=prior, itemloc=itemloc,
                        CUSTOM.IND=CUSTOM.IND, full=FALSE)
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
    for(j in 1L:item@ncat){
        for(i in 1L:nrow(Theta)){
            tmp[i, ] <- numDeriv::grad(P, x=Theta[i, , drop=FALSE], item=item, cat=j)
            tmp2[i, ] <- diag(numDeriv::hessian(P, x=Theta[i, , drop=FALSE], item=item, cat=j))
        }
        grad[[j]] <- tmp
        hess[[j]] <- tmp2
    }
    return(list(grad=grad, hess=hess))
}

numDeriv_dP <- function(item, Theta){
    P <- function(par, Theta, item, cat){
        item@par <- par
        sum(ProbTrace(item, Theta)[cat:item@ncat])
    }
    par <- item@par
    ret <- matrix(0, nrow(Theta), length(item@par))
    for(i in 1L:nrow(Theta)){
        tmp <- numeric(length(par))
        for(j in 1L:item@ncat)
            tmp <- tmp + numDeriv::grad(P, x=par, Theta=Theta[i, , drop=FALSE],
                              item=item, cat=j)
        ret[i, ] <- tmp
    }
    ret
}

# ----------------------------------------------------------------
# valid itemtype inputs

# flag to indicate an experimental item type (requires an S4 initializer in the definitions below)
Experimental_itemtypes <- function() c('experimental', 'grsmIRT')

Valid_iteminputs <- function() c('Rasch', '2PL', '3PL', '3PLu', '4PL', 'graded', 'grsm', 'gpcm',
                                'nominal', 'PC2PL','PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM',
                                'ideal', 'lca', 'spline', Experimental_itemtypes())

# Indicate which functions should use the R function instead of those written in C++
Use_R_ProbTrace <- function() c('custom', 'ideal', 'spline', Experimental_itemtypes())

Use_R_Deriv <- function() c('custom', 'rating', 'rsm', 'partcomp', 'nestlogit',
                            'ideal', 'spline', Experimental_itemtypes())

# ----------------------------------------------------------------
# Begin class and method definitions

setClass("GroupPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        parnum='numeric',
                        itemclass='integer',
                        dat='matrix',
                        nfact='integer',
                        gradient='numeric',
                        hessian='matrix',
                        itemtrace='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric',
                        rr='numeric',
                        rrb='numeric',
                        rrs='matrix',
                        bindex='integer',
                        sindex='matrix',
                        BFACTOR='logical',
                        theta='matrix',
                        Thetabetween='matrix',
                        density='numeric',
                        sig='matrix',
                        invsig='matrix',
                        mu='numeric',
                        meanTheta='numeric',
                        gen='function',
                        den='function',
                        safe_den='function',
                        derivType='character',
                        gr='function',
                        usegr='logical',
                        hss='function',
                        usehss='logical')
)

setMethod(
    f = "print",
    signature = signature(x = 'GroupPars'),
    definition = function(x, ...){
        cat('Object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'GroupPars'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'GroupPars'),
    definition = function(x){
        x
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta, CUSTOM.IND, EM = FALSE, pars = NULL, itemloc = NULL,
                          tabdata = NULL, freq = NULL, estHess=FALSE, prior = NULL){
        if(x@itemclass < 0L){
            LLfun <- function(par, obj, Theta){
                obj@par[obj@est] <- par
                den <- obj@safe_den(obj, Theta)
                LL <- sum(obj@rr * log(den))
                LL
            }
            grad <- rep(0, length(x@par))
            hess <- matrix(0, length(x@par), length(x@par))
            if(any(x@est)){
                if(x@usegr) grad <- x@gr(x, Theta)
                else grad[x@est] <- numerical_deriv(x@par[x@est], LLfun, obj=x, Theta=Theta,
                                                    type=x@derivType)
                if(estHess){
                    if(x@usehss) hess <- x@hss(x, Theta)
                    else hess[x@est, x@est] <-
                            numerical_deriv(x@par[x@est], LLfun, obj=x, Theta=Theta,
                                            gradient=FALSE, type=x@derivType)
                }
            }
            return(list(grad = grad, hess=hess))
        } else {
            if(EM){
                grad <- rep(0, length(x@par))
                hess <- matrix(0, length(x@par), length(x@par))
                if(estHess){
                    if(any(x@est)){
                        hess[x@est,x@est] <- numDeriv::hessian(EML2, x@par[x@est], Theta=Theta,
                                                               pars=pars, tabdata=tabdata, freq=freq,
                                                               itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
                    }
                }
                return(list(grad=grad, hess=hess))
            }
            return(.Call("dgroup", x, Theta, matrix(0), estHess, FALSE, FALSE, FALSE))
        }
    }
)

# ----------------------------------------------------------------

setClass("RandomPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        between='logical',
                        parnum='numeric',
                        ndim='integer',
                        gframe='data.frame',
                        gdesign='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        cand.t.var='numeric',
                        drawvals='matrix',
                        mtch='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric')
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'RandomPars'),
    definition = function(x){
        x
    }
)

setMethod(
    f = "DrawValues",
    signature = signature(x = 'RandomPars', Theta = 'matrix'),
    definition = function(x, Theta, pars, fulldata, itemloc, offterm0, CUSTOM.IND, LR = FALSE){
        J <- length(pars) - 1L
        theta0 <- x@drawvals
        total_0 <- attr(theta0, 'log.lik_full')
        N <- nrow(theta0)
        unif <- runif(N)
        prior.mu <- rep(0, ncol(theta0))
        prior.t.var <- matrix(0, ncol(theta0), ncol(theta0))
        prior.t.var[lower.tri(prior.t.var, diag=TRUE)] <- x@par
        d <- if(ncol(theta0) == 1) matrix(prior.t.var) else diag(diag(prior.t.var))
        prior.t.var <- prior.t.var + t(prior.t.var) - d
        sigma <- if(ncol(theta0) == 1L) matrix(x@cand.t.var) else diag(rep(x@cand.t.var,ncol(theta0)))
        theta1 <- theta0 + mirt_rmvnorm(N, prior.mu, sigma)
        if(is.null(total_0)) theta1 <- theta0 #for intial draw
        log_den1 <- mirt_dmvnorm(theta1,prior.mu,prior.t.var,log=TRUE)
        itemtrace1 <- matrix(0, ncol=ncol(fulldata), nrow=nrow(fulldata))
        if(LR){
            Theta2 <- Theta - rowSums(x@gdesign * theta0[x@mtch, , drop=FALSE]) +
                rowSums(x@gdesign * theta1[x@mtch, , drop=FALSE])
            itemtrace1 <- computeItemtrace(pars, Theta=Theta2, offterm=offterm0,
                                           itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
        } else {
            tmp1 <- rowSums(x@gdesign * theta1[x@mtch, , drop=FALSE])
            offterm1 <- matrix(tmp1, nrow(offterm0), ncol(offterm0))
            itemtrace1 <- computeItemtrace(pars, Theta=Theta, offterm=offterm1,
                                           itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
        }
        LL <- fulldata * log(itemtrace1)
        LL2 <- matrix(0, nrow(LL), J)
        if(LR){
            LL2 <- rowSums(LL)
        } else {
            for(i in 1L:J)
                LL2[,i] <- rowSums(LL[,itemloc[i]:(itemloc[i+1L] - 1L)])
        }
        total_1 <- tapply(LL2, x@mtch, sum) + log_den1
        if(is.null(total_0)){ #for intial draw
            attr(theta1, 'log.lik_full') <- total_1
            return(theta1)
        }
        diff <- total_1 - total_0
        accept <- unif < exp(diff)
        theta1[!accept, ] <- theta0[!accept, ]
        total_1[!accept] <- total_0[!accept]
        attr(theta1, "Proportion Accepted") <- sum(accept)/N
        attr(theta1, 'log.lik_full') <- total_1
        return(theta1)
    }
)

setMethod(
    f = "RandomDeriv",
    signature = signature(x = 'RandomPars'),
    definition = function(x){
        Theta <- x@drawvals
        estHess <- TRUE
        pick <- -c(1L:ncol(Theta))
        out <- .Call("dgroup", x, Theta, matrix(0L), estHess, TRUE, FALSE, FALSE)
        out$grad <- out$grad[pick]
        out$hess <- out$hess[pick, pick, drop=FALSE]
        diag(out$hess) <- -abs(diag(out$hess)) #hack for very small clusters
        out
    }
)

# ----------------------------------------------------------------

setClass("lrPars",
         representation(par='numeric',
                        SEpar='numeric',
                        est='logical',
                        beta='matrix',
                        sigma = 'matrix',
                        mus='matrix',
                        df='data.frame',
                        parnum='numeric',
                        nfact='integer',
                        nfixed='integer',
                        X='matrix',
                        tXX='matrix',
                        inv_tXX='matrix',
                        lbound='numeric',
                        ubound='numeric',
                        any.prior='logical',
                        prior.type='integer',
                        prior_1='numeric',
                        prior_2='numeric',
                        formula='list',
                        EM='logical')
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'lrPars'),
    definition = function(x){
        x
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'lrPars'),
    definition = function(x, cov, theta){
        inv_sigma <- solve(cov)
        tmp <- t(inv_sigma %*% t(theta - x@mus) %*% x@X)
        tmp2 <- -det(inv_sigma) * x@tXX
        list(grad=tmp, hess=tmp2)
    }
)

# ----------------------------------------------------------------

setClass("dich", contains = 'AllItemsClass')

setMethod(
    f = "print",
    signature = signature(x = 'dich'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'dich'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'dich'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'dich'),
    definition = function(x){
        par <- x@par
        d <- par[1L:x@nfact]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'dich'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 rnorm(1L),
                 rnorm(1L, -2, .5),
                 rnorm(1L, 2, .5))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'dich'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.mirt <- function(par, Theta, asMatrix = FALSE, ot = 0)
{
    return(.Call("traceLinePts", par, Theta, ot))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        P <- P.mirt(x@par, Theta=Theta, ot=ot)
        return(P)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        ret <- .Call('dparsDich', x@par, Theta, estHess, x@dat, offterm)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta){
        N <- nrow(Theta)
        nfact <- ncol(Theta)
        parlength <- length(x@par)
        u <- antilogit(x@par[parlength])
        g <- antilogit(x@par[parlength - 1L])
        d <- x@par[parlength - 2L]
        a <- x@par[1L:nfact]
        Pstar <- P.mirt(c(a, d, -999, 999), Theta)[,2L]
        grad <- hess <- vector('list', 2L)
        grad[[1L]] <- grad[[2L]] <- hess[[1L]] <- hess[[2L]] <- matrix(0, N, nfact)
        for(i in 1L:nfact){
            grad[[2L]][ ,i] <- (u-g) * a[i] * (Pstar * (1 - Pstar))
            grad[[1L]][ ,i] <- -1 * grad[[2L]][ ,i]
            hess[[2L]][ ,i] <- 2 * (u - g) * a[i]^2 * ((1 - Pstar)^2 * Pstar) -
                (u - g) * a[i]^2 * (Pstar * (1 - Pstar))
            hess[[1L]][ ,i] <- -1 * hess[[2L]][ ,i]
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta){
        P <- ProbTrace(x, Theta)
        PQ <- apply(P, 1L, prod)
        ret <- cbind(Theta * PQ, PQ, P)
        ret
    }
)

# ----------------------------------------------------------------

setClass("graded", contains = 'AllItemsClass')

setMethod(
    f = "print",
    signature = signature(x = 'graded'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'graded'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'graded'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'graded'),
    definition = function(x){
        par <- x@par
        d <- par[-(1L:x@nfact)]
        d
    }
)

setMethod(
    f = "CheckIntercepts",
    signature = signature(x = 'graded'),
    definition = function(x){
        z <- ExtractZetas(x)
        return(all(sort(z, decreasing = TRUE) == z))
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'graded'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'graded'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.poly <- function(par, Theta, itemexp = FALSE, ot = 0)
{
    return(.Call('gradedTraceLinePts', par, Theta, itemexp, ot, FALSE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.poly(x@par, Theta=Theta, itemexp=itemexp, ot=ot))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        ret <- .Call("dparsPoly", x@par, Theta, offterm, x@dat,
            length(x@par) - ncol(Theta), estHess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- ExtractLambdas(x)
        P <- ProbTrace(x, Theta, itemexp = FALSE)
        grad <- hess <- vector('list', x@ncat)
        for(i in 1L:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1L:x@nfact){
            for(i in 1L:(ncol(P)-1L)){
                w1 <- P[,i] * (1-P[,i]) * a[j]
                w2 <- P[,i+1L] * (1-P[,i+1L]) * a[j]
                grad[[i]][ ,j] <- w1 - w2
                hess[[i]][ ,j] <- a[j]^2 * (2 * P[ ,i] * (1 - P[,i])^2 -
                                                P[ ,i] * (1 - P[,i]) -
                                                2 * P[ ,i+1L] * (1 - P[,i+1L])^2 +
                                                P[ ,i+1L] * (1 - P[,i+1L]))
            }
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){
        P <- ProbTrace(x, Theta, itemexp=FALSE)
        P <- P[,-c(1, ncol(P)), drop=FALSE]
        PQ <- P * (1 - P)
        ret <- cbind(Theta * rowSums(PQ), PQ)
        ret
    }
)

# ----------------------------------------------------------------

setClass("rating", contains = 'AllItemsClass')

setMethod(
    f = "print",
    signature = signature(x = 'rating'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'rating'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'rating'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'rating'),
    definition = function(x){
        par <- x@par
        d <- par[-c(1L:x@nfact, length(par))]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'rating'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE),
                 rnorm(1L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'rating'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.rating <- function(par, Theta, itemexp = FALSE, ot = 0)
{
    return(.Call('gradedTraceLinePts', par, Theta, itemexp, ot, TRUE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.rating(x@par, Theta=Theta, itemexp=itemexp, ot=ot))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        hess <- matrix(0, length(x@par), length(x@par))
        dat <- x@dat
        nfact <- x@nfact
        a <- x@par[1L:nfact]
        d <- ExtractZetas(x)
        nzetas <- length(d)
        shiftind <- length(x@par)
        shift <- x@par[shiftind]
        nd <- length(d)
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        par <- c(a, d + shift)
        P <- P.poly(par=par, Theta=Theta, ot=offterm)
        ret <- .Call("dparsPoly", par, Theta, offterm, x@dat,
                     length(d), estHess)
        grad <- ret$grad
        hess <- ret$hess
        hess <- cbind(hess, rep(0, nrow(hess)))
        hess <- rbind(hess, rep(0, ncol(hess)))
        dc <- numeric(1L)
        Pfull <- P
        PQfull <- Pfull * (1-Pfull)
        P <- P.poly(c(a, d + shift), Theta, itemexp=TRUE, ot=offterm)
        rs <- dat
        for(i in 1L:ncol(rs))
            dc <- dc + rs[,i]/P[,i] * (PQfull[,i] - PQfull[,i+1L])
        dc <- sum(dc)
        grad <- c(grad, dc)
        if(estHess){
            cind <- ncol(hess)
            ddc <- ddd <- numeric(nrow(P))
            dda <- matrix(0, nrow(P), nfact)
            for(i in 1L:ncol(rs))
                ddc <- ddc + rs[,i]/P[,i]  * (Pfull[,i] - 3*Pfull[,i]^2 + 2*Pfull[,i]^3 -
                    Pfull[,i+1L] + 3*Pfull[,i+1L]^2 - 2*Pfull[,i+1L]^3) -
                    rs[,i]/P[,i]^2 * (PQfull[,i] - PQfull[,i+1L])^2
            hess[cind, cind] <- sum(ddc)
            for(i in 1L:nzetas)
                hess[cind, nfact + i] <- hess[nfact + i, cind] <-
                    sum((rs[,i]/P[,i] * (-Pfull[,i+1L] + 3*Pfull[,i+1L]^2 - 2*Pfull[,i+1L]^3) -
                    rs[,i]/P[,i]^2 * (PQfull[,i] - PQfull[,i+1L]) * (-PQfull[,i+1L]) +
                    rs[,i+1L]/P[,i+1L] * (Pfull[,i+1L] - 3*Pfull[,i+1L]^2 + 2*Pfull[,i+1L]^3) -
                    rs[,i+1L]/P[,i+1L]^2 * (PQfull[,i+1L] - PQfull[,i+2L]) * (PQfull[,i+1L])))
            for(j in 1L:nfact){
                tmp <- 0
                for(i in 1L:ncol(rs))
                        tmp <- tmp + (rs[,i]/P[,i] * Theta[,j] *
                                          (Pfull[,i] - 3*Pfull[,i]^2 + 2*Pfull[,i]^3 -
                                               Pfull[,i+1L] + 3*Pfull[,i+1L]^2 - 2*Pfull[,i+1L]^3) -
                                 rs[,i]/P[,i]^2 * (PQfull[,i] - PQfull[,i+1L]) * Theta[,j] *
                                      (PQfull[,i] - PQfull[,i+1L]))
                hess[cind, j] <- hess[j, cind] <- sum(tmp)
            }
        }
        ret <- list(grad=grad, hess=hess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- ExtractLambdas(x)
        P <- ProbTrace(x, Theta, itemexp = FALSE)
        grad <- hess <- vector('list', x@ncat)
        for(i in 1L:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1L:x@nfact){
            for(i in 1L:(ncol(P)-1L)){
                w1 <- P[,i] * (1-P[,i]) * a[j]
                w2 <- P[,i+1L] * (1-P[,i+1L]) * a[j]
                grad[[i]][ ,j] <- w1 - w2
                hess[[i]][ ,j] <- a[j]^2 * (2 * P[ ,i] * (1 - P[,i])^2 -
                                                P[ ,i] * (1 - P[,i]) -
                                                2 * P[ ,i+1L] * (1 - P[ ,i+1L])^2 +
                                                P[ ,i+1L] * (1 - P[ ,i+1L]))
            }
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_dP(x, Theta)
    }
)

# ----------------------------------------------------------------

setClass("gpcm", contains = 'AllItemsClass',
         representation = representation(mat='logical'))

setMethod(
    f = "print",
    signature = signature(x = 'gpcm'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'gpcm'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        par <- x@par
        d <- par[(length(par)-x@ncat+1L):length(par)]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        ns <- sum(grepl('ak', names(x@est)))
        par <- c(rlnorm(x@nfact, meanlog=-.2, sdlog=.5), rep(NA, ns), 0,
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.gpcm <- function(par, Theta, ot = 0, mat = FALSE, returnNum = FALSE)
{
    return(.Call("gpcmTraceLinePts", par, Theta, ot, FALSE, mat, returnNum))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.gpcm(x@par, Theta=Theta, ot=ot, mat=x@mat))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        ret <- .Call("dparsNominal", x, Theta, offterm, FALSE, estHess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- 0:(x@ncat - 1L)
        P <- P.nominal(c(a, ak, d), ncat=length(d), Theta=Theta)
        Num <- P.nominal(c(a, ak, d), ncat=length(d), Theta=Theta, returnNum = TRUE)
        Den <- rowSums(Num)
        grad <- hess <- vector('list', x@ncat)
        for(i in 1L:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1L:x@nfact){
            for(i in 1L:x@ncat){
                grad[[i]][ ,j] <- ak[i] * a[j] * P[ ,i] - P[ ,i] * (Num %*% (ak * a[j])) / Den
                hess[[i]][ ,j] <- ak[i]^2 * a[j]^2 * P[ ,i] -
                    2 * ak[i] * a[j] * P[,i] * (Num %*% (ak * a[j])) / Den +
                    2 * P[,i] * ((Num %*% (ak * a[j])) / Den)^2 -
                    P[,i] * ((Num %*% (ak^2 * a[j]^2)) / Den)
            }
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){
        nfact <- ncol(Theta)
        ncat <- x@ncat
        num <- P.gpcm(x@par, Theta=Theta, returnNum = TRUE, mat = x@mat)
        den <- rowSums(num)
        P <- num/den
        a <- x@par[1L:ncol(Theta)]
        if(!x@mat){
            ak <- matrix(x@par[(ncol(Theta)+1L):(ncol(Theta)+x@ncat)], nfact, ncat, byrow=TRUE)
        } else {
            ak <- matrix(x@par[(nfact+1L):(nfact+ncat*nfact)], nfact, ncat, byrow=TRUE)
        }
        dp <- matrix(0, nrow(Theta), length(x@par))
        aknum <- eak <- vector('list', nfact)
        e <- 0:(x@ncat-1)
        for(i in 1L:nfact){
            aknum[[i]] <- t(ak[i,] * t(num))
            eak[[i]] <- e*ak[i,]
        }
        aTheta <- as.numeric(a %*% t(Theta))
        for(i in 1L:nfact){
            for(j in 2L:ncat)
                dp[,i] <- dp[,i] + eak[[i]][j]*Theta[,i]*P[,j] -
                    e[j]*P[,j]*rowSums(Theta[,i] * aknum[[i]])
        }
        for(j in 1L:x@ncat){
            dp[,nfact + ncat + j] <- e[j] * P[,j] - e[j] * P[,j]^2 -
                as.numeric(e[-j] %*% t(P[,-j]*P[,j]))
        }
        return(dp)

    }
)

# ----------------------------------------------------------------

setClass("rsm", contains = 'AllItemsClass',
         representation = representation(mat='logical'))

setMethod(
    f = "print",
    signature = signature(x = 'rsm'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'rsm'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'rsm'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'rsm'),
    definition = function(x){
        par <- x@par
        d <- par[(length(par) - x@ncat):length(par)]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'rsm'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5), 0:(x@ncat-1), 0,
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE),
                 rnorm(1L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'rsm'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.rsm <- function(par, Theta, ot = 0, mat = FALSE, returnNum = FALSE)
{
    return(.Call("gpcmTraceLinePts", par, Theta, ot, TRUE, mat, returnNum))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.rsm(x@par, Theta=Theta, ot=ot))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        dat <- x@dat
        # nfact <- x@nfact
        # nzetas <- ncol(dat)
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        shift <- d[length(d)]
        dshift <- d <- d[-length(d)]
        dshift[-1L] <- d[-1L] + shift
        ak <- 0:(length(d)-1L)
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        # P <- ProbTrace(x=x, Theta=Theta, useDesign = FALSE, ot=offterm)
        tmp <- .Call("dparsNominal", x, Theta, offterm, TRUE, estHess)
        # num <- P.nominal(c(a, ak, dshift), ncat=length(ak), Theta=Theta, returnNum=TRUE, ot=offterm)
        grad <- tmp$grad
        hess <- tmp$hess

        #quick calcs for derivs
        # nfact <- length(a)
        # ncat <- length(d)
        # akind <- nfact
        # dind <- nfact + ncat*2 + 1L #go backwards
        # ak2 <- ak^2
        # P2 <- P^2
        # P3 <- P^3
        # aTheta <- as.vector(Theta %*% a)
        # aTheta2 <- aTheta^2
        # dat_num <- dat/num
        # numsum <- rowSums(num)
        # numD <- num %*% c(0, rep(1, ncol(num)-1L))
        # numak <- matrix(num %*% ak, nrow(Theta), ncol(Theta))
        # numakThetaD <- numak * Theta
        # numD2 <- num %*% c(0, rep(1, ncol(num)-1L))
        # numakThetaD2 <- numak * Theta
        # ak0 <- ak
        # ak0[1L] <- 0
        # cind <- length(grad)
        # if(estHess){
        #     tmp <- 0
        #     for(i in 1L:nzetas)
        #         tmp <- tmp + dat[,i]*numD^2 / numsum^2 - dat[,i]*numD2/numsum
        #     hess[cind, cind] <- sum(tmp)
        #     for(j in 1L:nzetas){
        #         tmp <- 0
        #         for(i in 1L:nzetas)
        #             tmp <- tmp + dat[,i]*P[,j]*numD/numsum - dat[,i]*P[,j]
        #         hess[cind, nfact+j] <- hess[nfact+j, cind] <- sum(tmp)
        #     }
        #     for(j in 1L:nfact){
        #         tmp <- 0
        #         for(i in 1L:nzetas)
        #             tmp <- tmp + dat[,i]*numD*numakThetaD[,j]/numsum^2 -
        #                 dat[,i]* (num %*% ak0*Theta[,j])/numsum
        #         hess[cind, j] <- hess[j, cind] <- sum(tmp)
        #     }
        # }
        ####
        #TODO - can't seem to get the last value of the gradient quite right for some reason....
        x2 <- x
        x2@est <- c(rep(FALSE, length(x2@est)-1L), TRUE)
        grad[x2@est] <- numerical_deriv(x@par[x2@est], EML, obj=x2, Theta=Theta,
                                        type='central')
        if(estHess && any(x@est)){
            hess[x@est, x@est] <- numerical_deriv(x@par[x@est], EML, obj=x,
                                                    Theta=Theta, type = 'Richardson',
                                                    gradient = FALSE)
        }
        ####
        ret <- list(grad=grad, hess=hess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        ret
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        t <- d[length(d)]
        d <- d[-length(d)]
        d[-1L] <- d[-1L] + t
        ak <- 0:(x@ncat - 1L)
        P <- P.nominal(c(a, ak, d), ncat=length(d), Theta=Theta)
        Num <- P.nominal(c(a, ak, d), ncat=length(d), Theta=Theta, returnNum = TRUE)
        Den <- rowSums(Num)
        grad <- hess <- vector('list', x@ncat)
        for(i in 1L:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1L:x@nfact){
            for(i in 1L:x@ncat){
                grad[[i]][ ,j] <- ak[i] * a[j] * P[ ,i] - P[ ,i] * (Num %*% (ak * a[j])) / Den
                hess[[i]][ ,j] <- ak[i]^2 * a[j]^2 * P[ ,i] -
                    2 * ak[i] * a[j] * P[,i] * (Num %*% (ak * a[j])) / Den +
                    2 * P[,i] * ((Num %*% (ak * a[j])) / Den)^2 -
                    P[,i] * ((Num %*% (ak^2 * a[j]^2)) / Den)
            }
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_dP(x, Theta)
    }
)

# ----------------------------------------------------------------

setClass("nominal", contains = 'AllItemsClass',
         representation = representation(mat='logical'))

setMethod(
    f = "print",
    signature = signature(x = 'nominal'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'nominal'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'nominal'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'nominal'),
    definition = function(x){
        d <- x@par[(length(x@par) - x@ncat + 1L):length(x@par)]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'nominal'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=-.2, sdlog=.5), 0,
                 abs(rnorm(x@ncat-1L, (x@ncat-1L) / 2, sd = 1)), 0,
                 rnorm(x@ncat-1L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'nominal'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x@est[(x@nfact+1L):(x@nfact + x@ncat)] <- FALSE
        x
    }
)

P.nominal <- function(par, ncat, Theta, returnNum = FALSE, ot = 0){
    return(.Call("nominalTraceLinePts", par, ncat, Theta, returnNum, ot))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nominal(par=x@par, ncat=x@ncat, Theta=Theta, ot=ot))
    }
)

nominalParDeriv <- function(a, ak, d, Theta, P, num, dat, estHess, gpcm = FALSE){
    nfact <- length(a)
    ncat <- length(d)
    akind <- nfact
    dind <- nfact + ncat
    ak2 <- ak^2
    P2 <- P^2
    P3 <- P^3
    aTheta <- as.vector(Theta %*% a)
    aTheta2 <- aTheta^2
    dat_num <- dat/num
    numsum <- rowSums(num)
    numakD <- num %*% ak
    numak2D2 <- num %*% ak2
    numakDTheta_numsum <- matrix(0, nrow(num), nfact)
    for(i in 1L:nfact)
        numakDTheta_numsum[,i] <- (num %*% ak * Theta[, i])/ numsum
    ret <- .Call('dparsNominal', a, ak, d, Theta, P, num, dat, nfact, ncat,
                 akind, dind, ak2, P2, P3, aTheta, aTheta2, dat_num, numsum, numakD,
                 numak2D2, numakDTheta_numsum, estHess)
    ret
}

setMethod(
    f = "Deriv",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        ret <- .Call("dparsNominal", x, Theta, offterm, FALSE, estHess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- x@par[(x@nfact+1):(x@nfact+x@ncat)]
        P <- P.nominal(c(a, ak, d), ncat=length(d), Theta=Theta)
        Num <- P.nominal(c(a, ak, d), ncat=length(d), Theta=Theta, returnNum = TRUE)
        Den <- rowSums(Num)
        grad <- hess <- vector('list', x@ncat)
        for(i in 1L:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1L:x@nfact){
            for(i in 1L:x@ncat){
                grad[[i]][ ,j] <- ak[i] * a[j] * P[ ,i] - P[ ,i] * (Num %*% (ak * a[j])) / Den
                hess[[i]][ ,j] <- ak[i]^2 * a[j]^2 * P[ ,i] -
                    2 * ak[i] * a[j] * P[,i] * (Num %*% (ak * a[j])) / Den +
                    2 * P[,i] * ((Num %*% (ak * a[j])) / Den)^2 -
                    P[,i] * ((Num %*% (ak^2 * a[j]^2)) / Den)
            }
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){
        num <- P.nominal(x@par, ncat=x@ncat, Theta=Theta, returnNum=TRUE)
        den <- rowSums(num)
        P <- num/den
        a <- x@par[1L:ncol(Theta)]
        ak <- x@par[(ncol(Theta)+1L):(ncol(Theta)+x@ncat)]
        dp <- matrix(0, nrow(Theta), length(x@par))
        aknum <- t(ak * t(num))
        aTheta <- as.numeric(a %*% t(Theta))
        nfact <- ncol(Theta)
        ncat <- x@ncat
        e <- 0:(x@ncat-1)
        eak <- e*ak
        for(i in 1L:nfact){
            for(j in 2L:ncat)
                dp[,i] <- dp[,i] + eak[j]*Theta[,i]*P[,j] - e[j]*P[,j]*rowSums(Theta[,i] * aknum)
        }
        for(j in 1L:x@ncat){
            dp[,nfact + ncat + j] <- e[j] * P[,j] - e[j] * P[,j]^2 - as.numeric(e[-j] %*% t(P[,-j]*P[,j]))
            dp[,nfact + j] <- dp[,nfact + ncat + j] * aTheta
        }
        return(dp)
    }
)

# ----------------------------------------------------------------

setClass("partcomp", contains = 'AllItemsClass')

setMethod(
    f = "print",
    signature = signature(x = 'partcomp'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'partcomp'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        d <- x@par[(x@nfact+1L):(length(x@par)-2L)]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 rlnorm(x@nfact, meanlog=.2, sdlog=.5),
                 rnorm(1L, -2, .5),
                 rnorm(1L, 2, .5))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.comp <- function(par, Theta, ot = 0)
{
    return(.Call('partcompTraceLinePts', par, Theta))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        return(P.comp(x@par, Theta=Theta))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        #local derivative from previous version with small mod
        #u and g in logit form
        dpars.comp <- function(lambda,zeta,g,r,f,Thetas,estHess)
        {
            nfact <- length(lambda)
            pars <- c(zeta,lambda,g)
            pgrad <- function(pars, r, thetas){
                nfact <- ncol(thetas)
                d <- pars[1L:nfact]
                a <- pars[(nfact+1L):(length(pars)-1L)]
                c <- pars[length(pars)]
                P <- P.comp(c(a,d,c,999), thetas)[,2L]
                Pstar <- P.comp(c(a,d,-999,999),thetas)[,2L]
                Qstar <- 1 - Pstar
                Q <- 1 - P
                c <- antilogit(c)
                g_1g <- c * (1 - c)
                const1 <- (r/P - (f-r)/Q)
                dd <- da <- rep(0,nfact)
                dc <- sum(r/P * (g_1g * (1 - Pstar)) + (f-r)/Q * (g_1g * (Pstar - 1)))
                for(i in 1L:nfact){
                    Pk <- P.mirt(c(a[i],d[i],-999,999),matrix(thetas[,i]))[,2L]
                    Qk <- 1 - Pk
                    dd[i] <- sum((1-c)*Pstar*Qk*const1)
                    da[i] <- sum((1-c)*Pstar*Qk*thetas[,i]*const1)
                }
                return(c(dd,da,dc))
            }
            phess <- function(pars, r, thetas){
                nfact <- ncol(thetas)
                d <- pars[1L:nfact]
                a <- pars[(nfact+1L):(length(pars)-1L)]
                c <- pars[length(pars)]
                P <- P.comp(c(a,d,c,999), thetas)[,2L]
                Pstar <- P.comp(c(a,d,-999,999),thetas)[,2L]
                Qstar <- 1 - Pstar
                Q <- 1 - P
                g <- c <- antilogit(c)
                g_1g <- c * (1 - c)
                const1 <- (r/P - (f-r)/Q)
                const2 <- (r/P^2 + (f-r)/Q^2)
                hess <- matrix(0,nfact*2+1,nfact*2+1)
                dNames <- paste("d",1:nfact,sep='_')
                aNames <- paste("a",1:nfact,sep='_')
                Names <- c(paste("d",1:nfact,sep='_'),paste("a",1:nfact,sep='_'),'c_0')
                for(i in 1L:(nfact*2+1L)){
                    for(j in 1L:(nfact*2+1L)){
                        if(i <= j){
                            d1 <- strsplit(Names[c(i,j)],"_")[[1L]]
                            d2 <- strsplit(Names[c(i,j)],"_")[[2L]]
                            k <- as.numeric(d1[2L])
                            m <- as.numeric(d2[2L])
                            Pk <- P.mirt(c(a[k],d[k],-999,999),matrix(thetas[,k]))[,2L]
                            Qk <- 1 - Pk
                            Pm <- P.mirt(c(a[m],d[m],-999,999),matrix(thetas[,m]))[,2L]
                            Qm <- 1 - Pm
                            if(i == j && d1[1L] == 'd'){
                                hess[i,i] <- sum(r/P * ((1-g)*(Pstar - 3*Pstar*Pk + 2*Pstar*Pk^2) ) -
                                                     r/P^2 * ((1-g) * (Pstar - Pstar*Pk))^2 +
                                                     (f-r)/Q * ((1-g)*(-Pstar + 3*Pstar*Pk - 2*Pstar*Pk^2)) -
                                                     (f-r)/Q^2 * ((1-g)*(-Pstar + Pstar*Pk))^2)
                                next
                            }
                            if(i == j && d1[1L] == 'a'){
                                hess[i,i] <- sum(r/P * ((1-g)*thetas[,k]^2*(Pstar - 3*Pstar*Pk + 2*Pstar*Pk^2) ) -
                                        r/P^2 * ((1-g)*thetas[,k] * (Pstar - Pstar*Pk))^2 +
                                        (f-r)/Q * ((1-g)*thetas[,k]^2*(-Pstar + 3*Pstar*Pk - 2*Pstar*Pk^2)) -
                                        (f-r)/Q^2 * ((1-g)*thetas[,k] * (-Pstar + Pstar*Pk))^2)
                                next
                            }
                            if(i == j && d1[1L] == 'c'){
                                hess[i,i] <- sum(r/P * (g_1g * (2.0*(1-g) - 1.0 - 2.0*(1-g)*Pstar + Pstar)) -
                                        r/P^2 * (g_1g * (1.0 - Pstar)) * (g_1g * (1.0 - Pstar)) +
                                        (f-r)/Q * (g_1g * (-2.0*(1-g) + 1.0 + 2.0*(1-g)*Pstar - Pstar)) -
                                        (f-r)/Q^2 * (g_1g * (-1.0 + Pstar)) * (g_1g * (-1.0 + Pstar)))
                                next
                            }
                            if(d1[1L] == 'a' && d2[1L] == 'a'){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*thetas[,m]*
                                                                  Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2))
                                next
                            }
                            if(d1[1L] == 'd' && d2[1L] == 'd'){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2))
                                next
                            }
                            if(d1[1L] == 'a' && d2[1L] == 'c'){
                                hess[i,j] <- hess[j,i] <- sum(r/P*(g_1g*thetas[,k]*(-Pstar + Pstar*Pk)) -
                                                                  r/P^2*((1-g)*thetas[,k]*(Pstar - Pstar*Pk)*(g_1g*(1 - Pstar))) +
                                                                 (f-r)/Q * (g_1g*thetas[,k]*(Pstar - Pstar*Pk))-
                                                                  (f-r)/Q^2 * ((1-g)*thetas[,k]*(-Pstar + Pstar*Pk)*(g_1g*(-1 + Pstar))))
                                next
                            }
                            if(d1[1L] == 'd' && d2[1L] == 'c'){
                                hess[i,j] <- hess[j,i] <- sum(r/P*(g_1g*(-Pstar + Pstar*Pk)) -
                                                                  r/P^2*((1-g)*(Pstar - Pstar*Pk)*(g_1g*(1 - Pstar))) +
                                                                  (f-r)/Q * (g_1g*(Pstar - Pstar*Pk))-
                                                                  (f-r)/Q^2 * ((1-g)*(-Pstar + Pstar*Pk)*(g_1g*(-1 + Pstar))))
                                next
                            }
                            if(d1[1L] == 'd' && d2[1L] == 'a' && d1[2] == d2[2]){
                                hess[i,j] <- hess[j,i] <- sum(
                                        r/P * ((1-g)*thetas[,k]*(Pstar - 3*Pstar*Pk + 2*Pstar*Pk^2) ) -
                                        r/P^2 * ((1-g)*thetas[,k] * (Pstar - Pstar*Pk) * (1-g)*(Pstar - Pstar*Pk)) +
                                        (f-r)/Q * ((1-g)*thetas[,k]*(-Pstar + 3*Pstar*Pk - 2*Pstar*Pk^2)) -
                                        (f-r)/Q^2 * ((1-g)*thetas[,k] * (-Pstar + Pstar*Pk) * (1-g)*(-Pstar + Pstar*Pk)))
                                next
                            }
                            if(d1[1L] == 'd' && d2[1L] == 'a' && d1[2] != d2[2]){
                                hess[i,j] <- hess[j,i] <- sum(r/P * ((1-g)*thetas[,m]*(Pstar - Pstar*Pk - Pstar*Pm + Pstar^2)) -
                                            r/P^2 * ((1-g)*thetas[,m] * (Pstar - Pstar*Pm) * (1-g)*(Pstar - Pstar*Pk))  +
                                            (f-r)/Q * ((1-g)*thetas[,m]*(-Pstar + Pstar*Pk + Pstar*Pm - Pstar^2)) -
                                            (f-r)/Q^2 * ((1-g)*thetas[,m] * (-Pstar + Pstar*Pm) * (1-g)*(-Pstar + Pstar*Pk)))
                                next
                            }
                        }
                    }
                }
                return(hess)
            }
            #old pars in the form d, a, g
            g <- pgrad(pars, r, Thetas)
            if(estHess) h <- phess(pars, r, Thetas)
            else h <- matrix(0, length(g), length(g))

            #translate into current version
            grad <- c(g[(nfact+1L):(nfact*2)], g[1L:nfact], g[length(g)], 0)
            hess <- matrix(0, ncol(h) + 1L, ncol(h) + 1L)
            if(estHess){
                hess[1L:nfact, 1L:nfact] <- h[(nfact+1L):(nfact*2),(nfact+1L):(nfact*2)] #a block
                hess[(nfact+1L):(nfact*2),(nfact+1L):(nfact*2)] <- h[1L:nfact, 1L:nfact] #d block
                hess[nfact*2 + 1L, nfact*2 + 1L] <- h[nfact*2 + 1L, nfact*2 + 1L] #g
                hess[nfact*2 + 1L, 1L:nfact] <- hess[1:nfact, nfact*2 + 1L] <-
                    h[nfact*2 + 1L, (nfact+1L):(nfact*2)] #ga
                hess[nfact*2 + 1L, (nfact+1L):(nfact*2)] <- hess[(nfact+1L):(nfact*2), nfact*2 + 1L] <-
                    h[nfact*2 + 1L, 1L:nfact] #gd
                hess[(nfact+1L):(nfact*2), 1L:nfact] <- t(h[(nfact+1L):(nfact*2), 1L:nfact])
                hess[1L:nfact, (nfact+1L):(nfact*2)] <- t(h[1L:nfact, (nfact+1L):(nfact*2)]) #ads
            }

            return(list(grad=grad, hess=hess))
        }
        #####
        f <- rowSums(x@dat)
        r <- x@dat[ ,2L]
        nfact <- x@nfact
        a <- x@par[1L:nfact]
        d <- x@par[(nfact+1L):(nfact*2L)]
        g <- x@par[length(x@par)-1L]
        tmp <- dpars.comp(lambda=ExtractLambdas(x),zeta=ExtractZetas(x),g=x@par[nfact*2L + 1L],r=r, f=f,
                          Thetas=Theta, estHess=estHess)
        ret <- list(grad=tmp$grad, hess=tmp$hess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){
        N <- nrow(Theta)
        nfact <- ncol(Theta)
        parlength <- length(x@par)
        u <- x@par[parlength]
        g <- x@par[parlength - 1L]
        d <- ExtractZetas(x)
        a <- ExtractLambdas(x)
        Pstar <- P.comp(c(a, d, -999, 999), Theta)[,2L]
        g <- antilogit(g)
        u <- antilogit(u)
        grad <- hess <- vector('list', 2L)
        grad[[1L]] <- grad[[2L]] <- hess[[1L]] <- hess[[2L]] <- matrix(0, N, nfact)
        for(j in 1L:nfact){
            Pn <- P.mirt(c(a[j], d[j], -999, 999), Theta)[,2L]
            grad[[2L]][ ,j] <- (u-g)*Pstar*a[j]*(1 - Pn)
            grad[[1L]][ ,j] <- -1*grad[[2L]][ ,j]
            hess[[2L]][ ,j] <- (u-g)*a[j]^2*Pstar*(1 - 3*Pn + 2*Pn^2)
            hess[[1L]][ ,j] <- -1*hess[[2L]][ ,j]
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){
        message('partcomp derivatives not optimized') ##TODO
        numDeriv_dP(x, Theta)
    }
)

# ----------------------------------------------------------------

setClass("nestlogit", contains = 'AllItemsClass',
         representation = representation(correctcat='integer'))

setMethod(
    f = "print",
    signature = signature(x = 'nestlogit'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'nestlogit'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'nestlogit'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'nestlogit'),
    definition = function(x){
    	stop('not written')
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'nestlogit'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 rnorm(1L),
                 rnorm(1L, -2, .5),
                 rnorm(1L, 2, .5), 0,
                 abs(rnorm(x@ncat-2L, (x@ncat-2L) / 2, sd = 1)), x@ncat,
                 rnorm(x@ncat-2L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'nestlogit'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x@est[(x@nfact+5L):(x@nfact + x@ncat + 1L)] <- FALSE
        x
    }
)

P.nestlogit <- function(par, Theta, correct, ncat)
{
    return(.Call("nestlogitTraceLinePts", par, Theta, correct, ncat))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        return(P.nestlogit(x@par, Theta=Theta, correct=x@correctcat, ncat=x@ncat))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        #u and g in logit
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        dat <- x@dat
        if(estHess && any(x@est))
            hess[x@est, x@est] <- numerical_deriv(x@par[x@est], EML, obj=x, Theta=Theta,
                                                  gradient = FALSE, type = 'Richardson')
        nfact <- x@nfact
        a <- x@par[1L:x@nfact]
        d <- x@par[x@nfact+1L]
        g <- x@par[x@nfact+2L]
        u <- x@par[x@nfact+3L]
        ak <- x@par[(x@nfact+4L):(x@nfact+4L+x@ncat-2L)]
        dk <- x@par[(length(x@par)-length(ak)+1):length(x@par)]
        correct <- x@correctcat
        Pd <- P.mirt(c(a, d, g, u), Theta=Theta)[,2L]
        g <- antilogit(g)
        u <- antilogit(u)
        Qd <- 1 - Pd
        Pstar <- P.mirt(c(a, d, -999, 999), Theta=Theta)[,2L]
        Qstar <- 1 - Pstar
        num <- P.nominal(c(rep(1, nfact), ak, dk), ncat=length(dk), Theta=Theta, returnNum=TRUE)
        den <- rowSums(num)
        Pn <- num/den
        cdat <- dat[,correct]
        idat <- dat[,-correct]
        nd <- ncol(idat)
        g_1g <- g * (1 - g)
        u_1u <- u * (1 - u)
        for(i in 1L:nfact)
            grad[i] <- sum( (u-g) * Theta[,i] * Qstar * Pstar * (
                cdat / Pd - rowSums(idat/Qd)) )
        grad[nfact+1L] <- sum( (u-g) * Qstar * Pstar * (
                cdat / Pd - rowSums(idat/Qd)) )
        grad[nfact+2L] <- sum( ((cdat * g_1g * (1-Pstar)/Pd) + rowSums(idat * g_1g * (Pstar - 1)/Qd)) )
        grad[nfact+3L] <- sum( (cdat * u_1u * Pstar / Pd - rowSums(idat * u_1u * Pstar / Qd) ))
        for(j in 1L:nd){
            grad[nfact+3L+j] <- sum((
                (idat[,j] * Qd * rowSums(Theta) * (Pn[,j] - Pn[,j]^2) * den) / (Qd * num[,j]) -
                    rowSums(idat[,-j, drop=FALSE]) * rowSums(Theta) * Pn[,j]))
            grad[nfact+3L+nd+j] <- sum((
                (idat[,j] * Qd * (Pn[,j] - Pn[,j]^2) * den) / (Qd * num[,j]) -
                    rowSums(idat[,-j, drop=FALSE]) * Pn[,j]))
        }
        ret <- list(grad=grad, hess=hess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- x@par[1:x@nfact]
        d <- x@par[x@nfact+1L]
        g <- x@par[x@nfact+2L]
        u <- x@par[x@nfact+3L]
        ak <- x@par[(x@nfact+4L):(x@nfact+4L+x@ncat-2L)]
        dk <- x@par[(length(x@par)-length(ak)+1):length(x@par)]
        Pn <- P.nominal(c(rep(1,ncol(Theta)), ak, dk), ncat=length(dk), Theta=Theta)
        Num <- P.nominal(c(rep(1,ncol(Theta)), ak, dk), ncat=length(dk), Theta=Theta, returnNum = TRUE)
        Den <- rowSums(Num)
        Pstar <- P.mirt(c(a, d, -999, 999), Theta)[,2L]
        Q <- P.mirt(c(a, d, g, u), Theta)[,1]
        g <- antilogit(g)
        u <- antilogit(u)
        Num2 <- P <- matrix(0, nrow(Theta), x@ncat)
        P[,-x@correctcat] <- Pn
        Num2[,-x@correctcat] <- Num
        Num <- Num2
        ak2 <- dk2 <- numeric(x@ncat)
        ak2[-x@correctcat] <- ak
        dk2[-x@correctcat] <- dk
        ak <- ak2
        dk <- dk2
        grad <- hess <- vector('list', x@ncat)
        for(i in 1L:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1L:x@nfact){
            for(i in 1L:x@ncat){
                if(i == x@correctcat){
                    grad[[i]][ ,j] <- (u-g) * a[j] * (Pstar * (1 - Pstar))
                    hess[[i]][ ,j] <- 2 * (u - g) * a[j]^2 * ((1 - Pstar)^2 * Pstar) -
                        (u - g) * a[j]^2 * (Pstar * (1 - Pstar))
                } else {
                    grad[[i]][ ,j] <- -(u-g) * a[j] * (Pstar * (1 - Pstar)) * P[,i] +
                        Q * (ak[i] * P[ ,i] - P[ ,i] * (Num %*% (ak)) / Den)
                    hess[[i]][ ,j] <-
                        -2 * (u - g) * a[j]^2 * (1 - Pstar)^2 * Pstar * P[,i] +
                        (u - g) * a[j]^2 * Pstar * (1 - Pstar) * P[,i] -
                        2 * (u - g) * a[j] * ak[i] * (1 - Pstar) * Pstar * P[,i] +
                        2 * a[j] *  (Pstar * (1 - Pstar)) * P[,i] * (Num %*% (ak)) / Den +
                        Q * ak[i]^2 * P[ ,i] -
                        2 * Q * ak[i] * P[,i] * (Num %*% (ak)) / Den +
                        2 * Q * P[,i] * ((Num %*% (ak)) / Den)^2 -
                        Q * P[,i] * ((Num %*% (ak^2)) / Den)
                }
            }
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta){
        message('nestlogit derivatives not optimized') ##TODO
        numDeriv_dP(x, Theta)
    }
)

# ----------------------------------------------------------------

setClass("ideal", contains = 'AllItemsClass')

setMethod(
    f = "print",
    signature = signature(x = 'ideal'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'ideal'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'ideal'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'ideal'),
    definition = function(x){
        d <- x@par[length(x@par)]
        d
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'ideal'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 -rlnorm(1, meanlog=0, sdlog=.5))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'ideal'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.ideal <- function(par, Theta, ot = 0)
{
    alpha <- par[length(par)]
    beta <- par[-length(par)]
    eta <- -0.5*(Theta %*% beta + alpha)^2
    eta <- ifelse(eta < -20, -20, eta)
    P <- exp(eta)
    P <- cbind((1-P), P)
    P <- ifelse(P < 1e-20, 1e-20, P)
    P <- ifelse(P > (1 - 1e-20), (1 - 1e-20), P)
    return(P)
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'ideal', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        return(P.ideal(x@par, Theta=Theta))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'ideal', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        d <- x@par[length(x@par)]
        a <- x@par[-length(x@par)]
        P <- ProbTrace(x, Theta=Theta)[,2L]
        Q <- (1 - P)
        Q <- ifelse(Q < 1e-7, 1e-7, Q)
        int <- as.numeric(Theta %*% a + d)
        for(i in 1L:ncol(Theta))
            grad[i] <- -sum( x@dat[,1] * int * Theta[,i] * -P / Q +
                               x@dat[,2] * int * Theta[,i])
        grad[i+1L] <- -sum(2 * x@dat[,1] * int * -P / Q +
                           2 * x@dat[,2] * int)/2
        if(estHess && any(x@est))
            hess[x@est, x@est] <- numerical_deriv(x@par[x@est], EML, obj=x,
                                                  Theta=Theta, gradient=FALSE, type = 'Richardson')
        return(list(grad = grad, hess=hess))
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'ideal', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_DerivTheta(x, Theta)
#         N <- nrow(Theta)
#         nfact <- ncol(Theta)
#         P <- ProbTrace(x, Theta=Theta)[,2L]
#         d <- x@par[length(x@par)]
#         a <- x@par[-length(x@par)]
#         int <- as.numeric(t(a %*% t(Theta)) + d)
#         grad <- hess <- vector('list', 2L)
#         grad[[1L]] <- grad[[2L]] <- hess[[1L]] <- hess[[2L]] <- matrix(0, N, nfact)
#         for(i in 1L:nfact){
#             grad[[2L]][ ,i] <- 2 * a[i] * int * P
#             grad[[1L]][ ,i] <- -1 * grad[[2L]][ ,i]
#             hess[[2L]][ ,i] <- 2 * a[i]^2 * int * P + 4 * int^2 * a[i]^2 * P
#             hess[[1L]][ ,i] <- -1 * hess[[2L]][ ,i]
#         }
#         return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'ideal', Theta = 'matrix'),
    definition = function(x, Theta){
        P <- ProbTrace(x, Theta=Theta)[,2L]
        dp <- matrix(0, nrow(Theta), length(x@par))
        d <- x@par[length(x@par)]
        a <- x@par[-length(x@par)]
        int <- as.numeric(t(a %*% t(Theta)) + d)
        for(i in 1L:ncol(Theta))
            dp[,i] <- -2 * Theta[,i] * int * P
        dp[,i+1L] <- -2 * int * P
        dp
    }
)

# ----------------------------------------------------------------

setClass("lca", contains = 'AllItemsClass',
         representation = representation(score='numeric'))

setMethod(
    f = "print",
    signature = signature(x = 'lca'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'lca'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'lca'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'lca'),
    definition = function(x){
    	stop('not written')
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'lca'),
    definition = function(x){
        par <- rnorm(length(x@par))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'lca'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

P.lca <- function(par, Theta, ncat, returnNum = FALSE){
    return(.Call('lcaTraceLinePts', par, Theta, ncat, returnNum))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'lca', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        P <- P.lca(x@par, Theta=Theta, ncat=x@ncat)
        return(P)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'lca', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        ret <- .Call('dparslca', x@par, Theta, FALSE, x@dat, offterm) #TODO change FALSE to estHess
        if(estHess && any(x@est)){
            ret$hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x,
                                                         Theta=Theta)
        }
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'lca', Theta = 'matrix'),
    definition = function(x, Theta){
        stop('not written')
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'lca', Theta = 'matrix'),
    definition = function(x, Theta){
        P <- ProbTrace(x, Theta)
        dp <- matrix(0, nrow(Theta), length(x@par))
        ind <- 1L
        for(j in 2L:x@ncat){
            for(i in 1:ncol(Theta)){
                dp[,ind] <- Theta[,i] * (P[,j] -
                        rowSums(P[,j,drop=FALSE] * P[,j,drop=FALSE]))
                ind <- ind + 1L
            }
        }
        dp
    }
)

# ----------------------------------------------------------------

setClass("spline", contains = 'AllItemsClass',
         representation = representation(stype='character',
                                         Theta_prime='matrix',
                                         sargs='list'))

setMethod(
    f = "print",
    signature = signature(x = 'spline'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'spline'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'spline'),
    definition = function(x){
        numeric(x@nfact)
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'spline'),
    definition = function(x){
        stop('not written')
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'spline'),
    definition = function(x){
        par <- rnorm(length(x@par), sd = 20)
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'spline'),
    definition = function(x){
        stop('spline null should not be run')
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'spline', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(Theta) != nrow(x@Theta_prime))
            x <- loadSplineParsItem(x, Theta)
        P <- P.lca(x@par, Theta=x@Theta_prime, ncat=x@ncat)
        return(P)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'spline', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        ret <- .Call('dparslca', x@par, x@Theta_prime, FALSE, x@dat, offterm)
        if(estHess && any(x@est)){
            ret$hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x,
                                                        Theta=x@Theta_prime)
        }
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'spline', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_DerivTheta(x, Theta)
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'spline', Theta = 'matrix'),
    definition = function(x, Theta){
        P <- ProbTrace(x, Theta)
        dp <- matrix(0, nrow(Theta), length(x@par))
        ind <- 1L
        for(j in 2L:x@ncat){
            for(i in 1:ncol(Theta)){
                dp[,ind] <- Theta[,i] * (P[,j] -
                                             rowSums(P[,j,drop=FALSE] * P[,j,drop=FALSE]))
                ind <- ind + 1L
            }
        }
        dp
    }
)

# ----------------------------------------------------------------

setClass('custom', contains = 'AllItemsClass',
         representation = representation(name='character',
                                         P='function',
                                         gr='function',
                                         usegr='logical',
                                         hss='function',
                                         gen='function',
                                         usehss='logical',
                                         userdata='matrix',
                                         useuserdata='logical',
                                         derivType='character'))

setMethod(
    f = "print",
    signature = signature(x = 'custom'),
    definition = function(x, ...){
        cat('Custom item object named:', x@name)
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'custom'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'custom'),
    definition = function(x){
        a <- rep(.001, x@nfact)
        a
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'custom'),
    definition = function(x){
    	stop('not written')
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'custom'),
    definition = function(x){
        x@par <- x@gen(x)
        x
    }
)

setMethod(
    f = "set_null_model",
    signature = signature(x = 'custom'),
    definition = function(x){
        stop('calculating null models for custom item types is ambiguous.', call. = FALSE)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(x@useuserdata) Theta <- cbind(Theta, x@userdata)
        return(x@P(x@par, Theta=Theta, ncat=x@ncat))
    }
)

setMethod("initialize",
          'custom',
          function(.Object, name, par, est, lbound, ubound, P, gr, hss, gen, userdata, derivType) {
              dummyfun <- function(...) return(NULL)
              names(est) <- names(par)
              usegr <- usehss <- useuserdata <- TRUE
              .Object@name <- name
              .Object@par <- par
              .Object@est <- est
              .Object@P <- P
              .Object@derivType <- derivType
              if(is.null(gr)){
                  .Object@gr <- dummyfun
                  usegr <- FALSE
              } else .Object@gr <- gr
              if(is.null(hss)){
                  .Object@hss <- dummyfun
                  usehss <- FALSE
              } else .Object@hss <- hss
              if(is.null(gen)){
                  .Object@gen <- function(object) object@par
              } else .Object@gen <- gen
              if(is.null(userdata)){
                  .Object@userdata <- matrix(NaN)
                  useuserdata <- FALSE
              } else .Object@userdata <- userdata
              .Object@usegr <- usegr
              .Object@usehss <- usehss
              .Object@useuserdata <- useuserdata
              .Object@lbound <- if(!is.null(lbound)) lbound  else rep(-Inf, length(par))
              .Object@ubound <- if(!is.null(ubound)) ubound  else rep(Inf, length(par))
              .Object
          }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        if(x@useuserdata) Theta <- cbind(Theta, x@userdata)
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(x@usegr) grad <- x@gr(x, Theta)
        else grad[x@est] <- numerical_deriv(x@par[x@est], EML, obj=x, Theta=Theta, type=x@derivType)
        if(estHess){
            if(x@usehss) hess <- x@hss(x, Theta)
            else hess[x@est, x@est] <- numerical_deriv(x@par[x@est], EML, obj=x,
                                                       Theta=Theta, type=x@derivType, gradient=FALSE)
        }
        return(list(grad = grad, hess=hess))
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_DerivTheta(x, Theta)
    }
)

setMethod(
    f = "dP",
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_dP(x, Theta)
    }
)

# ----------------------------------------------------------------

# itemtype='grsmIRT', 1-dim graded rating scale
# slope, thresholds and 1 location parameters
# ExtractZetas, GenRandomPar, ProbTrace, initialize MODIFIED


setClass("grsmIRT", contains = 'AllItemsClass',
         representation = representation())

setMethod(
  f = "print",
  signature = signature(x = 'grsmIRT'),
  definition = function(x, ...){
    cat('Item object of class:', class(x))
  }
)

setMethod(
  f = "show",
  signature = signature(object = 'grsmIRT'),
  definition = function(object){
    print(object)
  }
)

#extract the slopes (should be a vector of length nfact)
setMethod(
  f = "ExtractLambdas",
  signature = signature(x = 'grsmIRT'),
  definition = function(x){
    x@par[1L:x@nfact] #slopes
  }
)

#extract the intercepts
setMethod(
  f = "ExtractZetas",
  signature = signature(x = 'grsmIRT'),
  definition = function(x){
    #x@par[(x@nfact+1):(length(x@par))] #intercepts
    x@par[(x@nfact+1):(x@ncat+1)] #intercepts
    # if we set :(length(x@par)), when collapsed, may produce error?
    # MAYBE DOESNT MATTER?
  }
)

# generating random starting values (only called when, e.g., mirt(..., GenRandomPars = TRUE))
setMethod(
  f = "GenRandomPars",
  signature = signature(x = 'grsmIRT'),
  definition = function(x){
    par <- c(rlnorm(1,0,1),sort(rnorm(x@ncat-1, 0, 1), decreasing=TRUE), rnorm(1,0,0.5))
    x@par[x@est] <- par[x@est]
    x
  }
)

# how to set the null model to compute statistics like CFI and TLI (usually just fixing slopes to 0)
setMethod(
  f = "set_null_model",
  signature = signature(x = 'grsmIRT'),
  definition = function(x){
    x@par[1L:(x@nfact+x@ncat)] <- 0
    x@est[1L:(x@nfact+x@ncat)] <- FALSE
    x
  }
)

# probability trace line function. Must return a matrix with a trace line for each category
setMethod(
  f = "ProbTrace",
  signature = signature(x = 'grsmIRT', Theta = 'matrix'),
  definition = function(x, Theta){


    th1 = Theta[,1];

    ncat = x@ncat
    a1 = x@par[1]
    d = x@par[2:ncat]
    c = x@par[ncat+1]
    nr = nrow(Theta); nc=ncat-1

    if (all(d == sort(d, decreasing=TRUE))) {

      D.star = matrix(c+d, nrow=nr, ncol=nc, byrow=TRUE)
      TH1 = matrix(th1, nrow=nr, ncol=nc)
      A = matrix(a1, nrow=nr, ncol=nc)
      P.star = 1/(1+exp(-A*(D.star+TH1)))

      P=cbind(1, P.star)-cbind(P.star, 0)

      P <- ifelse(P < 1e-20, 1e-20, P)
      P <- ifelse(P > (1 - 1e-20), (1 - 1e-20), P)
    } else {

      P <- matrix(1e-20, nrow=nr, ncol=ncat)}

    return(P)

  }
)

# complete-data derivative used in parameter estimation
# Analytic Derivation is provided...
# and it fails, numerical derivation is used.

setMethod(
  f = "Deriv",
  signature = signature(x = 'grsmIRT', Theta = 'matrix'),
  definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
    grad <- rep(0, length(x@par))
    hess <- matrix(0, length(x@par), length(x@par))

    dat <- x@dat;
    th1 <- Theta[,1];

    ncat <- x@ncat
    a1 <- x@par[1]
    d <- x@par[2:ncat]

    c <- x@par[ncat+1]
    nr <- length(th1);
    nc <- ncat - 1

    beta.mean <- mean(d)
    beta.mat <- matrix(d, nrow=nr, ncol=nc, byrow=TRUE)

    Sigma.mat <- (beta.mat - beta.mean)
    th1.mat <- matrix(th1, nr, nc)
    B.mat <- th1.mat + c + beta.mat

    lam.star.mat <- exp(-a1*B.mat)
    lam.mat <- cbind(lam.star.mat, 1)-cbind(0, lam.star.mat)

    P.star <- 1/(1+lam.star.mat)
    Q.star <- 1- P.star
    P <- P.star[,-nc] - P.star[,-1]
    inv.P <- 1/P

    #when a1*B.mat >30, inv.P would be overflow...

    PQ.star <- P.star*Q.star

    P.star.a1 <- PQ.star*B.mat
    P.star.beta <- list()

    for (i in 1:(ncat-1))
      P.star.beta[[i]] <- matrix(0, nr, nc)

    for (i in 1:(ncat-1))
      P.star.beta[[i]][,i] <- (PQ.star*a1)[,i]

    P.star.c <- PQ.star*a1

    P.a1 <- matrix(NA, nr, ncat);
    P.c <- matrix(NA, nr, ncat);
    P.beta <- list(); for (i in 1:(ncat-1)) P.beta[[i]] <- matrix(NA, nr, ncat);

    P.a1[, 2:nc] <- inv.P*(P.star.a1[,-nc] - P.star.a1[, -1])

    for (i in 1:(ncat-1))
      P.beta[[i]][, 2:nc] <- inv.P*(P.star.beta[[i]][, -nc] - P.star.beta[[i]][, -1])

    P.c[, 2:nc] <- inv.P * (P.star.c[, -nc] - P.star.c[, -1])

    # the first response category
    P.a1[,1] <- -P.star[,1]*B.mat[,1]
    for (i in 2:(ncat-1))
      P.beta[[i]][,1] <- 0

    P.beta[[1]][,1] <- -P.star[,1]*a1
    P.c[,1] <- -P.star[,1]*a1

    # the last response category
    P.a1[,ncat]=Q.star[,ncat-1]*B.mat[,ncat-1]

    for (i in 1:(ncat-2))
      P.beta[[i]][,ncat] <- 0

    P.beta[[ncat-1]][,ncat] <- Q.star[,ncat-1]*a1
    P.c[,ncat] <- Q.star[,ncat-1]*a1


    grad[1] <- sum(P.a1 * dat)
    for (i in 1:(ncat-1))
      grad[i+1] <- sum(P.beta[[i]] *dat)
    grad[ncat+1] <- sum(P.c * dat)

    grad[!x@est] <- 0

    Nest <- is.nan(grad) | is.infinite(grad)

    if (sum(Nest)>0) {
      grad <- rep(0, length(x@par))
      grad[x@est & Nest] <- numDeriv::grad(EML, x@par[x@est & Nest], obj=x, Theta=Theta)
    }

    if(estHess && any(x@est))
      hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x,
                                              Theta=Theta)
    ret <- list(grad=grad, hess=hess)
    if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
    return(ret)
  }
)


# derivative of the model wft to the Theta values (done numerically here)
setMethod(
  f = "DerivTheta",
  signature = signature(x = 'grsmIRT', Theta = 'matrix'),
  definition = function(x, Theta){
    numDeriv_DerivTheta(x, Theta) #replace with analytical derivatives
  }
)

# derivative of the probability trace line function wrt Theta (done numerically here)
setMethod(
  f = "dP",
  signature = signature(x = 'grsmIRT', Theta = 'matrix'),
  definition = function(x, Theta){
    numDeriv_dP(x, Theta) #replace with analytical derivatives
  }
)

# defines how the item should be initiallized with the S4 function new()
#   before the parameter estimates are updated
#   (set starting values, lower/upper bounds, logicals indicating estimation, etc)
setMethod("initialize",
          'grsmIRT',
          function(.Object, nfact, ncat){
            if(nfact != 1L)
                stop('grsmIRT only possible for unidimensional models')
            stopifnot(ncat >= 2L)
            .Object@par <- c(rep(1, nfact),  seq(1, -1, length.out=ncat-1), 0)
            #.Object@par <- c(rep(1, nfact),  seq(-3, 3, length.out=ncat-1), 0)
            # -3 ~ 3 seems to be too far away
            names(.Object@par) = c(paste("a",1:nfact, sep=""),
                                   paste("b", 1:(ncat-1), sep=""), "c")
            .Object@est <- c(rep(TRUE, nfact), rep(TRUE,ncat-1), TRUE)
            .Object@lbound <- rep(-Inf, nfact+ncat)
            .Object@ubound <- rep(Inf, nfact+ncat)
            .Object
          }
)


# ----------------------------------------------------------------

# experimental itemtype (used as a template to create custom IRT models). Edit this if you want
#  to experiment with your own customized IRT models

setClass("experimental", contains = 'AllItemsClass',
         representation = representation())

setMethod(
    f = "print",
    signature = signature(x = 'experimental'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'experimental'),
    definition = function(object){
        print(object)
    }
)

#extract the slopes (should be a vector of length nfact)
setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'experimental'),
    definition = function(x){
        x@par[1L:x@nfact] #slopes
    }
)

#extract the intercepts
setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'experimental'),
    definition = function(x){
        x@par[length(x@par)] #intercepts
    }
)

# generating random starting values (only called when, e.g., mirt(..., GenRandomPars = TRUE))
setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'experimental'),
    definition = function(x){
        par <- c(rlnorm(1, .2, .2), rnorm(1))
        x@par[x@est] <- par[x@est]
        x
    }
)

# how to set the null model to compute statistics like CFI and TLI (usually just fixing slopes to 0)
setMethod(
    f = "set_null_model",
    signature = signature(x = 'experimental'),
    definition = function(x){
        x@par[1L:x@nfact] <- 0
        x@est[1L:x@nfact] <- FALSE
        x
    }
)

# probability trace line function. Must return a matrix with a trace line for each category
setMethod(
    f = "ProbTrace",
    signature = signature(x = 'experimental', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- x@par[1L]
        b <- x@par[2L]
        p <- exp(a * (Theta - b)) / (1 + exp(a * (Theta - b)))
        p <- ifelse(p < 1e-10, 1e-10, p) #numerical constraints to avoid log() problems
        p <- ifelse(p > 1 - 1e-10, 1 - 1e-10, p)
        P <- cbind(1-p, p)
        return(P)
    }
)

# complete-data derivative used in parameter estimation (here it is done numerically)
setMethod(
    f = "Deriv",
    signature = signature(x = 'experimental', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(any(x@est)){
            grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta)
            if(estHess){
                hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x,
                                                            Theta=Theta)
            }
        }
        return(list(grad=grad, hess=hess)) #replace with analytical derivatives
    }
)

# derivative of the model wft to the Theta values (done numerically here)
setMethod(
    f = "DerivTheta",
    signature = signature(x = 'experimental', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_DerivTheta(x, Theta) #replace with analytical derivatives
    }
)

# derivative of the probability trace line function wrt Theta (done numerically here)
setMethod(
    f = "dP",
    signature = signature(x = 'experimental', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_dP(x, Theta) #replace with analytical derivatives
    }
)

# defines how the item should be initiallized with the S4 function new()
#   before the parameter estimates are updated
#   (set starting values, lower/upper bounds, logicals indicating estimation, etc)
setMethod("initialize",
          'experimental',
          function(.Object, nfact, ncat){
              stopifnot(nfact == 1L)
              stopifnot(ncat == 2L)
              .Object@par <- c(a=1, b=0)
              .Object@est <- c(TRUE, TRUE)
              .Object@lbound <- rep(-Inf, 2L)
              .Object@ubound <- rep(Inf, 2L)
              .Object
          }
)
