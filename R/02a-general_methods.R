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
                               Theta=Theta, prior=prior, wmiss=rep(1, nrow(tabdata)),
                               Priorbetween=Priorbetween, specific=bfactor_info$specific,
                               sitems=sitems, itemloc=itemloc, CUSTOM.IND=CUSTOM.IND, omp_threads=1L)
    } else {
        rlist <- Estep.mirt(pars=pars, tabdata=tabdata, freq=freq, wmiss=rep(1, nrow(tabdata)),
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
    r_P[is.nan(r_P)] <- 0
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

Deriv.mix <- function(x, estHess=FALSE){
    phi2psi <- function(par){
        E <- exp(par)
        E / sum(E)
    }
    LL <- function(par, x){
        phi <- x$par
        phi[x$est] <- par
        psi <- phi2psi(phi)
        dmultinom(x$dat, prob=psi, log=TRUE)
    }
    ret <- list(grad=numeric(length(x$par)),
                hess=matrix(0, length(x$par), length(x$par)))
    ret$grad[x$est] <- numerical_deriv(par=x$par[x$est], f=LL, x=x)
    if(estHess && any(x$est))
        ret$hess[x$est, x$est] <- numerical_deriv(par=x$par[x$est],
                                                  f=LL, x=x, gradient=FALSE)
    ret
}

# ----------------------------------------------------------------
# global S3 methods


#' Print generic for customized data.frame console output
#'
#' Provides a nicer output for most printed \code{data.frame} objects
#' defined by functions in \code{mirt}.
#'
#' @method print mirt_df
#' @param x object of class \code{'mirt_df'}
#' @param digits number of digits to round
#' @param ... additional arguments passed to \code{print(...)}
#' @export
print.mirt_df <- function(x, digits = 3, ...){
    cls <- class(x)[2L]
    class(x) <- cls
    if(nrow(x) > 0){
        clsss <- sapply(x, class)
        for(i in 1:length(clsss))
            if(clsss[i] == 'numeric')
                x[,i] <- round(x[,i], digits=digits)
    }
    if(!is.null(x[['p']])){
        if(!is.null(x$X2)) x$X2 <- as.character(x$X2)
        if(!is.null(x$df)) x$df <- as.character(x$df)
        if(!is.null(x$p)) x$p <- as.character(x$p)
        print(x, na.print = " ", ...)
    } else {
        print(x, ...)
    }
}

#' Print generic for customized matrix console output
#'
#' Provides a nicer output for most printed \code{matrix}
#' objects defined by functions in \code{mirt}.
#'
#' @method print mirt_matrix
#' @param x object of class \code{'mirt_matrix'}
#' @param digits number of digits to round
#' @param ... additional arguments passed to \code{print(...)}
#' @export
print.mirt_matrix <- function(x, digits = 3, ...){
    cls <- class(x)[2L]
    class(x) <- cls
    x <- round(x, digits=digits)
    print(x, ...)
}

#' @method print matrix
#' @rdname print.mirt_matrix
#' @export
print.matrix <- function(x, ...){
    if(is(x, 'mirt_matrix')){
        class(x) <- 'mirt_matrix'
        print(x, ...)
    } else {
        print.default(x, ...)
    }
    invisible(NULL)
}

#' Head generic for customized matrix console output
#'
#' Provides a nicer output for most printed \code{matrix}
#' objects defined by functions in \code{mirt}.
#'
#' @method head mirt_matrix
#' @param x object of class \code{'mirt_matrix'}
#' @param digits number of digits to round
#' @param ... additional arguments passed to \code{print(...)}
#' @export
head.mirt_matrix <- function(x, digits = 3, ...){
    cls <- class(x)[2L]
    class(x) <- cls
    x <- round(x, digits=digits)
    head(x, ...)
}

#' Tail generic for customized matrix console output
#'
#' Provides a nicer output for most printed \code{matrix}
#' objects defined by functions in \code{mirt}.
#'
#' @method tail mirt_matrix
#' @param x object of class \code{'mirt_matrix'}
#' @param digits number of digits to round
#' @param ... additional arguments passed to \code{print(...)}
#' @export
tail.mirt_matrix <- function(x, digits = 3, ...){
    cls <- class(x)[2L]
    class(x) <- cls
    x <- round(x, digits=digits)
    tail(x, ...)
}

#' Print generic for customized list console output
#'
#' Provides a nicer output for most printed \code{list} objects
#' defined by functions in \code{mirt}.
#'
#' @method print mirt_list
#' @param x object of class \code{'mirt_list'}
#' @param digits number of digits to round
#' @param ... additional arguments passed to \code{print(...)}
#' @export
print.mirt_list <- function(x, digits = 3, ...){
    cls <- class(x)[2L]
    class(x) <- cls
    x <- lapply(x, function(x, digits){
        if(is.list(x) && !is.data.frame(x))
            lapply(x, function(y) round(y, digits=digits))
        else round(x, digits=digits)
    }, digits=digits)
    print(x, ...)
}

# ----------------------------------------------------------------
# Begin class and method definitions

setClass("GroupPars",
         representation(par='numeric',
                        SEpar='numeric',
                        parnames='character',
                        est='logical',
                        parnum='numeric',
                        itemclass='integer',
                        dat='matrix',
                        nfact='integer',
                        standardize='logical',
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
                        structure='matrix',
                        theta='matrix',
                        Thetabetween='matrix',
                        density='numeric',
                        dentype='character',
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
        if(x@dentype == "Davidian"){
            pick <- x@est
            pick[1L:2L] <- FALSE
            x@par[pick] <- runif(sum(pick), -pi/2, pi/2)
        }
        x
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta, CUSTOM.IND, EM = FALSE, pars = NULL, itemloc = NULL,
                          tabdata = NULL, freq = NULL,
                          estHess=FALSE, prior = NULL, bfactor_info=NULL){
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
                else grad[x@est] <- numerical_deriv(x@par[x@est], LLfun, obj=x,
                                                    Theta=Theta, type=x@derivType)
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
                        hess[x@est,x@est] <- numerical_deriv(x@par[x@est], EML2, Theta=Theta,
                                                             pars=pars, tabdata=tabdata, freq=freq,
                                                             itemloc=itemloc, CUSTOM.IND=CUSTOM.IND,
                                                             bfactor_info=bfactor_info,
                                                             type = if(ncol(Theta) > 2)  "central"
                                                                    else 'Richardson', gradient=FALSE)
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
                        parnames='character',
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
        sigma <- if(ncol(theta0) == 1L) matrix(x@cand.t.var)
           else diag(rep(x@cand.t.var,ncol(theta0)))
        theta1 <- theta0 + mirt_rmvnorm(N, prior.mu, sigma)
        if(is.null(total_0)) theta1 <- theta0 #for initial draw
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
            for(i in seq_len(J))
                LL2[,i] <- rowSums(LL[,itemloc[i]:(itemloc[i+1L] - 1L)])
        }
        total_1 <- tapply(LL2, x@mtch, sum) + log_den1
        if(is.null(total_0)){ #for intial draw
            attr(theta1, 'log.lik_full') <- total_1
            return(theta1)
        }
        diff <- total_1 - total_0
        accept <- log(unif) < diff
        theta1[!accept, ] <- theta0[!accept, ]
        total_1[!accept] <- total_0[!accept]
        attr(theta1, "Proportion Accepted") <- sum(accept)/N
        attr(theta1, 'log.lik_full') <- total_1
        return(theta1)
    }
)

RandomDeriv <- function(x, estHess = TRUE){
    Theta <- x@drawvals
    pick <- -seq_len(ncol(Theta))
    out <- .Call("dgroup", x, Theta, matrix(0L), estHess, TRUE, FALSE, FALSE)
    grad <- out$grad[pick]
    hess <- out$hess[pick, pick, drop=FALSE]
    diag(hess) <- -abs(diag(hess)) #hack for very small clusters
    list(grad=grad, hess=hess)
}

# ----------------------------------------------------------------

setClass("lrPars",
         representation(par='numeric',
                        SEpar='numeric',
                        parnames='character',
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
                        qr_XX='list',
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
    definition = function(x, cov, theta, estHess = TRUE){
        inv_sigma <- solve(cov)
        tmp <- t(inv_sigma %*% t(theta - x@mus) %*% x@X)
        tmp2 <- -det(inv_sigma) * x@tXX
        list(grad=tmp, hess=tmp2)
    }
)

