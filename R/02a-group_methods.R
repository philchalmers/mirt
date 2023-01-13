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
                        discrete='logical',
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
                                                             type = if(pars[[length(pars)]]@dentype == 'bfactor')
                                                                 "central" else "Richardson",
                                                             gradient=FALSE)
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

