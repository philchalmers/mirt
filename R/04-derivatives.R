#----------------------------------------------------------------------------
# Derivatives wrt item parameters

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
                    rs[,i]/P[,i]^2 * D2 * (PQfull[,i] - PQfull[,i+1L])^2
            hess[cind, cind] <- sum(ddc)
            for(i in 1L:nzetas)
                hess[cind, nfact + i] <- hess[nfact + i, cind] <-
                    sum((rs[,i]/P[,i] * D2 * (-Pfull[,i+1L] + 3*Pfull[,i+1L]^2 - 2*Pfull[,i+1L]^3) -
                    rs[,i]/P[,i]^2 * D2 * (PQfull[,i] - PQfull[,i+1L]) * (-PQfull[,i+1L]) +
                    rs[,i+1L]/P[,i+1L] * D2 * (Pfull[,i+1L] - 3*Pfull[,i+1L]^2 + 2*Pfull[,i+1L]^3) -
                    rs[,i+1L]/P[,i+1L]^2 * D2 * (PQfull[,i+1L] - PQfull[,i+2L]) * (PQfull[,i+1L])))
            for(j in 1L:nfact){
                tmp <- 0
                for(i in 1L:ncol(rs))
                        tmp <- tmp + (rs[,i]/P[,i] * D2 * Theta[,j] *
                                          (Pfull[,i] - 3*Pfull[,i]^2 + 2*Pfull[,i]^3 -
                                               Pfull[,i+1L] + 3*Pfull[,i+1L]^2 - 2*Pfull[,i+1L]^3) -
                                 rs[,i]/P[,i]^2 * D2 * (PQfull[,i] - PQfull[,i+1L]) * Theta[,j] *
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
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        #local derivative from previous version with small mod
        #u and g in logit form
        dpars.comp <- function(lambda,zeta,g,r,f,Thetas,D,estHess)
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
    f = "Deriv",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        #u and g in logit
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        dat <- x@dat
        if(estHess)
            hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta)
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
                    rowSums(idat[,-j]) * rowSums(Theta) * Pn[,j]))
            grad[nfact+3L+nd+j] <- sum((
                (idat[,j] * Qd * (Pn[,j] - Pn[,j]^2) * den) / (Qd * num[,j]) -
                    rowSums(idat[,-j]) * Pn[,j]))
        }
        ret <- list(grad=grad, hess=hess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        dat <- x@dat
        nfact <- x@nfact
        nzetas <- ncol(dat)
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        shift <- d[length(d)]
        dshift <- d <- d[-length(d)]
        dshift[-1L] <- d[-1L] + shift
        ak <- 0:(length(d)-1L)
        if(nrow(x@fixed.design) > 1L && ncol(x@fixed.design) > 0L)
            Theta <- cbind(x@fixed.design, Theta)
        P <- ProbTrace(x=x, Theta=Theta, useDesign = FALSE, ot=offterm)
        oldpar <- x@par
        tmp <- .Call("dparsNominal", x, Theta, offterm, TRUE, estHess)        
        num <- P.nominal(c(a, ak, dshift), ncat=length(ak), Theta=Theta, returnNum=TRUE, ot=offterm)
        grad <- tmp$grad
        hess <- tmp$hess

        #quick calcs for derivs
        nfact <- length(a)
        ncat <- length(d)
        akind <- nfact
        dind <- nfact + ncat*2 + 1L #go backwards
        ak2 <- ak^2
        P2 <- P^2
        P3 <- P^3
        aTheta <- as.vector(Theta %*% a)
        aTheta2 <- aTheta^2
        dat_num <- dat/num
        numsum <- rowSums(num)
        numD <- num %*% c(0, rep(1, ncol(num)-1L))
        numak <- matrix(num %*% ak, nrow(Theta), ncol(Theta))
        numakThetaD <- numak * Theta
        numD2 <- num %*% c(0, rep(1, ncol(num)-1L))
        numakThetaD2 <- numak * Theta
        ak0 <- ak
        ak0[1L] <- 0
        cind <- length(grad)
        if(estHess){
            tmp <- 0
            for(i in 1L:nzetas)
                tmp <- tmp + dat[,i]*numD^2 / numsum^2 - dat[,i]*numD2/numsum
            hess[cind, cind] <- sum(tmp)
            for(j in 1L:nzetas){
                tmp <- 0
                for(i in 1L:nzetas)
                    tmp <- tmp + dat[,i]*P[,j]*numD/numsum - dat[,i]*P[,j]
                hess[cind, nfact+j] <- hess[nfact+j, cind] <- sum(tmp)
            }
            for(j in 1L:nfact){
                tmp <- 0
                for(i in 1L:nzetas)
                    tmp <- tmp + dat[,i]*numD*numakThetaD[,j]/numsum^2 -
                        dat[,i]* (num %*% ak0*Theta[,j])/numsum
                hess[cind, j] <- hess[j, cind] <- sum(tmp)
            }
        }
        ####
        #FIXME - can't seem to get the last value of the gradient quite right for some reason....
        x2 <- x
        x2@est <- c(rep(FALSE, length(x2@est)-1L), TRUE)
        grad[x2@est] <- numDeriv::grad(EML, x@par[x2@est], obj=x2, Theta=Theta)
        ####
        ret <- list(grad=grad, hess=hess)
        if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess)
        ret
    }
)

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
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta, CUSTOM.IND, EM = FALSE, pars = NULL, itemloc = NULL,
                          tabdata = NULL, prior = NULL, estHess=FALSE){
        if(EM){
            grad <- rep(0, length(x@par))
            hess <- matrix(0, length(x@par), length(x@par))
            if(estHess){
                if(any(x@est)){
                    hess[x@est,x@est] <- numDeriv::hessian(EML2, x@par[x@est], Theta=Theta,
                                                           pars=pars, tabdata=tabdata,
                                                           itemloc=itemloc, CUSTOM.IND=CUSTOM.IND)
                }
                return(list(grad=grad, hess=hess))
            }
            J <- length(pars) - 1L
            nfact <- x@nfact
            scores <- matrix(0, nrow(tabdata), nfact)
            r <- tabdata[ ,ncol(tabdata)]
            N <- sum(r)
            tabdata <- tabdata[ ,-ncol(tabdata)]
            itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc, 
                                          CUSTOM.IND=CUSTOM.IND)
            mu <- x@par[1L:nfact]
            siglong <- x@par[-(1L:nfact)]
            sig <- matrix(0,nfact,nfact)
            selcov <- lower.tri(sig, diag=TRUE)
            scores2 <- matrix(0, nrow(tabdata), sum(selcov))
            thetas2 <- numeric(sum(selcov))
            ret <- .Call('EAPgroup', log(itemtrace), tabdata, Theta, prior, mu)
            tmp <- cbind(ret$scores, ret$scores2) * r
            newpars <- apply(tmp, 2, sum) / N
            if(nfact > 1L){
                x@par[x@est] <- newpars[x@est]
                cov <- ExtractGroupPars(x)$gcov
                ev <- eigen(cov)
                if(any(ev$values <= 0)){
                    eval <- ev$values
                    eval[eval < 0] <- 100*.Machine$double.eps
                    eval <- eval / sum(eval) * sum(ev$values)
                    cov <- ev$vectors %*% diag(eval) %*% t(ev$vectors)
                    newpars[(nfact+1L):length(newpars)] <- cov[lower.tri(cov, TRUE)]
                }
            }
            return(newpars[x@est])
        }
        tr <- function(y) sum(diag(y))
        nfact <- x@nfact
        N <- nrow(Theta)
        u <- x@par[1L:nfact]
        MU <- matrix(rep(u, N), N, byrow = TRUE)
        siglong <- x@par[-(1L:nfact)]
        sig <- matrix(0,nfact,nfact)
        selcov <- lower.tri(sig, diag=TRUE)
        sig[selcov] <- siglong
        if(nfact != 1L)
            sig <- sig + t(sig) - diag(diag(sig))
        npars <- length(sig) + nfact
        invSig <- solve(sig)
        Z <- t(Theta-MU) %*% (Theta-MU)
        g1 <- N * invSig %*% (colMeans(Theta) - u)
        tmp <- invSig %*% (Z - N * sig) %*% invSig
        diag(tmp) <- diag(tmp)/2 #correct for symmetry
        g2 <- tmp[selcov]
        grad <- c(g1, g2)
        sel <- 1L:npars
        cMeans <- N*(colMeans(Theta) - u)
        Zdif <- (Z - N * sig)
        hess <- .Call("dgroup",
                   as.numeric(invSig),
                   as.numeric(cMeans),
                   as.numeric(Zdif),
                   as.integer(N),
                   as.integer(nfact),
                   as.integer(npars))
        sel <- sel[c(rep(TRUE,nfact),as.logical(selcov))]
        hess <- hess[sel,sel]
        return(list(hess=hess,grad=grad))
    }
)

setMethod(
    f = "RandomDeriv",
    signature = signature(x = 'RandomPars'),
    definition = function(x){
        tr <- function(y) sum(diag(y))
        nfact <- x@ndim
        Theta <- x@drawvals
        N <- nrow(Theta)
        u <- rep(0, nfact)
        MU <- matrix(rep(u, N), N, byrow = TRUE)
        siglong <- x@par
        sig <- matrix(0,nfact,nfact)
        selcov <- lower.tri(sig, diag=TRUE)
        sig[selcov] <- siglong
        if(nfact != 1L)
            sig <- sig + t(sig) - diag(diag(sig))
        npars <- length(sig) + nfact
        invSig <- solve(sig)
        Z <- t(Theta-MU) %*% (Theta-MU)
        g1 <- N * invSig %*% (colMeans(Theta) - u)
        tmp <- invSig %*% (Z - N * sig) %*% invSig
        diag(tmp) <- diag(tmp)/2 #correct for symmetry
        g2 <- tmp[selcov]
        grad <- c(g1, g2)
        sel <- 1L:npars
        cMeans <- N*(colMeans(Theta) - u)
        Zdif <- (Z - N * sig)
        hess <- .Call("dgroup",
                      as.numeric(invSig),
                      as.numeric(cMeans),
                      as.numeric(Zdif),
                      as.integer(N),
                      as.integer(nfact),
                      as.integer(npars))
        sel <- sel[c(rep(TRUE,nfact),as.logical(selcov))]
        hess <- hess[sel,sel]
        hess <- hess[(nfact+1L):length(grad), (nfact+1L):length(grad), drop=FALSE]
        diag(hess) <- -abs(diag(hess))
        grad <- grad[(nfact+1L):length(grad)]
        return(list(hess=hess,grad=grad))
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
        else grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta)
        if(estHess){
            if(x@usehss) hess <- x@hss(x, Theta)
            else hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x,
                                                         Theta=Theta)
        }
        return(list(grad = grad, hess=hess))
    }
)

#----------------------------------------------------------------------------
EML <- function(par, obj, Theta){
    obj@par[obj@est] <- par
    itemtrace <- ProbTrace(x=obj, Theta=Theta)
    LL <- sum(obj@dat * log(itemtrace))
    LL <- LL.Priors(x=obj, LL=LL)
    return(LL)
}

EML2 <- function(x, Theta, pars, tabdata, itemloc, CUSTOM.IND){
    obj <- pars[[length(pars)]]
    obj@par[obj@est] <- x
    r <- tabdata[, ncol(tabdata)]
    gpars <- ExtractGroupPars(obj)
    mu <- gpars$gmeans
    sigma <- gpars$gcov
    prior <- mvtnorm::dmvnorm(Theta, mean=mu, sigma=sigma)
    prior <- prior/sum(prior)
    rlist <- Estep.mirt(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior, itemloc=itemloc,
                        CUSTOM.IND=CUSTOM.IND)
    LL <- sum(r*log(rlist$expected))
    LL <- LL.Priors(x=obj, LL=LL)
    return(LL)
}

difexp <- function(x) x * (1 - x)

dif2exp <- function(x) 2 * (x * (1 - x)^2)

#----------------------------------------------------------------------------
# Derivatives wrt Theta, returns list with number of categories, and
#    inside a matrix with the number of factors

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
        P <- P.comp(c(a, d, g, 999), Theta)[,2L]
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

###

setMethod(
    f = "dP",
    signature = signature(x = 'dich', Theta = 'matrix', prior = 'numeric'),
    definition = function(x, Theta, prior){
        P <- ProbTrace(x, Theta)
        tmp <- P * (1-P)
        tmp[,1] <- tmp[,1] * Theta
        return(colSums(tmp*prior)) #FIXME U2PL only
    }
)