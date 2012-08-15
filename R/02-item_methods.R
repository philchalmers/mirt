#Probability Traces
setMethod(
    f = "ProbTrace",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta){  
        u <- x@par[length(x@par)]
        g <- x@par[length(x@par)-1]
        d <- x@par[length(x@par)-2]
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        P <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u)
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){                  
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        d <- x@par[(x@nfact+1):length(x@par)]
        P <- P.poly(a=a, d=d, Theta=Theta, itemexp=TRUE)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){                  
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        d <- x@par[-(1:x@nfact)]
        P <- P.gpcm(a=a, d=d, Theta=Theta)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){         
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        ak <- x@par[(x@nfact+1):(x@nfact + x@ncat)]
        d <- x@par[length(x@par):(length(x@par) - x@ncat + 1)]
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){    
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(length(x@par)-2)]
        g <- x@par[length(x@par)-1]
        u <- x@par[length(x@par)]
        if(x@bfactor) a <- a[x@est[1:nfact]]
        if(x@bfactor) d <- d[x@est[(nfact+1):(nfact*2)]]        
        P <- P.comp(a=a, d=d, Theta=Theta, g=g, u=u)
        return(P)
    }
)

#----------------------------------------------------------------------------
#LogLik
setMethod(
    f = "LogLik",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        LL <- (-1) * sum(x@rs * log(itemtrace))
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        LL <- (-1) * sum(x@rs * log(itemtrace))
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        LL <- (-1) * sum(x@rs * log(itemtrace))
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        LL <- (-1) * sum(x@rs * log(itemtrace))
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        LL <- (-1) * sum(x@rs * log(itemtrace))
        return(LL)
    }
)

#----------------------------------------------------------------------------
setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'dich'),
    definition = function(x){          
        par <- x@par
        a <- par[1:x@nfact]
        a        
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'graded'),
    definition = function(x){          
        par <- x@par
        a <- par[1:x@nfact]
        a        
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'gpcm'),
    definition = function(x){          
        par <- x@par
        a <- par[1:x@nfact]
        a        
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'nominal'),
    definition = function(x){          
        par <- x@par
        a <- par[1:x@nfact]
        a        
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'partcomp'),
    definition = function(x){          
        par <- x@par
        a <- par[1:x@nfact]
        a        
    }
)

#----------------------------------------------------------------------------
setMethod(
    f = "Deriv",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta){          
        nfact <- x@nfact
        parlength <- length(x@par)
        u <- x@par[parlength]
        g <- x@par[parlength - 1]
        d <- x@par[parlength - 2]
        a <- x@par[1:nfact]        
        P <- P.mirt(a, d, Theta, g, u)
        Q <- 1 - P
        dat <- x@dat[,2] 
        d2L <- matrix(0,nfact+3, nfact+3)						
        if(!x@est[parlength] && !x@est[parlength-1]){ #'2PL'
            PQ <- P*Q
            L1 <- colSums((dat-P) * Theta)
            L2 <- sum(dat-P)
            dL <- c(L1,L2,0,0)    	
            L11 <- .Call("dichOuter", Theta, PQ, nrow(Theta))
            d2L[1:nfact, 1:nfact] <- -L11
            d2L[nfact+1, nfact+1] <- (-1)*sum(PQ)		
            d2L[nfact+1, 1:nfact] <- d2L[1:nfact, nfact+1] <- (-1)*colSums(PQ * Theta)
        } else if(!x@est[parlength] && x@est[parlength-1]){ #'3PL'
            f <- 1		
            Pstar <- P.mirt(a,d,Theta,0,1)		
            Qstar <- 1 - Pstar
            da <- rep(0,nfact)	
            dd <- sum((1-g)*Pstar*Qstar*(dat/P - (f-dat)/Q))
            dc <- sum(Qstar*(dat/P - (f-dat)/Q))
            for(i in 1:nfact){
                da[i] <- sum(Theta[,i]*Pstar*Qstar*(1-g)*(dat/P - (f-dat)/Q))
            }
            dL <- c(da,dd,dc,0)				
            gloc <- nfact+2
            const1 <- (dat/P - (f-dat)/Q)*(Qstar-Pstar)
            const2 <- (dat/P^2 + (f-dat)/Q^2)	
            d2L[nfact+1,nfact+1] <- sum((1-g)*Pstar*Qstar*(const1 - 
                Pstar*Qstar*(1-g)*const2))		
            d2L[gloc,gloc] <- -sum(Qstar^2 *(dat/P^2 + (f-dat)/Q^2))
            d2L[gloc,nfact+1] <- d2L[nfact+1,gloc] <- sum(-Pstar*Qstar*((dat/P - (f-dat)/Q) + Qstar*(1-g)*const2)) 
            for(i in 1:nfact){
                d2L[nfact+1,i] <- d2L[i,nfact+1] <- sum((1-g)*Theta[,i]*Pstar*Qstar*(const1 - 
                    Pstar*Qstar*(1-g)*const2))			
                d2L[gloc,i] <- d2L[i,gloc] <- sum(-Theta[,i]*Pstar*Qstar*((dat/P - (f-dat)/Q) + 
                    Qstar*(1-g)*const2))		
                for(j in 1:nfact){
                    if(i == j)
                        d2L[i,i] <- sum(Theta[,i]^2 *Pstar*Qstar*(1-g)*(const1 - 
                            (1-g)*Pstar*Qstar*const2))
                    if(i < j)
                        d2L[i,j] <- d2L[j,i] <- sum(Theta[,i]*Theta[,j] *Pstar*Qstar*(1-g)*
                            (const1 - (1-g)*Pstar*Qstar*const2))					
                }
            }	
        }  	
        return(list(grad = dL, hess = d2L))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){         
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[-(1:nfact)]
        nd <- length(d)    			
        P <- P.poly(a, d,Theta)			    	
        ret <- .Call("dparsPoly", P, Theta, x@dat, nd)
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(nfact*2)]
        g <- x@par(length(x@par)-1)
        u <- x@par(length(x@par))
        if(x@est[nfact*2 + 1] && !x@est[nfact*2+2]){
            grad <- function(a, d, g, u, r, Theta){
                f <- 1			
                P <- P.comp(a,d,Theta,g,1)		
                Pstar <- P.comp(a,d,Theta,0)		
                Qstar <- 1 - Pstar
                Q <- 1 - P
                const1 <- (r/P - (f-r)/Q)
                dd <- da <- rep(0,nfact)		
                dg <- sum(Qstar*const1)
                for(i in 1:nfact){
                    Pk <- P.mirt(a[i],d[i],Theta[ , i, drop=FALSE],0)
                    Qk <- 1 - Pk
                    dd[i] <- sum((1-g)*Pstar*Qk*const1)
                    da[i] <- sum((1-g)*Pstar*Qk*Theta[,i]*const1)
                }
                return(c(da,dd,dg,0))
            }		
            hess <- function(a, d, g, u, r, Theta){ 
                nfact <- length(a)
                d2L <- matrix(0, nfact*2 + 2, nfact*2 + 2)
                f <- 1			
                P <- P.comp(a,d,Theta,g, 1)		
                Pstar <- P.comp(a,d,Theta,0, 1)		
                Qstar <- 1 - Pstar
                Q <- 1 - P	
                const1 <- (r/P - (f-r)/Q)
                const2 <- (r/P^2 + (f-r)/Q^2)	
                Names <- c(paste("a",1:nfact,sep='_'),paste("d",1:nfact,sep='_'),'g_0')
                for(i in 1:(nfact*2+1)){		
                    for(j in 1:(nfact*2+1)){
                        if(i <= j){
                            d1 <- strsplit(Names[c(i,j)],"_")[[1]]
                            d2 <- strsplit(Names[c(i,j)],"_")[[2]]
                            k <- as.numeric(d1[2])
                            m <- as.numeric(d2[2])
                            Pk <- P.mirt(a[k],d[k],Theta[ , k, drop=FALSE],0)
                            Qk <- 1 - Pk	
                            Pm <- P.mirt(a[m],d[m],Theta[ , m, drop=FALSE],0)
                            Qm <- 1 - Pm									
                            if(i == j && d1[1] == 'd'){
                                d2L[i,i] <- sum((1-g)*Pstar*Qk*(const1*((1-g)*Qk - Pk) - Pstar*Qk*(1-g)*const2))
                                next
                            }
                            if(i == j && d1[1] == 'a'){
                                d2L[i,i] <- sum((1-g)*Theta[,k]^2*Pstar*Qk*(const1*((1-g)*Qk - Pk) - Pstar*Qk*
                                    (1-g)*const2))
                                next		
                            }
                            if(i == j && d1[1] == 'g'){
                                d2L[i,i] <- -sum(Qstar^2 * const2)
                                next		
                            }	
                            if(d1[1] == 'a' && d2[1] == 'a'){
                                d2L[i,j] <- d2L[j,i] <- sum((1-g)*Theta[,k]*Theta[,m]*Qk*Pstar*Qm*(const1 - 
                                    Pstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'd'){
                                d2L[i,j] <- d2L[j,i] <- sum((1-g)*Qk*Pstar*Qm*(const1 - Pstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'a' && d2[1] == 'g'){
                                d2L[i,j] <- d2L[j,i] <- -sum(Theta[,k]*Pstar*Qk*(const1 + Qstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'g'){
                                d2L[i,j] <- d2L[j,i] <- -sum(Pstar*Qk*(const1 + Qstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                                d2L[i,j] <- d2L[j,i] <- sum((1-g)*Theta[,k]*Pstar*Qk*(const1*((1-g)*Qk - Pk) - 
                                    Pstar*Qk*(1-g)*const2))
                                next	
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                                d2L[i,j] <- d2L[j,i] <- sum((1-g)*Qk*Theta[,m]*Pstar*Qm*(const1 - 
                                    Pstar*(1-g)*const2))
                                next
                            }						
                        }
                    }
                }	
                return(d2L)
            }		
            return(list(grad = grad(a, d, g, u, x@dat, Theta), hess = hess(a, d, g, u, x@dat, Theta)))
        }
        if(!x@est[nfact*2 + 1] && !x@est[nfact*2+2]){
            d2L <- matrix(0, nfact*2 + 2, nfact*2 + 2)
            P <- P.comp(a,d,Theta)	
            Q <- 1 - P	
            da <- dd <- rep(0,nfact)	
            for(i in 1:nfact){
                Pk <- P.mirt(a[i],d[i],Theta[ , i, drop=FALSE],0, 1)
                Qk <- 1 - Pk
                const <- (1 - dat)*P/Q
                dd[i] <- sum(Qk*(dat - const))
                da[i] <- sum(Theta[,i]*Qk*(dat - const))
            }
            Names <- c(paste("a",1:nfact,sep='_'),paste("d",1:nfact,sep='_'))
            f <- 1
            r <- dat
            for(i in 1:(nfact*2)){		
                for(j in 1:(nfact*2)){
                    if(i <= j){
                        d1 <- strsplit(Names[c(i,j)],"_")[[1]]
                        d2 <- strsplit(Names[c(i,j)],"_")[[2]]
                        k <- as.numeric(d1[2])
                        m <- as.numeric(d2[2])
                        Pk <- P.mirt(a[k],d[k],Theta[ , k, drop=FALSE],0)
                        Qk <- 1 - Pk	
                        Pm <- P.mirt(a[m],d[m],Theta[ , m, drop=FALSE],0)
                        Qm <- 1 - Pm									
                        if(i == j && d1[1] == 'd'){
                            d2L[k,k] <- sum(-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2)
                            next
                        }
                        if(i == j && d1[1] == 'a'){
                            d2L[k+nfact,k+nfact] <- sum(Theta[,k]^2 *
                                (-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2))
                            next		
                        }				
                        if(d1[1] == 'a' && d2[1] == 'a'){
                            d2L[i,j] <- d2L[j,i] <- -sum(Theta[,k]*Theta[,m]*Qk*Qm*(f-r)*P/Q^2) 
                            next
                        }
                        if(d1[1] == 'd' && d2[1] == 'd'){
                            d2L[i,j] <- d2L[j,i] <- -sum(Qk*Qm*(f-r)*P/Q^2)
                            next
                        }	
                        if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                            d2L[i,j] <- d2L[j,i] <- sum(Theta[,k]*Qk*(-Pk*(r - (f-r)*P/Q) - 
                                Qk*(f-r)*P/Q^2))
                            next	
                        }
                        if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                            d2L[i,j] <- d2L[j,i] <- -sum(Qk*Qm*Theta[,m]*(f-r)*P/Q^2)
                            next
                        }						
                    }
                }
            }		
            return(list(grad = c(da,dd,0,0), hess = d2L)) 
        }	
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta){
        tr <- function(y) sum(diag(y))
        nfact <- x@nfact
        N <- nrow(Theta)
        u <- x@par[1:nfact]
        siglong <- x@par[-(1:nfact)]
        sig <- matrix(0,nfact,nfact)
        selcov <- lower.tri(sig, diag=TRUE) 
        sig[selcov] <- siglong
        sig <- sig + t(sig) - diag(diag(sig))
        npars <- length(sig) + nfact	
        g <- rep(0,nfact + nfact*(nfact+1)/2)	
        invSig <- solve(sig)	
        Z <- t(Theta-u) %*% (Theta-u)
        g[1:nfact] <- N * invSig %*% (colMeans(Theta) - u) 		
        tmp <- .5 * invSig %*% (Z - N * sig) %*% invSig  
        g[(nfact+1):length(g)] <- tmp[selcov]
        h <- matrix(0,npars,npars)
        sel <- 1:npars		
        cMeans <- N*(colMeans(Theta) - u)
        Zdif <- (Z - N * sig)		
        h <- .Call("dgroup",				
                   as.numeric(invSig),
                   as.numeric(cMeans),				
                   as.numeric(Zdif),
                   as.integer(N),
                   as.integer(nfact),
                   as.integer(npars))				
        sel <- sel[c(rep(TRUE,nfact),as.logical(selcov))]	
        h <- h[sel,sel] 
        return(list(hess=h,grad=g))
    }
)

#----------------------------------------------------------------------------
##Function passes
#trace lines for polymirt
P.poly <- function(a, d, Theta, itemexp = FALSE)
{    
    ncat <- length(d) + 1
    nfact <- length(a)
    Pk <- matrix(0,nrow(Theta),ncat+1)
    Pk[,1] <- 1	
    for(i in 1:(ncat-1))			
        Pk[ ,i+1] <- P.mirt(a, d[i], Theta, 0, 1)		
    if(itemexp){
        P <- matrix(0,nrow(Theta),ncat)		
        for(i in ncat:1)
            P[ ,i] <- Pk[ ,i] - Pk[ ,i+1]						
        Pk <- P
    }	
    return(Pk)
}

# Trace lines for mirt models
P.mirt <- function(a, d, Theta, g, u)
{ 		
    traces <- .Call("traceLinePts", a, d, g, u, Theta)
    return(traces)
}

# Trace lines for partially compensetory models
P.comp <- function(a, d, Theta, g = 0, u = 1)
{
    nfact <- length(a)
    P <- rep(1,nrow(Theta))
    for(i in 1:nfact)
        P <- P * P.mirt(a[i], d[i], Theta[ ,i, drop=FALSE], 0, 1)
    P <- g + (u - g) * P
    P	
}

# Trace lines for bfactor
P.bfactor <- function(a, d, Theta, g, u, patload)
{ 
    a <- a[patload]	
    if(length(d) > 1){
        ncat <- length(d) + 1
        nfact <- length(a)
        Pk <- matrix(0,nrow(Theta),ncat+1)
        Pk[,1] <- 1	
        for(i in 1:(ncat-1))			
            Pk[ ,i+1] <- P.mirt(a, d[i], Theta, 0)				
        P <- matrix(0,nrow(Theta),ncat)		
        for(i in ncat:1)
            P[ ,i] <- Pk[ ,i] - Pk[ ,i+1]		
    } else P <- .Call("traceLinePts", a, d, g, u, Theta)		
    return(P)
}

#d[1] == 0, ak[1] == 0, ak[length(ak)] == length(ak) - 1 
P.nominal <- function(a, ak, d, Theta){
    ncat <- length(d)
    nfact <- ncol(Theta)    
    a <- matrix(a)    
    numerator <- matrix(0, nrow(Theta), ncat)            
    for(i in 1:ncat)
        numerator[ ,i] <- exp(1.702 * ak[i] * (Theta %*% a) + 1.702 * d[i])
    P <- numerator/rowSums(numerator)
    return(P)   
}

#d[1] == 0
P.gpcm <- function(a, d, Theta){ 
    ncat <- length(d)
    nfact <- ncol(Theta)            
    k <- 0:(ncat-1)
    numerator <- matrix(0, nrow(Theta), ncat)    
    a <- matrix(a)    
    for(i in 1:ncat)
        numerator[ ,i] <- exp(1.702 * k[i] * (Theta %*% a) + 1.702 * d[i])
    P <- numerator/rowSums(numerator)
    return(P)   
}

