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
    signature = signature(x = 'grad', Theta = 'matrix'),
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
    signature = signature(x = 'nom', Theta = 'matrix'),
    definition = function(x, Theta){         
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        ak <- x@par[(x@nfact+1):(x@nfact + x@ncat)]
        d <- x@par[length(x@par):(length(x@par) - x@ncat + 1)]
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta)
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
    signature = signature(x = 'grad', Theta = 'matrix'),
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
    signature = signature(x = 'nom', Theta = 'matrix'),
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
    signature = signature(x = 'grad'),
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
    signature = signature(x = 'nom'),
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
        dat <- x@dat[,1, drop=FALSE] #FIXME change order now that lamabdas come first
        if(!x@est[parlength] && !x@est[parlength-1]){ '2PL'
            PQ <- P*(1-P)
            L1 <- sum(dat-P)
            L2 <- colSums((dat-P) * Theta)
            dL <- c(L1,L2)    	
            d2L <- matrix(0,nfact+1, nfact+1)						
            L11 <- .Call("dichOuter", Theta, PQ, nrow(Theta))
            if(nfact > 1) d2L[1:nfact+1, 1:nfact+1] <- -L11
                else d2L[nfact+1, nfact+1] <- -L11 				
            d2L[1, 1] <- (-1)*sum(PQ)		
            d2L[1, 1:nfact+1] <- d2L[1:nfact+1, 1] <- (-1)*colSums(PQ * Theta)
        } else if(!x@est[parlength] && x@est[parlength-1]){ '3PL'
            r <- dat
            f <- 1		
            c <- g            			
            thetas <- Theta
            Pstar <- P.mirt(lambda,zeta,Theta,0)		
            Qstar <- 1 - Pstar
            Q <- 1 - P
            da <- rep(0,nfact)	
            dd <- sum((1-g)*Pstar*Qstar*(r/P - (f-r)/Q))
            dc <- sum(Qstar*(r/P - (f-r)/Q))
            for(i in 1:nfact){
                da[i] <- sum(Theta[,i]*Pstar*Qstar*(1-g)*(r/P - (f-r)/Q))
            }
            dL <- c(dd,da,dc)				
            hess <- matrix(0,nfact + 2,nfact + 2)	
            aNames <- paste("a",1:nfact,sep='_')
            Names <- c('d',paste("a",1:nfact,sep='_'),'c')
            colnames(hess) <- rownames(hess) <- Names
            hsize <- nfact+2
            const1 <- (r/P - (f-r)/Q)*(Qstar-Pstar)
            const2 <- (r/P^2 + (f-r)/Q^2)	
            hess[1,1] <- sum((1-c)*Pstar*Qstar*(const1 - 
                Pstar*Qstar*(1-c)*const2))		
            hess[hsize,hsize] <- -sum(Qstar^2 *(r/P^2 + (f-r)/Q^2))
            hess[hsize,1] <- hess[1,hsize] <- sum(-Pstar*Qstar*((r/P - (f-r)/Q) + Qstar*(1-c)*const2)) 
            for(i in 1:nfact){
                hess[1,1+i] <- hess[1+i,1] <- sum((1-c)*thetas[,i]*Pstar*Qstar*(const1 - 
                    Pstar*Qstar*(1-c)*const2))			
                hess[hsize,1+i] <- hess[1+i,hsize] <- sum(-thetas[,i]*Pstar*Qstar*((r/P - (f-r)/Q) + 
                    Qstar*(1-c)*const2))		
                for(j in 1:nfact){
                    if(i == j)
                        hess[1+i,1+i] <- sum(thetas[,i]^2 *Pstar*Qstar*(1-c)*(const1 - 
                            (1-c)*Pstar*Qstar*const2))
                    if(i < j)
                        hess[1+i,1+j] <- hess[1+j,1+i] <- sum(thetas[,i]*thetas[,j] *Pstar*Qstar*(1-c)*
                            (const1 - (1-c)*Pstar*Qstar*const2))					
                }
            }	
            d2L <- hess			
        }  	
        return(list(grad = dL, hess = d2L))
    }
)



###### HERE #########


setMethod(
    f = "Deriv",
    signature = signature(x = 'grad', Theta = 'matrix'),
    definition = function(x, Theta){
        nzeta <- length(zeta)    			
        P <- P.poly(lambda,zeta,Thetas)			    	
        ret <- .Call("dparsPoly", P, Thetas, dat, nzeta)	
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta){
        tr <- function(x) sum(diag(x))
        x <- theta
        u <- grouplist$u    
        sig <- grouplist$sig
        N <- nrow(x)
        nfact <- length(u)
        selcov <- matrix(FALSE,nfact,nfact)
        selcov <- lower.tri(selcov) 
        diag(selcov) <- TRUE
        npars <- length(sig) + nfact	
        g <- rep(0,nfact + nfact*(nfact+1)/2)	
        invSig <- solve(sig)	
        Z <- t(x-u) %*% (x-u)
        g[1:nfact] <- N * invSig %*% (colMeans(x) - u) 		
        tmp <- .5 * invSig %*% (Z - N * sig) %*% invSig  
        g[(nfact+1):length(g)] <- tmp[selcov]
        h <- matrix(0,npars,npars)
        sel <- 1:npars		
        cMeans <- N*(colMeans(x) - u)
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
        list(h=h,g=g)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){
        nfact <- length(lambda)    
        pars <- c(zeta,lambda,g)
        if(estg){
            grad <- function(pars, r, thetas){
                f <- 1			
                d <- pars[1:nfact]	
                a <- pars[(nfact+1):(length(pars)-1)]
                c <- pars[length(pars)]
                P <- P.comp(a,d,thetas,c)		
                Pstar <- P.comp(a,d,thetas,0)		
                Qstar <- 1 - Pstar
                Q <- 1 - P
                const1 <- (r/P - (f-r)/Q)
                dd <- da <- rep(0,nfact)		
                dc <- sum(Qstar*const1)
                for(i in 1:nfact){
                    Pk <- P.mirt(a[i],d[i],thetas[ , i, drop=FALSE],0)
                    Qk <- 1 - Pk
                    dd[i] <- sum((1-c)*Pstar*Qk*const1)
                    da[i] <- sum((1-c)*Pstar*Qk*thetas[,i]*const1)
                }
                return(c(dd,da,dc))
            }		
            hess <- function(pars, r, thetas){
                f <- 1			
                d <- pars[1:nfact]	
                a <- pars[(nfact+1):(length(pars)-1)]
                c <- pars[length(pars)]
                P <- P.comp(a,d,thetas,c)		
                Pstar <- P.comp(a,d,thetas,0)		
                Qstar <- 1 - Pstar
                Q <- 1 - P	
                const1 <- (r/P - (f-r)/Q)
                const2 <- (r/P^2 + (f-r)/Q^2)	
                hess <- matrix(0,nfact*2+1,nfact*2+1)
                dNames <- paste("d",1:nfact,sep='_')
                aNames <- paste("a",1:nfact,sep='_')
                Names <- c(paste("d",1:nfact,sep='_'),paste("a",1:nfact,sep='_'),'c_0')
                for(i in 1:(nfact*2+1)){		
                    for(j in 1:(nfact*2+1)){
                        if(i <= j){
                            d1 <- strsplit(Names[c(i,j)],"_")[[1]]
                            d2 <- strsplit(Names[c(i,j)],"_")[[2]]
                            k <- as.numeric(d1[2])
                            m <- as.numeric(d2[2])
                            Pk <- P.mirt(a[k],d[k],thetas[ , k, drop=FALSE],0)
                            Qk <- 1 - Pk	
                            Pm <- P.mirt(a[m],d[m],thetas[ , m, drop=FALSE],0)
                            Qm <- 1 - Pm									
                            if(i == j && d1[1] == 'd'){
                                hess[i,i] <- sum((1-c)*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*(1-c)*const2))
                                next
                            }
                            if(i == j && d1[1] == 'a'){
                                hess[i,i] <- sum((1-c)*thetas[,k]^2*Pstar*Qk*(const1*((1-c)*Qk - Pk) - Pstar*Qk*
                                    (1-c)*const2))
                                next		
                            }
                            if(i == j && d1[1] == 'c'){
                                hess[i,i] <- -sum(Qstar^2 * const2)
                                next		
                            }	
                            if(d1[1] == 'a' && d2[1] == 'a'){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*thetas[,m]*Qk*Pstar*Qm*(const1 - 
                                    Pstar*(1-c)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'd'){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2))
                                next
                            }
                            if(d1[1] == 'a' && d2[1] == 'c'){
                                hess[i,j] <- hess[j,i] <- -sum(thetas[,k]*Pstar*Qk*(const1 + Qstar*(1-c)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'c'){
                                hess[i,j] <- hess[j,i] <- -sum(Pstar*Qk*(const1 + Qstar*(1-c)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*Pstar*Qk*(const1*((1-c)*Qk - Pk) - 
                                    Pstar*Qk*(1-c)*const2))
                                next	
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*thetas[,m]*Pstar*Qm*(const1 - 
                                    Pstar*(1-c)*const2))
                                next
                            }						
                        }
                    }
                }	
                return(hess)
            }		
            return(list(grad = grad(pars, dat, Thetas), hess = hess(pars, dat, Thetas)))
        } else {			
            P <- P.comp(lambda,zeta,Thetas)	
            Q <- 1 - P	
            da <- dd <- rep(0,nfact)	
            for(i in 1:nfact){
                Pk <- P.mirt(lambda[i],zeta[i],Thetas[ , i, drop=FALSE],0)
                Qk <- 1 - Pk
                const <- (1 - dat)*P/Q
                dd[i] <- sum(Qk*(dat - const))
                da[i] <- sum(Thetas[,i]*Qk*(dat - const))
            }
            hess <- matrix(0,nfact*2,nfact*2)
            dNames <- paste("d",1:nfact,sep='_')
            aNames <- paste("a",1:nfact,sep='_')
            Names <- c(paste("d",1:nfact,sep='_'),paste("a",1:nfact,sep='_'))
            f <- 1
            r <- dat
            for(i in 1:(nfact*2)){		
                for(j in 1:(nfact*2)){
                    if(i <= j){
                        d1 <- strsplit(Names[c(i,j)],"_")[[1]]
                        d2 <- strsplit(Names[c(i,j)],"_")[[2]]
                        k <- as.numeric(d1[2])
                        m <- as.numeric(d2[2])
                        Pk <- P.mirt(lambda[k],zeta[k],Thetas[ , k, drop=FALSE],0)
                        Qk <- 1 - Pk	
                        Pm <- P.mirt(lambda[m],zeta[m],Thetas[ , m, drop=FALSE],0)
                        Qm <- 1 - Pm									
                        if(i == j && d1[1] == 'd'){
                            hess[k,k] <- sum(-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2)
                            next
                        }
                        if(i == j && d1[1] == 'a'){
                            hess[k+nfact,k+nfact] <- sum(Thetas[,k]^2 *
                                (-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2))
                            next		
                        }				
                        if(d1[1] == 'a' && d2[1] == 'a'){
                            hess[i,j] <- hess[j,i] <- -sum(Thetas[,k]*Thetas[,m]*Qk*Qm*(f-r)*P/Q^2) 
                            next
                        }
                        if(d1[1] == 'd' && d2[1] == 'd'){
                            hess[i,j] <- hess[j,i] <- -sum(Qk*Qm*(f-r)*P/Q^2)
                            next
                        }	
                        if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                            hess[i,j] <- hess[j,i] <- sum(Thetas[,k]*Qk*(-Pk*(r - (f-r)*P/Q) - 
                                Qk*(f-r)*P/Q^2))
                            next	
                        }
                        if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                            hess[i,j] <- hess[j,i] <- -sum(Qk*Qm*Thetas[,m]*(f-r)*P/Q^2)
                            next
                        }						
                    }
                }
            }		
        }	
        return(list(grad = c(dd,da), hess = hess)) 
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
P.comp <- function(a, d, Theta, c = 0, u = 1)
{
    nfact <- length(a)
    P <- rep(1,nrow(Theta))
    for(i in 1:nfact)
        P <- P * P.mirt(a[i], d[i], Theta[ ,i, drop=FALSE],0)
    P <- c + (u - c) * P
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

