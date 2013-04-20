#----------------------------------------------------------------------------
# Derivatives wrt item parameters

setMethod(
    f = "Deriv",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){                 
        f <- 1
        dat <- x@dat[ ,2]
        Prior <- rep(1, length(dat))
        if(EM){            
            dat <- x@rs[,2] 
            f <- rowSums(x@rs)
            Prior <- rep(1, length(dat))
            if(BFACTOR) Prior <- prior
        }
        nfact <- x@nfact
        parlength <- length(x@par)
        u <- x@par[parlength]
        g <- x@par[parlength - 1]
        d <- x@par[parlength - 2]
        a <- x@par[1:nfact]        
        P <- P.mirt(a,d,Theta,g,u, D=x@D)  
        Q <- 1 - P        
        hess <- matrix(0,nfact+3, nfact+3)						                
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))       
        Pstar <- P.mirt(a,d,Theta,0,1, D=x@D)    	
        Qstar <- 1 - Pstar
        D <- x@D
        da <- rep(0,nfact)	
        r1 <- dat
        r2 <- f - dat
        r1_P <- r1/P; r1_P2 <- r1/P^2;
        r2_Q <- r2/Q; r2_Q2 <- r2/Q^2;
        dd <- sum((u-g)*D*Pstar*Qstar*(r1_P - r2_Q)*Prior)
        dc <- sum(Qstar*(r1_P - r2_Q)*Prior)
        du <- sum(Pstar*(r1_P - r2_Q)*Prior)
        for(i in 1:nfact)
            da[i] <- sum(Theta[,i]*D*Pstar*Qstar*(u-g)*(r1_P - r2_Q)*Prior)            
        grad <- c(da,dd,dc,du)
        #handy collections
        gloc <- nfact+2; uloc <- nfact+3
        ugD2 <- (u-g) * D^2; ugD <- (u-g) * D; Pstar2 <- Pstar^2; Pstar3 <- Pstar^3        
        ###
        hess[nfact+1,nfact+1] <- sum((r1_P * (ugD2 * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                                          r1_P2 * (ugD * (Pstar - Pstar2))^2 +  
                                          r2_Q * (ugD2 * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                                          r2_Q2 * (ugD * (-Pstar + Pstar2))^2)*Prior)
        hess[gloc,gloc] <- -sum(Qstar^2 *(r1_P2 + r2_Q2)*Prior) 
        hess[uloc,uloc] <- -sum(Pstar2 *(r1_P2 + r2_Q2)*Prior) 
        hess[nfact + 1, gloc] <- hess[gloc, nfact + 1] <- 
            sum((r1_P * (D * (-Pstar + Pstar2)) -
                     r1_P2 * (ugD * (Pstar - Pstar2)) * Qstar +  
                     r2_Q * (D * (Pstar - Pstar2)) -
                     r2_Q2 * (ugD * (-Pstar + Pstar2)) * -Qstar) * Prior)
        hess[nfact + 1, uloc] <- hess[uloc, nfact + 1] <-
            sum((r1_P * (D * (Pstar - Pstar2)) -
                     r1_P2 * (ugD * (Pstar - Pstar2)) * Pstar +  
                     r2_Q * (D * (-Pstar + Pstar2)) +
                     r2_Q2 * (ugD * (-Pstar + Pstar2)) * Pstar) * Prior)
        hess[gloc, uloc] <- hess[uloc, gloc] <-
            sum((-r1_P2 * Pstar * Qstar + r2_Q2 * Pstar * (-1 + Pstar )) * Prior)
        for(i in 1:nfact)
            for(j in 1:nfact)
                if(i <= j)
                    hess[i, j] <- sum((r1_P * (ugD2 * Theta[,i] * Theta[,j] * 
                                                     (Pstar - 3*Pstar2 + 2*Pstar3)) -
                                           r1_P2 * (ugD * Theta[,i] * (Pstar - Pstar2)) * 
                                              (ugD * Theta[,j] * (Pstar - Pstar2)) +  
                                           r2_Q * (ugD2 * Theta[,i] * Theta[,j] * 
                                               (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                                           r2_Q2 * (ugD * Theta[,i] * (-Pstar + Pstar2)) * 
                                              (ugD * Theta[,j] * (-Pstar + Pstar2))) * Prior)
        for(i in 1:nfact){
            hess[i, nfact + 1] <- hess[nfact + 1, i] <- 
                sum((r1_P * (ugD2 * Theta[,i] * (Pstar - 3*Pstar2 + 2*Pstar3)) -
                         r1_P2 * (ugD * Theta[,i] * (Pstar - Pstar2)) * 
                            (ugD * (Pstar - Pstar2)) +  
                         r2_Q * (ugD2 * Theta[,i] * (-Pstar + 3*Pstar2 - 2*Pstar3)) -
                         r2_Q2 * (ugD * Theta[,i] * (-Pstar + Pstar2)) * 
                         (ugD * (-Pstar + Pstar2))) * Prior)
            hess[i, gloc] <- hess[gloc, i] <- 
                sum((r1_P * (D * Theta[,i] * (-Pstar + Pstar2)) -
                         r1_P2 * (ugD * Theta[,i] * (Pstar - Pstar2)) * Qstar +  
                         r2_Q * (D * Theta[,i] * (Pstar - Pstar2)) -
                         r2_Q2 * (ugD * Theta[,i] * (-Pstar + Pstar2) ) * (-1 + Pstar)) * Prior)
            hess[i, uloc] <- hess[uloc, i] <- 
                sum((r1_P * (D * Theta[,i] * (Pstar - Pstar2)) -
                         r1_P2 * (ugD * Theta[,i] * (Pstar - Pstar2)) * Pstar +  
                         r2_Q * (D * Theta[,i] * (-Pstar + Pstar2)) +
                         r2_Q2 * (ugD * Theta[,i] * (-Pstar + Pstar2) ) * Pstar) * Prior)
        }       
           
        ret <- DerivativePriors(x=x, grad=grad, hess=hess)       
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))        
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs
            Prior <- rep(1, nrow(dat))
            if(BFACTOR) Prior <- prior
        } 
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[-(1:nfact)]
        nd <- length(d)    			
        P <- P.poly(a, d,Theta, D=x@D)			    	
        ret <- .Call("dparsPoly", P, Theta, Prior * x@D, dat, nd) 
        grad <- c(ret$grad[-(1:nd)], ret$grad[1:nd])
        hess <- matrix(0,nfact+nd,nfact+nd)
        hess[1:nfact,1:nfact] <- ret$hess[-(1:nd),-(1:nd)]
        hess[(nfact+1):ncol(hess), (nfact+1):ncol(hess)] <- ret$hess[1:nd, 1:nd]
        hess[1:nfact, (nfact+1):ncol(hess)] <- hess[(nfact+1):ncol(hess), 1:nfact] <- 
            ret$hess[(nd+1):ncol(hess), 1:nd]
        hess <- hess * x@D
        ret <- DerivativePriors(x=x, grad=grad, hess=hess)       
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){         
        hess <- matrix(0, length(x@par), length(x@par))        
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs
            Prior <- rep(1, nrow(dat))
            if(BFACTOR) Prior <- prior
        } 
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- ExtractZetas(x)
        nzetas <- length(d)
        shiftind <- length(x@par)
        shift <- x@par[shiftind]        
        nd <- length(d)    			
        P <- P.poly(a, d + shift, Theta, D=x@D)			    	
        ret <- .Call("dparsPoly", P, Theta, Prior * x@D, dat, nd) 
        grad <- c(ret$grad[-(1:nd)], ret$grad[1:nd])
        hess <- matrix(0,nfact+nd,nfact+nd)
        hess[1:nfact,1:nfact] <- ret$hess[-(1:nd),-(1:nd)]
        hess[(nfact+1):ncol(hess), (nfact+1):ncol(hess)] <- ret$hess[1:nd, 1:nd]
        hess[1:nfact, (nfact+1):ncol(hess)] <- hess[(nfact+1):ncol(hess), 1:nfact] <- 
            ret$hess[(nd+1):ncol(hess), 1:nd]
        hess <- hess * x@D        
        dc <- numeric(1)
        D <- x@D
        D2 <- D^2
        Pfull <- P
        PQfull <- Pfull * (1-Pfull)
        P <- P.poly(a, d + shift, Theta, D=x@D, itemexp=TRUE)   
        rs <- dat
        for(i in 1:ncol(rs))
            dc <- dc + rs[,i]/P[,i] * D * (PQfull[,i] - PQfull[,i+1])         
        dc <- sum(dc * Prior)
        grad <- c(grad, dc)
        hess <- cbind(hess, rep(0, nrow(hess)))
        hess <- rbind(hess, rep(0, ncol(hess)))
        cind <- ncol(hess)
        ddc <- ddd <- numeric(length(Prior))        
        dda <- matrix(0, length(Prior), nfact)        
        for(i in 1:ncol(rs))
            ddc <- ddc + rs[,i]/P[,i] * D2 * (Pfull[,i] - 3*Pfull[,i]^2 + 2*Pfull[,i]^3 -  
                Pfull[,i+1] + 3*Pfull[,i+1]^2 - 2*Pfull[,i+1]^3) - 
                rs[,i]/P[,i]^2 * D2 * (PQfull[,i] - PQfull[,i+1])^2                
        hess[cind, cind] <- sum(ddc * Prior)
        for(i in 1:nzetas)
            hess[cind, nfact + i] <- hess[nfact + i, cind] <- 
                sum((rs[,i]/P[,i] * D2 * (-Pfull[,i+1] + 3*Pfull[,i+1]^2 - 2*Pfull[,i+1]^3) - 
                rs[,i]/P[,i]^2 * D2 * (PQfull[,i] - PQfull[,i+1]) * (-PQfull[,i+1]) +
                rs[,i+1]/P[,i+1] * D2 * (Pfull[,i+1] - 3*Pfull[,i+1]^2 + 2*Pfull[,i+1]^3) - 
                rs[,i+1]/P[,i+1]^2 * D2 * (PQfull[,i+1] - PQfull[,i+2]) * (PQfull[,i+1]))*Prior)
        for(j in 1:nfact){
            tmp <- 0
            for(i in 1:ncol(rs))
                    tmp <- tmp + (rs[,i]/P[,i] * D2 * Theta[,j] * 
                                      (Pfull[,i] - 3*Pfull[,i]^2 + 2*Pfull[,i]^3 -  
                                           Pfull[,i+1] + 3*Pfull[,i+1]^2 - 2*Pfull[,i+1]^3) - 
                             rs[,i]/P[,i]^2 * D2 * (PQfull[,i] - PQfull[,i+1]) * Theta[,j] *
                                  (PQfull[,i] - PQfull[,i+1]))
            hess[cind, j] <- hess[j, cind] <- sum(tmp * Prior)
        }
        ret <- DerivativePriors(x=x, grad=grad, hess=hess)       
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){
        #local derivative from previous version with small mod
        dpars.comp <- function(lambda,zeta,g,r,f,Thetas,D,Prior)
        {    
            nfact <- length(lambda)    
            pars <- c(zeta,lambda,g)
            pgrad <- function(pars, r, thetas, D, Prior){        
                nfact <- ncol(thetas)
                d <- pars[1:nfact]	
                a <- pars[(nfact+1):(length(pars)-1)]
                c <- pars[length(pars)]
                P <- P.comp(a,d,thetas,c,D=D)		
                Pstar <- P.comp(a,d,thetas,0,1,D=D)		
                Qstar <- 1 - Pstar
                Q <- 1 - P
                const1 <- (r/P - (f-r)/Q)
                dd <- da <- rep(0,nfact)		
                dc <- sum(Qstar*const1)
                for(i in 1:nfact){
                    Pk <- P.mirt(a[i],d[i],matrix(thetas[,i]),0,1,D=D)
                    Qk <- 1 - Pk
                    dd[i] <- sum((1-c)*Pstar*Qk*const1*Prior)
                    da[i] <- sum((1-c)*Pstar*Qk*thetas[,i]*const1*Prior)
                }
                return(c(dd,da,dc))
            }		
            phess <- function(pars, r, thetas, D, Prior){        
                nfact <- ncol(thetas)
                d <- pars[1:nfact]	
                a <- pars[(nfact+1):(length(pars)-1)]
                c <- pars[length(pars)]
                P <- P.comp(a,d,thetas,c, 1, D=D)		
                Pstar <- P.comp(a,d,thetas,0,1,D=D)		
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
                            Pk <- P.mirt(a[k],d[k],matrix(thetas[,k]),0,1,D=D)
                            Qk <- 1 - Pk	
                            Pm <- P.mirt(a[m],d[m],matrix(thetas[,m]),0,1,D=D)
                            Qm <- 1 - Pm									
                            if(i == j && d1[1] == 'd'){
                                hess[i,i] <- sum((1-c)*Pstar*Qk*(const1*((1-c)*Qk - Pk) - 
                                                                     Pstar*Qk*(1-c)*const2)*Prior)
                                next
                            }
                            if(i == j && d1[1] == 'a'){
                                hess[i,i] <- sum((1-c)*thetas[,k]^2*Pstar*Qk*(const1*((1-c)*Qk - Pk)
                                                                              - Pstar*Qk*(1-c)*const2)*Prior)
                                next		
                            }
                            if(i == j && d1[1] == 'c'){
                                hess[i,i] <- -sum(Qstar^2 * const2*Prior)
                                next		
                            }	
                            if(d1[1] == 'a' && d2[1] == 'a'){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*thetas[,m]*
                                                                  Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2)*Prior)
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'd'){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*Pstar*Qm*(const1 - Pstar*(1-c)*const2)*Prior)
                                next
                            }
                            if(d1[1] == 'a' && d2[1] == 'c'){
                                hess[i,j] <- hess[j,i] <- -sum(thetas[,k]*Pstar*Qk*(const1 + Qstar*(1-c)*const2)*Prior)
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'c'){
                                hess[i,j] <- hess[j,i] <- -sum(Pstar*Qk*(const1 + Qstar*(1-c)*const2)*Prior)
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*thetas[,k]*Pstar*Qk*(const1*((1-c)*Qk - Pk) - 
                                                                                             Pstar*Qk*(1-c)*const2)*Prior)
                                next	
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                                hess[i,j] <- hess[j,i] <- sum((1-c)*Qk*thetas[,m]*Pstar*Qm*(const1 - 
                                                                                                Pstar*(1-c)*const2)*Prior)
                                next
                            }						
                        }
                    }
                }	
                return(hess)
            }
            #old pars in the form d, a, g
            g <- pgrad(pars, r, Thetas, D=D, Prior=Prior)
            h <- phess(pars, r, Thetas, D=D, Prior=Prior*D)    
            
            #translate into current version
            grad <- c(g[(nfact+1):(nfact*2)], g[1:nfact], g[length(g)], 0)
            hess <- matrix(0, ncol(h) + 1, ncol(h) + 1)
            hess[1:nfact, 1:nfact] <- h[(nfact+1):(nfact*2),(nfact+1):(nfact*2)] #a block
            hess[(nfact+1):(nfact*2),(nfact+1):(nfact*2)] <- h[1:nfact, 1:nfact] #d block
            hess[nfact*2 + 1, nfact*2 + 1] <- h[nfact*2 + 1, nfact*2 + 1] #g
            hess[nfact*2 + 1, 1:nfact] <- hess[1:nfact, nfact*2 + 1] <- 
                h[nfact*2 + 1, (nfact+1):(nfact*2)] #ga
            hess[nfact*2 + 1, (nfact+1):(nfact*2)] <- hess[(nfact+1):(nfact*2), nfact*2 + 1] <- 
                h[nfact*2 + 1, 1:nfact] #gd
            hess[(nfact+1):(nfact*2), 1:nfact] <- t(h[(nfact+1):(nfact*2), 1:nfact])
            hess[1:nfact, (nfact+1):(nfact*2)] <- t(h[1:nfact, (nfact+1):(nfact*2)]) #ads
            
            return(list(grad=grad, hess=hess))    
        }
        #####
        f <- 1
        r <- x@dat[ ,2]
        Prior <- rep(1, length(r))
        if(EM){       
            r <- x@rs[,2]
            f <- rowSums(x@rs)
            Prior <- rep(1, length(r))
            if(BFACTOR) Prior <- prior
        }
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(nfact*2)]
        g <- x@par[length(x@par)-1]        
                        
        tmp <- dpars.comp(lambda=ExtractLambdas(x),zeta=ExtractZetas(x),g=x@par[nfact*2 + 1],r=r, f=f,
                                   Thetas=Theta, D=x@D, Prior=Prior*x@D)
        ret <- DerivativePriors(x=x, grad=tmp$grad, hess=tmp$hess)
        return(ret)                   
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){     
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs
            Prior <- rep(1, nrow(dat))
            if(BFACTOR) Prior <- prior
        } 
        nfact <- x@nfact
        nzetas <- ncol(dat)
        a <- ExtractLambdas(x)                
        d <- ExtractZetas(x)
        ak <- 0:(length(d)-1)
        D <- x@D
        P <- ProbTrace(x=x, Theta=Theta)
        num <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D, returnNum=TRUE)         
        tmp <- nominalParDeriv(a=a, ak=ak, d=d, Theta=Theta, 
                               D=D, Prior=Prior, P=P, num=num, dat=dat, gpcm=TRUE)
        keep <- rep(TRUE, length(tmp$grad))
        keep[(nfact+1):(nfact+length(d))] <- FALSE        
        ret <- DerivativePriors(x=x, grad=tmp$grad[keep], hess=tmp$hess[keep, keep])          
        return(ret)    
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){          
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))              
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs
            Prior <- rep(1, nrow(dat))
            if(BFACTOR) Prior <- prior
            if(estHess) 
                hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], 
                                                        obj=x, Theta=Theta, prior=Prior)     
        } else {
            hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)    
        } 
        nfact <- x@nfact        
        D <- x@D
        a <- x@par[1:x@nfact]
        d <- x@par[x@nfact+1]
        g <- x@par[x@nfact+2]
        u <- x@par[x@nfact+3]        
        ak <- x@par[(x@nfact+4):(x@nfact+4+x@ncat-2)]
        dk <- x@par[(length(x@par)-length(ak)+1):length(x@par)]        
        correct <- x@correctcat        
        Pd <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u, D=D)
        Qd <- 1 - Pd
        Pstar <- P.mirt(a=a, d=d, Theta=Theta, g=0, u=1, D=D)
        Qstar <- 1 - Pstar
        num <- P.nominal(a=rep(1, nfact), ak=ak, d=dk, Theta=Theta, D=D, returnNum=TRUE)  
        den <- rowSums(num)
        Pn <- num/den
        cdat <- dat[,correct]
        idat <- dat[,-correct]         
        nd <- ncol(idat)
        for(i in 1:nfact)
            grad[i] <- sum( (u-g) * D * Theta[,i] * Qstar * Pstar * Prior * (
                cdat / Pd - rowSums(idat/Qd)) )
        grad[nfact+1] <- sum( (u-g) * D * Qstar * Pstar * Prior * (
                cdat / Pd - rowSums(idat/Qd)) )
        grad[nfact+2] <- sum( Prior * ((cdat * (1-Pstar)/Pd) + rowSums(idat * (Pstar - 1)/Qd)) )
        grad[nfact+3] <- sum( Prior * (cdat * Pstar / Pd - rowSums(idat * Pstar / Qd) ))
        for(j in 1:nd){
            grad[nfact+3+j] <- sum(Prior *(
                (idat[,j] * Qd * D * rowSums(Theta) * (Pn[,j] - Pn[,j]^2) * den) / (Qd * num[,j]) - 
                    rowSums(idat[,-j]) * D * rowSums(Theta) * Pn[,j]))
            grad[nfact+3+nd+j] <- sum(Prior *(
                (idat[,j] * Qd * D * (Pn[,j] - Pn[,j]^2) * den) / (Qd * num[,j]) - 
                    rowSums(idat[,-j]) * D * Pn[,j]))
        }
        ret <- DerivativePriors(x=x, grad=grad, hess=hess)
        return(ret)
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs
            Prior <- rep(1, nrow(dat))
            if(BFACTOR) Prior <- prior
        } 
        nfact <- x@nfact
        nzetas <- ncol(dat)
        a <- ExtractLambdas(x)                
        d <- ExtractZetas(x)
        shift <- d[length(d)]
        dshift <- d <- d[-length(d)]
        dshift[-1] <- d[-1] + shift
        ak <- 0:(length(d)-1)
        D <- x@D
        P <- ProbTrace(x=x, Theta=Theta)
        num <- P.nominal(a=a, ak=ak, d=dshift, Theta=Theta, D=D, returnNum=TRUE)         
        tmp <- nominalParDeriv(a=a, ak=ak, d=dshift, Theta=Theta, 
                               D=D, Prior=Prior, P=P, num=num, dat=dat)
        keep <- rep(TRUE, length(tmp$grad))
        keep[(nfact+1):(nfact+length(d))] <- FALSE        
        grad <- c(tmp$grad[keep], 0)
        hess <- tmp$hess[keep, keep]
        hess <- cbind(hess, rep(0, nrow(hess)))
        hess <- rbind(hess, rep(0, ncol(hess)))
        
        #quick calcs for derivs
        nfact <- length(a)
        ncat <- length(d)
        akind <- nfact 
        dind <- nfact + ncat*2 + 1 #go backwards
        D2 <- D^2
        ak2 <- ak^2
        P2 <- P^2
        P3 <- P^3
        aTheta <- as.vector(Theta %*% a)
        aTheta2 <- aTheta^2
        dat_num <- dat/num
        numsum <- rowSums(num) 
        numD <- num %*% c(0, rep(1, ncol(num)-1)) * D                
        numakThetaD <- num %*% ak * D * Theta
        numD2 <- num %*% c(0, rep(1, ncol(num)-1)) * D2  
        numakThetaD2 <- num %*% ak * D2 * Theta
        ak0 <- ak
        ak0[1] <- 0
        cind <- length(grad)
        tmp <- 0
        for(i in 1:nzetas)
            tmp <- tmp + dat[,i]*numD^2 / numsum^2 - dat[,i]*numD2/numsum 
        hess[cind, cind] <- sum(tmp*Prior)
        for(j in 1:nzetas){
            tmp <- 0
            for(i in 1:nzetas)
                tmp <- tmp + dat[,i]*P[,j]*D*numD/numsum - dat[,i]*D2*P[,j]
            hess[cind, nfact+j] <- hess[nfact+j, cind] <- sum(tmp*Prior)                
        }
        for(j in 1:nfact){
            tmp <- 0
            for(i in 1:nzetas)
                tmp <- tmp + dat[,i]*numD*numakThetaD[,j]/numsum^2 - 
                    dat[,i]* (num %*% ak0*D2*Theta[,j])/numsum
            hess[cind, j] <- hess[j, cind] <- sum(tmp*Prior) 
        }
        ####
        #TEMP - can't seem to get the last value of the gradient quite right for some reason....
        x2 <- x
        x2@est <- c(rep(FALSE, length(x2@est)-1), TRUE)
        grad[x2@est] <- numDeriv::grad(EML, x@par[x2@est], obj=x2, Theta=Theta, prior=Prior)        
        ####
        ret <- DerivativePriors(x=x, grad=grad, hess=hess)
        ret
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){                 
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs
            Prior <- rep(1, nrow(dat))
            if(BFACTOR) Prior <- prior
        } 
        nfact <- x@nfact
        nzetas <- ncol(dat)
        a <- ExtractLambdas(x)        
        ak <- x@par[(nfact+1):(nzetas + nfact)]
        d <- ExtractZetas(x)
        D <- x@D
        P <- ProbTrace(x=x, Theta=Theta)                        
        num <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D, returnNum=TRUE)                 
        tmp <- nominalParDeriv(a=a, ak=ak, d=d, Theta=Theta, 
                               D=D, Prior=Prior, P=P, num=num, dat=dat)        
        ret <- DerivativePriors(x=x, grad=tmp$grad, hess=tmp$hess)          
        return(ret)
    }
)

nominalParDeriv <- function(a, ak, d, Theta, D, Prior, P, num, dat, gpcm = FALSE){     
    nfact <- length(a)
    ncat <- length(d)
    akind <- nfact 
    dind <- nfact + ncat
    D2 <- D^2
    ak2 <- ak^2
    P2 <- P^2
    P3 <- P^3
    aTheta <- as.vector(Theta %*% a)
    aTheta2 <- aTheta^2
    dat_num <- dat/num
    numsum <- rowSums(num) 
    numakD <- num %*% ak * D
    numak2D2 <- num %*% ak2 * D2
    numakDTheta_numsum <- matrix(0, nrow(num), nfact)
    for(i in 1:nfact)
        numakDTheta_numsum[,i] <- (num %*% ak * D * Theta[, i])/ numsum    
    ret <- .Call('dparsNominal', a, ak, d, Theta, D, Prior, P, num, dat, nfact, ncat, 
                 akind, dind, ak2, P2, P3, aTheta, aTheta2, dat_num, numsum, numakD, 
                 numak2D2, numakDTheta_numsum)      
    ret
}

setMethod(
    f = "Deriv",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        Prior <- rep(1, nrow(x@rs))        
        if(BFACTOR) Prior <- prior
        if(EM){                
            grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=Prior)
            if(estHess) hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=Prior)     
            return(list(grad = grad, hess=hess))            
        }        
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, pars = NULL, itemloc = NULL, 
                          tabdata = NULL, prior=NULL, estHess=FALSE){
        if(EM){
            grad <- rep(0, length(x@par))
            hess <- matrix(0, length(x@par), length(x@par))            
            if(estHess){
                if(any(x@est)){
                    #grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, pars=pars, tabdata=tabdata,
                    #                              itemloc=itemloc)
                    hess[x@est,x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, pars=pars, tabdata=tabdata,
                                                           itemloc=itemloc)                                      
                }            
                return(list(grad=grad, hess=hess))
            }
            J <- length(pars) - 1
            nfact <- x@nfact
            scores <- matrix(0, nrow(tabdata), nfact)                 
            r <- tabdata[ ,ncol(tabdata)]
            N <- sum(r)
            tabdata <- tabdata[ ,-ncol(tabdata)]            
            itemtrace <- computeItemtrace(pars=pars, Theta=Theta, itemloc=itemloc)            
            mu <- x@par[1:nfact]
            siglong <- x@par[-(1:nfact)]
            sig <- matrix(0,nfact,nfact)                
            selcov <- lower.tri(sig, diag=TRUE)
            scores2 <- matrix(0, nrow(tabdata), sum(selcov))
            thetas2 <- numeric(sum(selcov))
            log_itemtrace <- log(itemtrace)                
            ret <- .Call('EAPgroup', log_itemtrace, tabdata, Theta, prior, mu)                
            tmp <- cbind(ret$scores, ret$scores2) * r
            newpars <- apply(tmp, 2, sum) / N                
            return(newpars[x@est])                                                        
        }
        tr <- function(y) sum(diag(y))         
        nfact <- x@nfact
        N <- nrow(Theta)
        u <- x@par[1:nfact]
        MU <- matrix(rep(u, N), N, byrow = TRUE)
        siglong <- x@par[-(1:nfact)]
        sig <- matrix(0,nfact,nfact)
        selcov <- lower.tri(sig, diag=TRUE) 
        sig[selcov] <- siglong
        if(nfact != 1)
            sig <- sig + t(sig) - diag(diag(sig))
        npars <- length(sig) + nfact	        
        invSig <- solve(sig)	
        Z <- t(Theta-MU) %*% (Theta-MU)
        g1 <- N * invSig %*% (colMeans(Theta) - u) 		
        tmp <- invSig %*% (Z - N * sig) %*% invSig         
        diag(tmp) <- diag(tmp)/2 #correct for symmetry
        g2 <- tmp[selcov]        
        grad <- c(g1, g2)
        sel <- 1:npars		
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
        #prior parameter constraints
        if(any(!is.nan(x@n.prior.mu))){            
            ind <- !is.na(x@n.prior.mu)            
            val <- x@par[ind]
            mu <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            h <- g <- rep(0, length(val))
            for(i in 1:length(val)){
                g[i] <- -(val[i] - mu[i])/(s[i]^2)
                h[i] <- -1/(s[i]^2)
            }
            grad[ind] <- grad[ind] + g
            if(length(val) == 1) hess[ind, ind] <- hess[ind, ind] + h
            else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + h           
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.na(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            bphess <- bpgrad <- rep(0, length(val))
            for(i in 1:length(val)){
                tmp <- betaprior(val[i], a[i], b[i])
                bpgrad[i] <- tmp$grad
                bphess[i] <- tmp$hess
            }
            if(length(val) == 1) hess[ind, ind] <- hess[ind, ind] + bpgrad
            else diag(hess[ind, ind]) <- diag(hess[ind, ind]) + bphess                
        }
        return(list(hess=hess,grad=grad))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL, estHess = FALSE){
        if(x@useuserdata) Theta <- cbind(Theta, x@userdata)
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        Prior <- rep(1, nrow(x@rs))        
        if(BFACTOR) Prior <- prior
        if(EM){
            if(x@usegr) grad <- x@gr(x, Theta, BFACTOR = FALSE, prior = NULL)
            else grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=Prior)
            if(estHess){ 
                if(x@usehss) hess <- x@hss(x, Theta, BFACTOR = FALSE, prior = NULL)
                else hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, 
                                                             Theta=Theta, prior=Prior)     
            }
            return(list(grad = grad))            
        }        
        if(x@usegr) grad <- x@gr(x, Theta, BFACTOR = FALSE, prior = NULL)
        else grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)     
        if(x@usehss) hess <- x@hss(x, Theta, BFACTOR = FALSE, prior = NULL)
        else hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)
        return(list(grad=grad, hess=hess))
    }
)       

#----------------------------------------------------------------------------
#TEMPORARY, until i calculate the analytical derivatives sometime
L <- function(par, obj, Theta){
    obj@par[obj@est] <- par
    P <- ProbTrace(obj, Theta)
    LL <- obj@dat * P
    LL[LL < .Machine$double.eps] <- 1
    LL <- sum(log(LL))
    if(any(!is.nan(obj@n.prior.mu))){
        ind <- !is.nan(obj@n.prior.mu)
        val <- obj@par[ind]
        u <- obj@n.prior.mu[ind]
        s <- obj@n.prior.sd[ind]
        for(i in 1:length(val))            
            LL <- LL + log(dnorm(val[i], u[i], s[i]))
    }
    if(any(!is.nan(obj@b.prior.alpha))){
        ind <- !is.nan(obj@b.prior.alpha)
        val <- obj@par[ind]
        a <- obj@b.prior.alpha[ind]
        b <- obj@b.prior.beta[ind]
        for(i in 1:length(val)){            
            tmp <- dbeta(val[i], a[i], b[i])
            LL <- LL + log(ifelse(tmp == 0, 1, tmp))
        }
    }    
    return(LL)
}

EML <- function(par, obj, Theta, ...){    
    obj@par[obj@est] <- par
    L <- (-1)*LogLik(x=obj, Theta=Theta, EM=TRUE, ...)
    return(L)
}


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
        u <- x@par[parlength]
        g <- x@par[parlength - 1]
        d <- x@par[parlength - 2]
        a <- x@par[1:nfact]   
        D <- x@D        
        Pstar <- P.mirt(a, d, Theta, g=0, u=u, D=x@D) 
        grad <- hess <- vector('list', 2)
        grad[[1]] <- grad[[2]] <- hess[[1]] <- hess[[2]] <- matrix(0, N, nfact)
        for(i in 1:nfact){            
            grad[[2]][ ,i] <- (u-g) * D * a[i] * (Pstar * (1 - Pstar))
            grad[[1]][ ,i] <- -1 * grad[[2]][ ,i]            
            hess[[2]][ ,i] <- 2 * (u - g) * D^2 * a[i]^2 * ((1 - Pstar)^2 * Pstar) - 
                (u - g) * D^2 * a[i]^2 * (Pstar * (1 - Pstar))
            hess[[1]][ ,i] <- -1 * hess[[2]][ ,i]
        }
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- ExtractLambdas(x)
        D <- x@D
        P <- ProbTrace(x, Theta, itemexp = FALSE) 
        grad <- hess <- vector('list', x@ncat)
        for(i in 1:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1:x@nfact){
            for(i in 1:(ncol(P)-1)){
                w1 <- P[,i] * (1-P[,i]) * a[j] * D
                w2 <- P[,i+1] * (1-P[,i+1]) * a[j] * D
                grad[[i]][ ,j] <- w1 - w2                
                hess[[i]][ ,j] <- D^2 * a[j]^2 * (2 * P[ ,i] * (1 - P[,i])^2 - 
                                                P[ ,i] * (1 - P[,i]) - 
                                                2 * P[ ,i+1] * (1 - P[,i+1])^2 + 
                                                P[ ,i+1] * (1 - P[,i+1]))
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
        D <- x@D
        P <- ProbTrace(x, Theta, itemexp = FALSE) 
        grad <- hess <- vector('list', x@ncat)
        for(i in 1:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1:x@nfact){
            for(i in 1:(ncol(P)-1)){
                w1 <- P[,i] * (1-P[,i]) * a[j] * D 
                w2 <- P[,i+1] * (1-P[,i+1]) * a[j] * D
                grad[[i]][ ,j] <- w1 - w2                
                hess[[i]][ ,j] <- D^2 * a[j]^2 * (2 * P[ ,i] * (1 - P[,i])^2 - 
                                                P[ ,i] * (1 - P[,i]) - 
                                                2 * P[ ,i+1] * (1 - P[ ,i+1])^2 + 
                                                P[ ,i+1] * (1 - P[ ,i+1]))
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
        g <- x@par[parlength - 1]
        d <- ExtractZetas(x)
        a <- ExtractLambdas(x)
        D <- x@D
        P <- P.comp(a, d, Theta, g, 1, D=D)
        Pdich <- P.mirt(a, d, Theta, 0, 1, D=D)
        Pstar <- P - g
        grad <- hess <- vector('list', 2)
        grad[[1]] <- grad[[2]] <- hess[[1]] <- hess[[2]] <- matrix(0, N, nfact)
        for(j in 1:nfact){            
            grad[[2]][ ,j] <- (u - g) * D * a[j] * Pstar * (1 - Pdich) 
            grad[[1]][ ,j] <- 1 - grad[[2]][ ,j]
            hess[[2]][ ,j] <- (u - g) * D^2 * a[j]^2 * ( 2 * (1 - Pdich)^2 * Pstar - 
                                                             (1 - Pdich) * Pstar)
            hess[[1]][ ,j] <- 1 - hess[[2]][ ,j]                
        }
        stop('DerivTheta for class \'', class(x), '\' not yet written.')    
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){
        D <- x@D
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- 0:(x@ncat - 1)        
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D)
        Num <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D, returnNum = TRUE)
        Den <- rowSums(Num)        
        grad <- hess <- vector('list', x@ncat)
        for(i in 1:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1:x@nfact){            
            for(i in 1:x@ncat){                
                grad[[i]][ ,j] <- D * ak[i] * a[j] * P[ ,i] - P[ ,i] * (Num %*% (D * ak * a[j])) / Den               
                hess[[i]][ ,j] <- D^2 * ak[i]^2 * a[j]^2 * P[ ,i] - 
                    2 * D * ak[i] * a[j] * P[,i] * (Num %*% (D * ak * a[j])) / Den + 
                    2 * P[,i] * ((Num %*% (D * ak * a[j])) / Den)^2 - 
                    P[,i] * ((Num %*% (D^2 * ak^2 * a[j]^2)) / Den)
            }    
        }               
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta){
        D <- x@D
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        t <- d[length(d)]
        d <- d[-length(d)]
        d[-1] <- d[-1] + t
        ak <- 0:(x@ncat - 1)        
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D)
        Num <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D, returnNum = TRUE)
        Den <- rowSums(Num)        
        grad <- hess <- vector('list', x@ncat)
        for(i in 1:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1:x@nfact){            
            for(i in 1:x@ncat){                
                grad[[i]][ ,j] <- D * ak[i] * a[j] * P[ ,i] - P[ ,i] * (Num %*% (D * ak * a[j])) / Den               
                hess[[i]][ ,j] <- D^2 * ak[i]^2 * a[j]^2 * P[ ,i] - 
                    2 * D * ak[i] * a[j] * P[,i] * (Num %*% (D * ak * a[j])) / Den + 
                    2 * P[,i] * ((Num %*% (D * ak * a[j])) / Den)^2 - 
                    P[,i] * ((Num %*% (D^2 * ak^2 * a[j]^2)) / Den)
            }    
        }               
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){
        D <- x@D
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- x@par[(x@nfact+1):(x@nfact+x@ncat)]
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D)
        Num <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=D, returnNum = TRUE)
        Den <- rowSums(Num)        
        grad <- hess <- vector('list', x@ncat)
        for(i in 1:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)
        for(j in 1:x@nfact){            
            for(i in 1:x@ncat){                
                grad[[i]][ ,j] <- D * ak[i] * a[j] * P[ ,i] - P[ ,i] * (Num %*% (D * ak * a[j])) / Den               
                hess[[i]][ ,j] <- D^2 * ak[i]^2 * a[j]^2 * P[ ,i] - 
                    2 * D * ak[i] * a[j] * P[,i] * (Num %*% (D * ak * a[j])) / Den + 
                    2 * P[,i] * ((Num %*% (D * ak * a[j])) / Den)^2 - 
                    P[,i] * ((Num %*% (D^2 * ak^2 * a[j]^2)) / Den)
            }    
        }              
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta){        
        stop('Derivatives for mcm items not yet written')
    }
)

setMethod(
    f = "DerivTheta",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta){                  
        D <- x@D
        a <- x@par[1:x@nfact]
        d <- x@par[x@nfact+1]
        g <- x@par[x@nfact+2]
        u <- x@par[x@nfact+3]        
        ak <- x@par[(x@nfact+4):(x@nfact+4+x@ncat-2)]
        dk <- x@par[(length(x@par)-length(ak)+1):length(x@par)]
        Pn <- P.nominal(a=rep(1,ncol(Theta)), ak=ak, d=dk, Theta=Theta, D=D)
        Num <- P.nominal(a=rep(1,ncol(Theta)), ak=ak, d=dk, Theta=Theta, D=D, returnNum = TRUE)
        Den <- rowSums(Num)        
        Pstar <- P.mirt(a, d, Theta, g=0, u=u, D=x@D)         
        Q <- 1 - P.mirt(a, d, Theta, g=g, u=u, D=x@D)         
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
        for(i in 1:x@ncat)
            grad[[i]] <- hess[[i]] <- matrix(0, nrow(Theta), x@nfact)            
        for(j in 1:x@nfact){            
            for(i in 1:x@ncat){
                if(i == x@correctcat){
                    grad[[i]][ ,j] <- (u-g) * D * a[j] * (Pstar * (1 - Pstar))
                    hess[[i]][ ,j] <- 2 * (u - g) * D^2 * a[j]^2 * ((1 - Pstar)^2 * Pstar) - 
                        (u - g) * D^2 * a[j]^2 * (Pstar * (1 - Pstar))
                } else {                    
                    grad[[i]][ ,j] <- -(u-g) * D * a[j] * (Pstar * (1 - Pstar)) * P[,i] + 
                        Q * (D * ak[i] * P[ ,i] - P[ ,i] * (Num %*% (D * ak)) / Den)
                    hess[[i]][ ,j] <-
                        -2 * (u - g) * D^2 * a[j]^2 * (1 - Pstar)^2 * Pstar * P[,i] +
                        (u - g) * D^2 * a[j]^2 * Pstar * (1 - Pstar) * P[,i] - 
                        2 * (u - g) * D^2 * a[j] * ak[i] * (1 - Pstar) * Pstar * P[,i] +
                        2 * D * a[j] *  (Pstar * (1 - Pstar)) * P[,i] * (Num %*% (D * ak)) / Den +
                        Q * D^2 * ak[i]^2 * P[ ,i] - 
                        2 * Q * D * ak[i] * P[,i] * (Num %*% (D * ak)) / Den + 
                        2 * Q * P[,i] * ((Num %*% (D * ak)) / Den)^2 - 
                        Q * P[,i] * ((Num %*% (D^2 * ak^2)) / Den)
                }
            }    
        }  
        return(list(grad=grad, hess=hess))                
    }
)

