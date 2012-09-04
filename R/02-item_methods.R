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
    definition = function(x, Theta, itemexp = TRUE){                  
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        d <- x@par[(x@nfact+1):length(x@par)]
        P <- P.poly(a=a, d=d, Theta=Theta, itemexp=itemexp)
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
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta){    
        a <- x@par[1:x@nfact]
        if(x@bfactor) a <- a[x@est[1:x@nfact]]
        ind <- x@nfact + 1
        ak <- x@par[ind:(ind + x@ncat - 1)]
        ind <- ind + x@ncat
        d <- x@par[ind:(ind + x@ncat - 1)]
        ind <- ind + x@ncat
        t <- x@par[ind:length(x@par)]        
        P <- P.mcm(a=a, ak=ak, d=d, t=t, Theta=Theta)
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
        itemtrace[itemtrace < 1e-8] <- 1e-8
        LL <- (-1) * sum(x@rs * log(itemtrace))
        if(any(!is.nan(x@n.prior.mu))){
            ind <- !is.nan(x@n.prior.mu)
            val <- x@par[ind]
            u <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            for(i in 1:length(val))            
                LL <- LL - log(dnorm(val[i], u[i], s[i]))
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.nan(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            for(i in 1:length(val)){            
                tmp <- dbeta(val[i], a[i], b[i])
                LL <- LL - log(ifelse(tmp == 0, 1, tmp))
            }
        }        
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        LL <- (-1) * sum(x@rs * log(itemtrace))
        if(any(!is.nan(x@n.prior.mu))){
            ind <- !is.nan(x@n.prior.mu)
            val <- x@par[ind]
            u <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            for(i in 1:length(val))            
                LL <- LL - log(dnorm(val[i], u[i], s[i]))
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.nan(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            for(i in 1:length(val)){            
                tmp <- dbeta(val[i], a[i], b[i])
                LL <- LL - log(ifelse(tmp == 0, 1, tmp))
            }
        }
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        LL <- (-1) * sum(x@rs * log(itemtrace))
        if(any(!is.nan(x@n.prior.mu))){
            ind <- !is.nan(x@n.prior.mu)
            val <- x@par[ind]
            u <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            for(i in 1:length(val))            
                LL <- LL - log(dnorm(val[i], u[i], s[i]))
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.nan(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            for(i in 1:length(val)){            
                tmp <- dbeta(val[i], a[i], b[i])
                LL <- LL - log(ifelse(tmp == 0, 1, tmp))
            }
        }
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        LL <- (-1) * sum(x@rs * log(itemtrace))
        if(any(!is.nan(x@n.prior.mu))){
            ind <- !is.nan(x@n.prior.mu)
            val <- x@par[ind]
            u <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            for(i in 1:length(val))            
                LL <- LL - log(dnorm(val[i], u[i], s[i]))
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.nan(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            for(i in 1:length(val)){            
                tmp <- dbeta(val[i], a[i], b[i])
                LL <- LL - log(ifelse(tmp == 0, 1, tmp))
            }
        }
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        LL <- (-1) * sum(x@rs * log(itemtrace))
        if(any(!is.nan(x@n.prior.mu))){
            ind <- !is.nan(x@n.prior.mu)
            val <- x@par[ind]
            u <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            for(i in 1:length(val))            
                LL <- LL - log(dnorm(val[i], u[i], s[i]))
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.nan(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            for(i in 1:length(val)){            
                tmp <- dbeta(val[i], a[i], b[i])
                LL <- LL - log(ifelse(tmp == 0, 1, tmp))
            }
        }
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        LL <- (-1) * sum(x@rs * log(itemtrace))
        if(any(!is.nan(x@n.prior.mu))){
            ind <- !is.nan(x@n.prior.mu)
            val <- x@par[ind]
            u <- x@n.prior.mu[ind]
            s <- x@n.prior.sd[ind]
            for(i in 1:length(val))            
                LL <- LL - log(dnorm(val[i], u[i], s[i]))
        }
        if(any(!is.nan(x@b.prior.alpha))){
            ind <- !is.nan(x@b.prior.alpha)
            val <- x@par[ind]
            a <- x@b.prior.alpha[ind]
            b <- x@b.prior.beta[ind]
            for(i in 1:length(val)){            
                tmp <- dbeta(val[i], a[i], b[i])
                LL <- LL - log(ifelse(tmp == 0, 1, tmp))
            }
        }
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

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'mcm'),
    definition = function(x){          
        par <- x@par
        a <- par[1:x@nfact]
        a        
    }
)

#----------------------------------------------------------------------------
setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'dich'),
    definition = function(x){          
        par <- x@par
        d <- par[1:x@nfact]
        d        
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'graded'),
    definition = function(x){          
        par <- x@par
        d <- par[-(1:x@nfact)]
        d        
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'gpcm'),
    definition = function(x){          
        par <- x@par
        d <- par[-(1:x@nfact)]
        d        
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'nominal'),
    definition = function(x){          
        par <- x@par
        d <- x@par[length(x@par):(length(x@par) - x@ncat + 1)]
        d        
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'partcomp'),
    definition = function(x){          
        par <- x@par
        d <- x@par[(nfact+1):(length(x@par)-2)]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'mcm'),
    definition = function(x){          
        par <- x@par
        d <- x@par[(x@nfact + x@ncat +1):(x@nfact - x@ncat*2)]
        d        
    }
)

#----------------------------------------------------------------------------
setMethod(
    f = "ItemInfo",
    signature = signature(x = 'dich', A = 'numeric', Theta = 'matrix'),
    definition = function(x, A, Theta){          
        P <- ProbTrace(x, Theta)[,2]
        nfact <- ncol(Theta)
        Pstar <- P.mirt(x@par[1:nfact], x@par[nfact + 1], Theta, 0, 1)
        info <- A^2 * P * (1-P) * Pstar/P 
        info    
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'graded', A = 'numeric', Theta = 'matrix'),
    definition = function(x, A, Theta){          
        P <- ProbTrace(x, Theta, itemexp = FALSE)  
        info <- 0
        for(i in 1:(ncol(P)-1)){
            w1 <- P[,i]*(1-P[,i])*A
            w2 <- P[,i+1]*(1-P[,i+1])*A
            info <- info + ((w1 - w2)^2) / (P[,i] - P[,i+1])             
        }    
        info
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'gpcm', A = 'numeric', Theta = 'matrix'),
    definition = function(x, A, Theta){
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- seq(0, x@ncat-1, by = 1)
        info <- Info.nominal(Theta=Theta, a=a, ak=ak, A=A, d=d)
        info
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'nominal', A = 'numeric', Theta = 'matrix'),
    definition = function(x, A, Theta){          
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- x@par[(length(a)+1):(length(a) + length(d))]
        info <- Info.nominal(Theta=Theta, a=a, ak=ak, A=A, d=d)
        info
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'partcomp', A = 'numeric', Theta = 'matrix'),
    definition = function(x, A, Theta){          
        stop('Information functions not yet written for ', class(x))
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'mcm', A = 'numeric', Theta = 'matrix'),
    definition = function(x, A, Theta){          
        stop('Information functions not yet written for ', class(x))
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
        hess <- matrix(0,nfact+3, nfact+3)						
        if(!x@est[parlength] && !x@est[parlength-1]){ #'2PL'
            PQ <- P*Q
            L1 <- colSums((dat-P) * Theta)
            L2 <- sum(dat-P)
            grad <- c(L1,L2,0,0)    	
            L11 <- .Call("dichOuter", Theta, PQ, nrow(Theta))
            hess[1:nfact, 1:nfact] <- -L11
            hess[nfact+1, nfact+1] <- (-1)*sum(PQ)		
            hess[nfact+1, 1:nfact] <- hess[1:nfact, nfact+1] <- (-1)*colSums(PQ * Theta)
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
            grad <- c(da,dd,dc,0)				
            gloc <- nfact+2
            const1 <- (dat/P - (f-dat)/Q)*(Qstar-Pstar)
            const2 <- (dat/P^2 + (f-dat)/Q^2)	
            hess[nfact+1,nfact+1] <- sum((1-g)*Pstar*Qstar*(const1 - 
                Pstar*Qstar*(1-g)*const2))		
            hess[gloc,gloc] <- -sum(Qstar^2 *(dat/P^2 + (f-dat)/Q^2))
            hess[gloc,nfact+1] <- hess[nfact+1,gloc] <- sum(-Pstar*Qstar*((dat/P - (f-dat)/Q) + Qstar*(1-g)*const2)) 
            for(i in 1:nfact){
                hess[nfact+1,i] <- hess[i,nfact+1] <- sum((1-g)*Theta[,i]*Pstar*Qstar*(const1 - 
                    Pstar*Qstar*(1-g)*const2))			
                hess[gloc,i] <- hess[i,gloc] <- sum(-Theta[,i]*Pstar*Qstar*((dat/P - (f-dat)/Q) + 
                    Qstar*(1-g)*const2))		
                for(j in 1:nfact){
                    if(i == j)
                        hess[i,i] <- sum(Theta[,i]^2 *Pstar*Qstar*(1-g)*(const1 - 
                            (1-g)*Pstar*Qstar*const2))
                    if(i < j)
                        hess[i,j] <- hess[j,i] <- sum(Theta[,i]*Theta[,j] *Pstar*Qstar*(1-g)*
                            (const1 - (1-g)*Pstar*Qstar*const2))					
                }
            }	
        } else { #4PL
            grad <- rep(0, length(x@par))
            hess <- matrix(0, length(x@par), length(x@par))
            grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
            hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)            
            return(list(grad=grad, hess=hess))
        }       
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
        return(list(grad = grad, hess = hess))
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
        ret <- .Call("dparsPoly", P, Theta, x@dat, nd) #switch the order
        grad <- c(ret$grad[-(1:nd)], ret$grad[1:nd])
        hess <- matrix(0,nfact+nd,nfact+nd)
        hess[1:nfact,1:nfact] <- ret$hess[-(1:nd),-(1:nd)]
        hess[(nfact+1):ncol(hess), (nfact+1):ncol(hess)] <- ret$hess[1:nd, 1:nd]
        hess[1:nfact, (nfact+1):ncol(hess)] <- hess[(nfact+1):ncol(hess), 1:nfact] <- 
            ret$hess[(nd+1):ncol(hess), 1:nd]
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
        return(list(grad=grad, hess=hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta){
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(nfact*2)]
        g <- x@par[length(x@par)-1]
        u <- x@par[length(x@par)]
        if(x@est[nfact*2 + 1] && !x@est[nfact*2+2]){#3PL
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
                hess <- matrix(0, nfact*2 + 2, nfact*2 + 2)
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
                                hess[i,i] <- sum((1-g)*Pstar*Qk*(const1*((1-g)*Qk - Pk) - Pstar*Qk*(1-g)*const2))
                                next
                            }
                            if(i == j && d1[1] == 'a'){
                                hess[i,i] <- sum((1-g)*Theta[,k]^2*Pstar*Qk*(const1*((1-g)*Qk - Pk) - Pstar*Qk*
                                    (1-g)*const2))
                                next		
                            }
                            if(i == j && d1[1] == 'g'){
                                hess[i,i] <- -sum(Qstar^2 * const2)
                                next		
                            }	
                            if(d1[1] == 'a' && d2[1] == 'a'){
                                hess[i,j] <- hess[j,i] <- sum((1-g)*Theta[,k]*Theta[,m]*Qk*Pstar*Qm*(const1 - 
                                    Pstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'd'){
                                hess[i,j] <- hess[j,i] <- sum((1-g)*Qk*Pstar*Qm*(const1 - Pstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'a' && d2[1] == 'g'){
                                hess[i,j] <- hess[j,i] <- -sum(Theta[,k]*Pstar*Qk*(const1 + Qstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'g'){
                                hess[i,j] <- hess[j,i] <- -sum(Pstar*Qk*(const1 + Qstar*(1-g)*const2))
                                next
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                                hess[i,j] <- hess[j,i] <- sum((1-g)*Theta[,k]*Pstar*Qk*(const1*((1-g)*Qk - Pk) - 
                                    Pstar*Qk*(1-g)*const2))
                                next	
                            }
                            if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                                hess[i,j] <- hess[j,i] <- sum((1-g)*Qk*Theta[,m]*Pstar*Qm*(const1 - 
                                    Pstar*(1-g)*const2))
                                next
                            }						
                        }
                    }
                }	
                return(hess)
            }		
            grad <- grad(a, d, g, u, x@dat, Theta)
            hess <- hess(a, d, g, u, x@dat, Theta)
        }
        if(!x@est[nfact*2 + 1] && !x@est[nfact*2+2]){ #2PL
            dat <- x@dat
            hess <- matrix(0, nfact*2 + 2, nfact*2 + 2)
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
            grad <- c(da,dd,0,0)
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
                        Pk <- P.mirt(a[k],d[k],Theta[ , k, drop=FALSE],0,1)
                        Qk <- 1 - Pk	
                        Pm <- P.mirt(a[m],d[m],Theta[ , m, drop=FALSE],0,1)
                        Qm <- 1 - Pm									
                        if(i == j && d1[1] == 'd'){
                            hess[k,k] <- sum(-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2)
                            next
                        }
                        if(i == j && d1[1] == 'a'){
                            hess[k+nfact,k+nfact] <- sum(Theta[,k]^2 *
                                (-Pk*Qk*(r - (f-r)*P/Q) - Qk^2 * (f-r)*P/Q^2))
                            next		
                        }				
                        if(d1[1] == 'a' && d2[1] == 'a'){
                            hess[i,j] <- hess[j,i] <- -sum(Theta[,k]*Theta[,m]*Qk*Qm*(f-r)*P/Q^2) 
                            next
                        }
                        if(d1[1] == 'd' && d2[1] == 'd'){
                            hess[i,j] <- hess[j,i] <- -sum(Qk*Qm*(f-r)*P/Q^2)
                            next
                        }	
                        if(d1[1] == 'd' && d2[1] == 'a' && d1[2] == d2[2]){
                            hess[i,j] <- hess[j,i] <- sum(Theta[,k]*Qk*(-Pk*(r - (f-r)*P/Q) - 
                                Qk*(f-r)*P/Q^2))
                            next	
                        }
                        if(d1[1] == 'd' && d2[1] == 'a' && d1[2] != d2[2]){
                            hess[i,j] <- hess[j,i] <- -sum(Qk*Qm*Theta[,m]*(f-r)*P/Q^2)
                            next
                        }						
                    }
                }
            }                        
        }	
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
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)        
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)
        browser()
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)
        return(list(grad = grad, hess = hess))
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

#ak[1] == 0, ak[length(ak)] == length(ak) - 1, d[1] ==0, sum(t) == 1 
P.mcm <- function(a, ak, d, t, Theta){
    ncat <- length(d)
    nfact <- ncol(Theta)    
    a <- matrix(a)    
    P <- numerator <- matrix(0, nrow(Theta), ncat)  
    
    for(i in 1:ncat)
        numerator[ ,i] <- exp(1.702 * ak[i] * (Theta %*% a) + 1.702 * d[i])
    denominator <- rowSums(numerator)
    C0 <- 1 / denominator
    t[1] <- 1 - sum(t[2:length(t)])
    T <- matrix(t, nrow(P), ncat, byrow = TRUE)    
    P <- C0 * T + (1 - C0) * numerator/denominator
    return(P)   
}

#nominal/gpcm item info
Info.nominal <- function(Theta, a, ak, A, d){
    P <- P.nominal(a, ak, d, Theta)    
    AK <- matrix(ak, nrow(Theta), length(ak), byrow = TRUE)
    M <- AK * P
    M2 <- AK^2 * P
    d2P <- dP <- matrix(0,nrow(AK), ncol(AK))
    for(i in 1:ncol(dP))        
        dP[,i] <- A * P[,i] * (ak[i] - rowSums(M)) 
    for(i in 1:ncol(dP))        
        d2P[,i] <- ak[i]^2 * A^2 * P[,i] - 
            2 * ak[i] * A^2 * P[,i] * rowSums(M) + 
            2 * A^2 * P[,i] * rowSums(M^2) - 
            A^2 * P[,i] * rowSums(M2)
    info <- rowSums((dP)^2 / P - d2P)    
    info
}

#TEMPORARY, until i calculate the analytical derivatives sometime
L <- function(par, obj, Theta){
    obj@par[obj@est] <- par
    P <- ProbTrace(obj, Theta)
    LL <- obj@dat * P
    LL[LL < 1e-8] <- 1
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