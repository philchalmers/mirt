#----------------------------------------------------------------------------
setMethod(
    f = "print",
    signature = signature(x = 'dich'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "print",
    signature = signature(x = 'graded'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "print",
    signature = signature(x = 'rating'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "print",
    signature = signature(x = 'gpcm'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "print",
    signature = signature(x = 'nominal'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "print",
    signature = signature(x = 'partcomp'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "print",
    signature = signature(x = 'mcm'),
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
    f = "show",
    signature = signature(object = 'graded'),
    definition = function(object){
        print(object)
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
    f = "show",
    signature = signature(object = 'gpcm'),
    definition = function(object){
        print(object)
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
    f = "show",
    signature = signature(object = 'partcomp'),
    definition = function(object){
        print(object)
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'mcm'),
    definition = function(object){
        print(object)
    }
)

#----------------------------------------------------------------------------
#Probability Traces
setMethod(
    f = "ProbTrace",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta){  
        u <- x@par[length(x@par)]
        g <- x@par[length(x@par)-1]
        d <- x@par[length(x@par)-2]
        a <- x@par[1:x@nfact]        
        P <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u)
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE){                  
        a <- x@par[1:x@nfact]        
        d <- x@par[(x@nfact+1):length(x@par)]
        P <- P.poly(a=a, d=d, Theta=Theta, itemexp=itemexp)
        return(P)
    }
)


setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE){
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(length(x@par)-1)]
        t <- x@par[length(x@par)]
        P <- P.poly(a=a, d=(d + t), Theta=Theta, itemexp=itemexp)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){                  
        a <- x@par[1:x@nfact]        
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
        P <- P.comp(a=a, d=d, Theta=Theta, g=g, u=u)
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta){    
        a <- x@par[1:x@nfact]        
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
    definition = function(x, Theta, EM=FALSE, prior=NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)        
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        itemtrace[itemtrace < 1e-8] <- 1e-8
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
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
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta, pars, tabdata, itemloc, EM = TRUE){
        r <- tabdata[, ncol(tabdata)]
        gpars <- ExtractGroupPars(x)
        mu <- gpars$gmeans
        sigma <- gpars$gcov
        prior <- mvtnorm::dmvnorm(Theta, mean=mu, sigma=sigma)
        prior <- prior/sum(prior)
        rlist <- Estep.mirt(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior, itemloc=itemloc, 
                            debug='')
        LL <- (-1)*sum(r*log(rlist$expected))                 
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
    signature = signature(x = 'rating'),
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
    signature = signature(x = 'rating'),
    definition = function(x){          
        par <- x@par
        d <- par[-c(1:x@nfact, length(par))]
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
    signature = signature(x = 'dich', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){          
        P <- ProbTrace(x, Theta)[,2]
        nfact <- ncol(Theta)
        a <- ExtractLambdas(x)
        A <- sum((a * cosangle)^2)
        Pstar <- P.mirt(x@par[1:nfact], x@par[nfact + 1], Theta, 0, 1)
        info <- A * P * (1-P) * Pstar/P 
        info    
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'graded', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){
        P <- ProbTrace(x, Theta, itemexp = FALSE) 
        a <- ExtractLambdas(x)
        A <- sum((a * cosangle)^2)
        info <- 0
        for(i in 1:(ncol(P)-1)){
            w1 <- P[,i]*(1-P[,i]) * A
            w2 <- P[,i+1]*(1-P[,i+1]) * A
            info <- info + ((w1 - w2)^2) / (P[,i] - P[,i+1])             
        }    
        info
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'rating', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(length(x@par)-1)]
        t <- x@par[length(x@par)]
        ak <- seq(0, x@ncat-1, by = 1)
        info <- Info.nominal(Theta=Theta, a=a, ak=ak, d=(d+t), cosangle=cosangle)
        info        
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'gpcm', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- seq(0, x@ncat-1, by = 1)
        info <- Info.nominal(Theta=Theta, a=a, ak=ak, d=d, cosangle=cosangle)
        info
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'nominal', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){
        a <- ExtractLambdas(x)
        d <- ExtractZetas(x)
        ak <- x@par[(length(a)+1):(length(a) + length(d))]
        info <- Info.nominal(Theta=Theta, a=a, ak=ak, d=d, cosangle=cosangle)
        info
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'partcomp', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){
        stop('Information functions not yet written for ', class(x))
    }
)

setMethod(
    f = "ItemInfo",
    signature = signature(x = 'mcm', Theta = 'matrix', cosangle = 'numeric'),
    definition = function(x, Theta, cosangle = 1){
        stop('Information functions not yet written for ', class(x))
    }
)

#----------------------------------------------------------------------------
setMethod(
    f = "Deriv",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){        
        f <- 1
        dat <- x@dat[ ,2]
        Prior <- rep(1, length(dat))
        if(EM){
            Prior <- prior
            dat <- x@rs[,2] 
            f <- rowSums(x@rs)
        }
        nfact <- x@nfact
        parlength <- length(x@par)
        u <- x@par[parlength]
        g <- x@par[parlength - 1]
        d <- x@par[parlength - 2]
        a <- x@par[1:nfact]        
        P <- P.mirt(a, d, Theta, g, u)
        P[P < 1e-8] <- 1e-8
        P[P > .99999999] <- .99999999
        Q <- 1 - P        
        hess <- matrix(0,nfact+3, nfact+3)						
        if(x@par[parlength] == 1){ #'3PL'            	
            Pstar <- P.mirt(a,d,Theta,0,1)		
            Pstar[Pstar < 1e-8] <- 1e-8
            Pstar[Pstar > .99999999] <- .99999999
            Qstar <- 1 - Pstar
            da <- rep(0,nfact)	
            dd <- sum((1-g)*Pstar*Qstar*(dat/P - (f-dat)/Q)*Prior)
            dc <- sum(Qstar*(dat/P - (f-dat)/Q)*Prior)
            for(i in 1:nfact){
                da[i] <- sum(Theta[,i]*Pstar*Qstar*(1-g)*(dat/P - (f-dat)/Q)*Prior)
            }
            grad <- c(da,dd,dc,0)				
            gloc <- nfact+2
            const1 <- (dat/P - (f-dat)/Q)*(Qstar-Pstar)
            const2 <- (dat/P^2 + (f-dat)/Q^2)	
            hess[nfact+1,nfact+1] <- sum((1-g)*Pstar*Qstar*(const1 - 
                Pstar*Qstar*(1-g)*const2)*Prior)		
            hess[gloc,gloc] <- -sum(Qstar^2 *(dat/P^2 + (f-dat)/Q^2)*Prior)
            hess[gloc,nfact+1] <- hess[nfact+1,gloc] <- sum(-Pstar*Qstar*((dat/P - (f-dat)/Q) + 
                Qstar*(1-g)*const2)*Prior) 
            for(i in 1:nfact){
                hess[nfact+1,i] <- hess[i,nfact+1] <- sum((1-g)*Theta[,i]*Pstar*Qstar*(const1 - 
                    Pstar*Qstar*(1-g)*const2)*Prior)			
                hess[gloc,i] <- hess[i,gloc] <- sum(-Theta[,i]*Pstar*Qstar*((dat/P - (f-dat)/Q) + 
                    Qstar*(1-g)*const2)*Prior)		
                for(j in 1:nfact){
                    if(i == j)
                        hess[i,i] <- sum(Theta[,i]^2 *Pstar*Qstar*(1-g)*(const1 - 
                            (1-g)*Pstar*Qstar*const2)*Prior)
                    if(i < j)
                        hess[i,j] <- hess[j,i] <- sum(Theta[,i]*Theta[,j] *Pstar*Qstar*(1-g)*
                            (const1 - (1-g)*Pstar*Qstar*const2)*Prior)					
                }
            }	
        } else { #4PL 
            grad <- rep(0, length(x@par))
            hess <- matrix(0, length(x@par), length(x@par))
            if(EM){                
                grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)
                hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)       
                return(list(grad = grad, hess = hess))            
            }            
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){         
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))        
        dat <- x@dat 
        Prior <- rep(1, nrow(dat))
        if(EM){
            dat <- x@rs / sum(x@rs)                          
            Prior <- prior
        } 
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[-(1:nfact)]
        nd <- length(d)    			
        P <- P.poly(a, d,Theta)			    	
        ret <- .Call("dparsPoly", P, Theta, Prior, dat, nd) 
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
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){ 
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(EM){                
            grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)
            hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)       
            return(list(grad = grad, hess = hess))            
        }        
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)        
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){
        f <- 1
        r <- x@dat
        if(EM){       
            r <- x@rs[,2]
            f <- rowSums(x@rs)            
        }
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(nfact*2)]
        g <- x@par[length(x@par)-1]
        u <- x@par[length(x@par)]
        if(TRUE){#3PL
            grad <- function(a, d, g, u, r, f, Theta){                			
                P <- P.comp(a,d,Theta,g,1)		
                Pstar <- P.comp(a,d,Theta,0)		
                P[P < 1e-8] <- 1e-8
                P[P > .99999999] <- .99999999
                Pstar[Pstar < 1e-8] <- 1e-8
                Pstar[Pstar > .99999999] <- .99999999
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
            hess <- function(a, d, g, u, r, f, Theta){ 
                nfact <- length(a)
                hess <- matrix(0, nfact*2 + 2, nfact*2 + 2)                		
                P <- P.comp(a,d,Theta,g, 1)		
                Pstar <- P.comp(a,d,Theta,0, 1)
                P[P < 1e-8] <- 1e-8
                Pstart[Pstart < 1e-8] <- 1e-8
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
            grad <- grad(a=a, d=d, g=g, u=u, r=r, f=f, Theta=Theta)
            hess <- hess(a=a, d=d, g=g, u=u, r=r, f=f, Theta=Theta)
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
    definition = function(x, Theta, EM = FALSE, prior = NULL){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(EM){            
            grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)
            hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)       
            return(list(grad = grad, hess = hess))            
        }        
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)        
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(EM){                
            grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)
            hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)       
            return(list(grad = grad, hess = hess))            
        }        
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)        
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(EM){                
            grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)
            hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=prior)     
            return(list(grad = grad, hess = hess))            
        }        
        grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
        hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)
        return(list(grad = grad, hess = hess))
    }
)

setMethod(
    f = "Deriv",
    signature = signature(x = 'GroupPars', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, pars = NULL, itemloc = NULL, tabdata = NULL){
        if(EM){        
            grad <- rep(0, length(x@par))
            hess <- matrix(0, length(x@par), length(x@par))            
            if(any(x@est)){
                grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, pars=pars, tabdata=tabdata,
                                       itemloc=itemloc)
                hess[x@est,x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, pars=pars, tabdata=tabdata,
                                          itemloc=itemloc)                  
            }            
            return(list(grad = grad, hess = hess))            
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

#nominal and gpcm
Info.nominal <- function(Theta, a, ak, d, cosangle = 1){
    P <- P.nominal(a, ak, d, Theta)
    A <- a * cosangle
    nfact <- ncol(Theta)
    AK <- matrix(ak, nrow(Theta), length(ak), byrow = TRUE)
    M <- AK * P
    M2 <- AK^2 * P
    d2P <- dP <- matrix(0,nrow(AK), ncol(AK))
    info <- 0
    for(n in 1:nfact){
        for(i in 1:ncol(dP))        
            dP[,i] <- A[n] * P[,i] * (ak[i] - rowSums(M)) 
        for(i in 1:ncol(dP))        
            d2P[,i] <- ak[i]^2 * A[n]^2 * P[,i] - 
                2 * ak[i] * A[n]^2 * P[,i] * rowSums(M) + 
                2 * A[n]^2 * P[,i] * rowSums(M^2) - 
                A[n]^2 * P[,i] * rowSums(M2)
        info <- info + rowSums((dP)^2 / P - d2P)    
    }
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

EML <- function(par, obj, Theta, ...){    
    obj@par[obj@est] <- par
    L <- (-1)*LogLik(x=obj, Theta=Theta, EM=TRUE, ...)
    return(L)
}
