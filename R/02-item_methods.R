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
        LL <- LL.Priors(x=x, LL=LL)        
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
        LL <- LL.Priors(x=x, LL=LL)
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
        LL <- LL.Priors(x=x, LL=LL)
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
        LL <- LL.Priors(x=x, LL=LL)
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
        LL <- LL.Priors(x=x, LL=LL)
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
        LL <- LL.Priors(x=x, LL=LL)
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
        LL <- LL.Priors(x=x, LL=LL)
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
        LL <- LL.Priors(x=x, LL=LL)
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
        d <- x@par[length(x@par):(length(x@par) - x@ncat + 1)]
        d        
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'partcomp'),
    definition = function(x){          
        d <- x@par[(x@nfact+1):(length(x@par)-2)]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'mcm'),
    definition = function(x){          
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
        P <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u, D=x@D)
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE){                  
        a <- x@par[1:x@nfact]        
        d <- x@par[(x@nfact+1):length(x@par)]
        P <- P.poly(a=a, d=d, Theta=Theta, itemexp=itemexp, D=x@D)
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
        P <- P.poly(a=a, d=(d + t), Theta=Theta, itemexp=itemexp, D=x@D)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta){                  
        a <- x@par[1:x@nfact]        
        d <- x@par[-(1:x@nfact)]
        P <- P.gpcm(a=a, d=d, Theta=Theta, D=x@D)
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
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=x@D)
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
        P <- P.comp(a=a, d=d, Theta=Theta, g=g, u=u, D=x@D)
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
        P <- P.mcm(a=a, ak=ak, d=d, t=t, Theta=Theta, D=x@D)
        return(P)
    }
)

##Function passes
P.poly <- function(a, d, Theta, itemexp = FALSE, D)
{    
    ncat <- length(d) + 1
    nfact <- length(a)
    Pk <- matrix(0,nrow(Theta),ncat+1)
    Pk[,1] <- 1	
    for(i in 1:(ncat-1))			
        Pk[ ,i+1] <- P.mirt(a, d[i], Theta, g=0, u=1, D=D)		
    if(itemexp){
        P <- matrix(0,nrow(Theta),ncat)		
        for(i in ncat:1)
            P[ ,i] <- Pk[ ,i] - Pk[ ,i+1]						
        Pk <- P
    }	
    return(Pk)
}

# Trace lines for mirt models
P.mirt <- function(a, d, Theta, g = 0, u = 1, D)
{ 		
    traces <- .Call("traceLinePts", a, d, g, u, Theta, D) ###FIX RCPP CODE
    return(traces)
}

# Trace lines for partially compensetory models
P.comp <- function(a, d, Theta, g, u = 1, D)
{
    nfact <- length(a)
    P <- rep(1,nrow(Theta))
    for(i in 1:nfact)
        P <- P * P.mirt(a[i], d[i], Theta[ ,i, drop=FALSE], g=0, u=1, D=D)
    P <- g + (u - g) * P
    P	
}

#d[1] == 0, ak[1] == 0, ak[length(ak)] == length(ak) - 1 
P.nominal <- function(a, ak, d, Theta, D){
    ncat <- length(d)
    nfact <- ncol(Theta)    
    a <- matrix(a)    
    numerator <- matrix(0, nrow(Theta), ncat)            
    for(i in 1:ncat)
        numerator[ ,i] <- exp(D * ak[i] * (Theta %*% a) + D * d[i])
    P <- numerator/rowSums(numerator)
    return(P)   
}

#d[1] == 0
P.gpcm <- function(a, d, Theta, D){ 
    ncat <- length(d)
    nfact <- ncol(Theta)            
    k <- 0:(ncat-1)
    numerator <- matrix(0, nrow(Theta), ncat)    
    a <- matrix(a)    
    for(i in 1:ncat)
        numerator[ ,i] <- exp(D * k[i] * (Theta %*% a) + D * d[i])
    P <- numerator/rowSums(numerator)
    return(P)   
}

#ak[1] == 0, ak[length(ak)] == length(ak) - 1, d[1] ==0, sum(t) == 1 
P.mcm <- function(a, ak, d, t, Theta, D){
    ncat <- length(d)
    nfact <- ncol(Theta)    
    a <- matrix(a)    
    P <- numerator <- matrix(0, nrow(Theta), ncat)      
    for(i in 1:ncat)
        numerator[ ,i] <- exp(D * ak[i] * (Theta %*% a) + D * d[i])
    denominator <- rowSums(numerator)
    C0 <- 1 / denominator
    t[1] <- 1 - sum(t[2:length(t)])
    T <- matrix(t, nrow(P), ncat, byrow = TRUE)    
    P <- C0 * T + (1 - C0) * numerator/denominator
    return(P)   
}

