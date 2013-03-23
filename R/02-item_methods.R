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
    signature = signature(x = 'rsm'),
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
    f = "print",
    signature = signature(x = 'GroupPars'),
    definition = function(x, ...){
        cat('Object of class:', class(x))
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
    signature = signature(object = 'rsm'),
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

setMethod(
    f = "show",
    signature = signature(object = 'GroupPars'),
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
        Prior <- rep(1, nrow(itemtrace))
        if(EM) Prior <- prior
        LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
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
    signature = signature(x = 'rsm'),
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
    signature = signature(x = 'rsm'),
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
#Probability Traces
setMethod(
    f = "ProbTrace",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){  
        u <- x@par[length(x@par)]
        g <- x@par[length(x@par)-1]
        d <- x@par[length(x@par)-2]
        a <- x@par[1:x@nfact]   
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u, D=x@D)
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, fixed.design = NULL){                  
        a <- x@par[1:x@nfact]        
        d <- x@par[(x@nfact+1):length(x@par)]
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.poly(a=a, d=d, Theta=Theta, itemexp=itemexp, D=x@D)
        return(P)
    }
)


setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, fixed.design = NULL){
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(length(x@par)-1)]
        t <- x@par[length(x@par)]
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.poly(a=a, d=(d + t), Theta=Theta, itemexp=itemexp, D=x@D)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){                  
        a <- x@par[1:x@nfact]        
        d <- x@par[-(1:x@nfact)]
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.nominal(a=a, ak=0:(length(d)-1), d=d, Theta=Theta, D=x@D)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){                  
        a <- x@par[1:x@nfact]        
        d <- x@par[(x@nfact+1):(length(x@par)-1)]
        t <- x@par[length(x@par)]
        d[-1] <- d[-1] + t
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)        
        P <- P.nominal(a=a, ak=0:(length(d)-1), d=d, Theta=Theta, D=x@D)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){         
        a <- x@par[1:x@nfact]        
        ak <- x@par[(x@nfact+1):(x@nfact + x@ncat)]
        d <- x@par[length(x@par):(length(x@par) - x@ncat + 1)]
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=x@D)
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){    
        nfact <- x@nfact
        a <- x@par[1:nfact]
        d <- x@par[(nfact+1):(length(x@par)-2)]
        g <- x@par[length(x@par)-1]
        u <- x@par[length(x@par)]      
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.comp(a=a, d=d, Theta=Theta, g=g, u=u, D=x@D)
        return(cbind(1.0 - P, P))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){        
        a <- x@par[1:x@nfact]        
        ind <- x@nfact + 1
        ak <- x@par[ind:(ind + x@ncat)]
        ind <- ind + x@ncat + 1
        d <- x@par[ind:(ind + x@ncat)]
        ind <- ind + x@ncat + 1
        t <- x@par[ind:length(x@par)]  
        if(!is.null(fixed.design))
            Theta <- cbind(fixed.design, Theta)
        P <- P.mcm(a=a, ak=ak, d=d, t=t, Theta=Theta, D=x@D)
        return(P)
    }
)

##Function passes
P.poly <- function(a, d, Theta, itemexp = FALSE, D)
{   
    traces <- .Call('gradedTraceLinePts', a, d, Theta, D, itemexp)
    return(traces)        
}

# Trace lines for mirt models
P.mirt <- function(a, d, Theta, g = 0, u = 1, D)
{ 		
    traces <- .Call("traceLinePts", a, d, g, u, Theta, D) 
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
    P[P < 1e-8] <- 1e-8
    P[(1 - P) < 1e-8] <- 1 - 1e-8
    P	
}

#d[1] == 0, ak[1] == 0, ak[length(ak)] == length(ak) - 1 
P.nominal <- function(a, ak, d, Theta, D, returnNum = FALSE){
    ncat <- length(d)
    nfact <- ncol(Theta)    
    a <- matrix(a)    
    numerator <- matrix(0, nrow(Theta), ncat)            
    for(i in 1:ncat)
        numerator[ ,i] <- exp(D * ak[i] * (Theta %*% a) + D * d[i])
    P <- numerator/rowSums(numerator)
    P[P < 1e-8] <- 1e-8
    P[(1 - P) < 1e-8] <- 1 - 1e-8
    if(returnNum) return(numerator)
    return(P)   
}

#ak[1] and d[1] are latent process
P.mcm <- function(a, ak, d, t, Theta, D){
    ncat <- length(t)
    nfact <- ncol(Theta)    
    a <- matrix(a)    
    P <- numerator <- matrix(0, nrow(Theta), ncat)      
    numerator0 <- cbind(numerator, 0)
    for(i in 1:ncat)
        numerator0[ ,i+1] <- numerator[ ,i] <- exp(D * ak[i+1] * (Theta %*% a) + D * d[i+1])
    numerator0[, 1] <- exp(D * ak[1] * (Theta %*% a) + D * d[1])
    denominator <- rowSums(numerator)
    denominator0 <- rowSums(numerator0)
    C0 <- numerator0[,1] / denominator0
    C0 <- matrix(C0, nrow(P), ncat)
    t[1] <- 1 - sum(t[2:length(t)])
    T <- matrix(t, nrow(P), ncat, byrow = TRUE)    
    P <- C0 * T + (1 - C0) * numerator/denominator
    P[P < 1e-8] <- 1e-8
    P[(1 - P) < 1e-8] <- 1 - 1e-8
    return(P)   
}

