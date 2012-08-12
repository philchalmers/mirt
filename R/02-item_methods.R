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

