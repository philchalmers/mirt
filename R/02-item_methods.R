#----------------------------------------------------------------------------
setMethod(
    f = "print",
    signature = signature(x = 'custom'),
    definition = function(x, ...){
        cat('Custom item object named:', x@name)
    }
)

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
    signature = signature(object = 'custom'),
    definition = function(object){
        print(object)
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
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta, EM=FALSE, prior=NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)                
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)        
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, EM=FALSE, prior=NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)                        
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)        
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)        
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)        
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)
    }
)

setMethod(
    f = "LogLik",
    signature = signature(x = 'mcm', Theta = 'matrix'),
    definition = function(x, Theta, EM = FALSE, prior = NULL){          
        itemtrace <- ProbTrace(x=x, Theta=Theta)
        if(EM) LL <- (-1) * sum(x@rs * log(itemtrace) * prior)
        else LL <- (-1) * sum(x@rs * log(itemtrace))
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
        rlist <- Estep.mirt(pars=pars, tabdata=tabdata, Theta=Theta, prior=prior, itemloc=itemloc)
        LL <- (-1)*sum(r*log(rlist$expected)) 
        LL <- LL.Priors(x=x, LL=LL)
        return(LL)        
    }
)

#----------------------------------------------------------------------------

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'custom'),
    definition = function(x){             
        a <- rep(.001, x@nfact)
        a        
    }
)

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
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta, fixed.design = NULL){          
        if(x@useuserdata) Theta <- cbind(Theta, x@userdata)
        x@P(x@par, Theta=Theta, ncat=x@ncat)            
    }
)

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
    s.eps <- 1e-10
    P[P < s.eps] <- s.eps
    P[(1 - P) < s.eps] <- 1 - s.eps
    P	
}

#d[1] == 0, ak[1] == 0, ak[length(ak)] == length(ak) - 1 
P.nominal <- function(a, ak, d, Theta, D, returnNum = FALSE){    
    traces <- .Call("nominalTraceLinePts", a, ak, d, Theta, D, returnNum) 
    return(traces)    
}

#ak[1] and d[1] are latent process
P.mcm <- function(a, ak, d, t, Theta, D){    
    t[1] <- 1 - sum(t[2:length(t)])
    traces <- .Call("mcmTraceLinePts", a, ak, d, t, Theta, D) 
    return(traces)        
}

#----------------------------------------------------------------------------
## initialize for custom items
setMethod("initialize",
          'custom',
          function(.Object, name, par, est, lbound, ubound, P, gr, hss, userdata) {                  
              dummyfun <- function(...) return(NULL)
              names(est) <- names(par)
              usegr <- usehss <- useuserdata <- TRUE
              .Object@name <- name
              .Object@par <- par
              .Object@est <- est                      
              .Object@P <- P
              if(is.null(gr)){
                  .Object@gr <- dummyfun
                  usegr <- FALSE                  
              } else .Object@gr <- gr
              if(is.null(hss)){
                  .Object@hss <- dummyfun
                  usehss <- FALSE                  
              } else .Object@hss <- hss
              if(is.null(userdata)){
                  .Object@userdata <- matrix(NaN)
                  useuserdata <- FALSE                  
              } else .Object@userdata <- userdata
              .Object@usegr <- usegr
              .Object@usehss <- usehss
              .Object@useuserdata <- useuserdata              
              .Object@lbound <- if(!is.null(lbound)) lbound  else rep(-Inf, length(par)) 
              .Object@ubound <- if(!is.null(ubound)) ubound  else rep(Inf, length(par))                   
              .Object
          })
