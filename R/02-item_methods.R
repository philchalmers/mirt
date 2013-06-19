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
    signature = signature(x = 'nestlogit'),
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
    signature = signature(object = 'nestlogit'),
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
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
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
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'nestlogit'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'graded'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'rating'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'rsm'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'nominal'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        x@par[1L:x@nfact]
    }
)

#----------------------------------------------------------------------------
setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'dich'),
    definition = function(x){
        par <- x@par
        d <- par[1L:x@nfact]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'graded'),
    definition = function(x){
        par <- x@par
        d <- par[-(1L:x@nfact)]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'rating'),
    definition = function(x){
        par <- x@par
        d <- par[-c(1L:x@nfact, length(par))]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        par <- x@par
        d <- par[-(1L:x@nfact)]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'rsm'),
    definition = function(x){
        par <- x@par
        d <- par[-(1L:x@nfact)]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'nominal'),
    definition = function(x){
        d <- x@par[(length(x@par) - x@ncat + 1L):length(x@par)]
        d
    }
)

setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        d <- x@par[(x@nfact+1L):(length(x@par)-2L)]
        d
    }
)

#----------------------------------------------------------------------------
#Probability Traces
setMethod(
    f = "ProbTrace",
    signature = signature(x = 'custom', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(x@useuserdata) Theta <- cbind(Theta, x@userdata)
        return(x@P(x@par, Theta=Theta, ncat=x@ncat))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'dich', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        u <- x@par[length(x@par)]
        g <- x@par[length(x@par)-1L]
        d <- x@par[length(x@par)-2L]
        a <- x@par[1L:x@nfact]        
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        P <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u, D=x@D, asMatrix=TRUE, ot=ot)        
        return(P)
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        a <- x@par[1L:x@nfact]
        d <- x@par[x@nfact+1L]
        g <- x@par[x@nfact+2L]
        u <- x@par[x@nfact+3L]
        ak <- x@par[(x@nfact+4L):(x@nfact+4L+x@ncat-2L)]
        dk <- x@par[(length(x@par)-length(ak)+1):length(x@par)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nestlogit(a=a, d=d, Theta=Theta, g=g, u=u,
                         ak=ak, dk=dk, correct=x@correctcat, D=x@D))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, useDesign = TRUE, ot=0){
        a <- x@par[1L:x@nfact]
        d <- x@par[(x@nfact+1L):length(x@par)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.poly(a=a, d=d, Theta=Theta, itemexp=itemexp, D=x@D))
    }
)


setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, useDesign = TRUE, ot=0){
        nfact <- x@nfact
        a <- x@par[1L:nfact]
        d <- x@par[(nfact+1L):(length(x@par)-1L)]
        t <- x@par[length(x@par)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.poly(a=a, d=(d + t), Theta=Theta, itemexp=itemexp, D=x@D, ot=ot))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        a <- x@par[1L:x@nfact]
        d <- x@par[-(1L:x@nfact)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nominal(a=a, ak=0:(length(d)-1), d=d, Theta=Theta, D=x@D, ot=ot))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        a <- x@par[1L:x@nfact]
        d <- x@par[(x@nfact+1L):(length(x@par)-1L)]
        t <- x@par[length(x@par)]
        d[-1L] <- d[-1L] + t
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nominal(a=a, ak=0:(length(d)-1), d=d, Theta=Theta, D=x@D))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        a <- x@par[1L:x@nfact]
        ak <- x@par[(x@nfact+1L):(x@nfact + x@ncat)]
        d <- x@par[(length(x@par) - x@ncat + 1L):length(x@par)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nominal(a=a, ak=ak, d=d, Theta=Theta, D=x@D))
    }
)

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        nfact <- x@nfact
        a <- x@par[1L:nfact]
        d <- x@par[(nfact+1L):(length(x@par)-2L)]
        g <- x@par[length(x@par)-1L]
        u <- x@par[length(x@par)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.comp(a=a, d=d, Theta=Theta, g=g, u=u, D=x@D, asMatrix=TRUE))
    }
)

##Function passes
P.poly <- function(a, d, Theta, itemexp = FALSE, D)
{
    return(.Call('gradedTraceLinePts', a, d, Theta, D, itemexp))
}

# Trace lines for mirt models
P.mirt <- function(a, d, Theta, g = 0, u = 1, D, asMatrix = FALSE, ot = 0)
{
    return(.Call("traceLinePts", a, d, g, u, Theta, D, asMatrix, ot))
}

# Trace lines for partially compensetory models
P.comp <- function(a, d, Theta, g, u = 1, D, asMatrix = FALSE)
{
    nfact <- length(a)
    P <- rep(1,nrow(Theta))
    for(i in 1L:nfact)
        P <- P * P.mirt(a[i], d[i], Theta[ ,i, drop=FALSE], g=0, u=1, D=D)
    P <- g + (u - g) * P
    s.eps <- 1e-10
    P[P < s.eps] <- s.eps
    P[(1 - P) < s.eps] <- 1 - s.eps
    if(asMatrix) return(cbind(1-P, P))
    else return(P)
}

#d[1] == 0, ak[1] == 0, ak[length(ak)] == length(ak) - 1
P.nominal <- function(a, ak, d, Theta, D, returnNum = FALSE, ot = 0){
    return(.Call("nominalTraceLinePts", a, ak, d, Theta, D, returnNum, ot))
}

P.nestlogit <- function(a, d, Theta, g, u, ak, dk, correct, D)
{
    traces <- matrix(0, nrow(Theta), length(ak)+1L)
    traces[ ,correct] <- P.mirt(a=a, d=d, Theta=Theta, g=g, u=u, D=D)
    Q <- 1 - traces[ ,correct]
    Pn <- P.nominal(a=rep(1,ncol(Theta)), ak=ak, d=dk, Theta=Theta, D=D)
    traces[ ,-correct] <- Q * Pn
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


#----------------------------------------------------------------------------
## Random variable methods

setMethod(
    f = "DrawValues",
    signature = signature(x = 'RandomPars', Theta = 'matrix'),
    definition = function(x, Theta, pars, fulldata, itemloc, offterm0){
        tol <- .Machine$double.eps
        J <- length(pars) - 1L
        theta0 <- x@drawvals
        N <- nrow(theta0)
        unif <- runif(N)
        prior.mu <- rep(0, ncol(theta0))
        prior.t.var <- matrix(0, ncol(theta0), ncol(theta0))
        prior.t.var[lower.tri(prior.t.var, diag=TRUE)] <- x@par
        d <- if(ncol(theta0) == 1) matrix(prior.t.var) else diag(diag(prior.t.var))
        prior.t.var <- prior.t.var + t(prior.t.var) - d
        sigma <- if(ncol(theta0) == 1L) matrix(x@cand.t.var) else diag(rep(x@cand.t.var,ncol(theta0)))
        theta1 <- theta0 + mvtnorm::rmvnorm(N, prior.mu, sigma)
        log_den0 <- mvtnorm::dmvnorm(theta0,prior.mu,prior.t.var,log=TRUE)
        log_den1 <- mvtnorm::dmvnorm(theta1,prior.mu,prior.t.var,log=TRUE)
        itemtrace0 <- itemtrace1 <- matrix(0, ncol=ncol(fulldata), nrow=nrow(fulldata))
        if(x@between){
            offterm1 <- matrix(0, nrow(itemtrace0), J)            
            tmp1 <- rowSums(x@gdesign * theta1[x@mtch, , drop=FALSE])
            for(i in 1L:J) offterm1[,i] <- tmp1            
        } else {            
            tmp1 <- rowSums(x@gdesign * theta1[x@mtch, , drop=FALSE])            
            offterm1 <- matrix(tmp1, nrow(itemtrace0), J, byrow = TRUE)            
        }
        offterm1 <- offterm1 + offterm0
        for (i in 1L:J){
            itemtrace0[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <-
                ProbTrace(x=pars[[i]], Theta=Theta, ot=offterm0[,i])
            itemtrace1[ ,itemloc[i]:(itemloc[i+1L] - 1L)] <-
                ProbTrace(x=pars[[i]], Theta=Theta, ot=offterm1[,i])
        }
        if(x@between){
            total_0 <- rowSums(log(itemtrace0)*fulldata)
            total_1 <- rowSums(log(itemtrace1)*fulldata)
            total_0 <- tapply(total_0, x@mtch, sum) + log_den0
            total_1 <- tapply(total_1, x@mtch, sum) + log_den1
        } else {
            tmp0 <- colSums(log(itemtrace0)*fulldata)
            tmp1 <- colSums(log(itemtrace1)*fulldata)
            LL0 <- LL1 <- numeric(J)
            for(i in 1L:J){
                LL0[i] <- sum(tmp0[itemloc[i]:(itemloc[i+1L] - 1L)])
                LL1[i] <- sum(tmp1[itemloc[i]:(itemloc[i+1L] - 1L)])
            }
            total_0 <- LL0[x@mtch] + log_den0
            total_1 <- LL1[x@mtch] + log_den1
        }
        diff <- total_1 - total_0
        accept <- diff > 0
        accept[unif < exp(diff)] <- TRUE
        theta1[!accept, ] <- theta0[!accept, ]
        total_1[!accept] <- total_0[!accept]
        attr(theta1, "Proportion Accepted") <- sum(accept)/N
        return(theta1)
    }
)

setMethod(
    f = "RandomDeriv",
    signature = signature(x = 'RandomPars'),
    definition = function(x){
        browser()
        
    }
)
