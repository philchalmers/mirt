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
setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'custom'),
    definition = function(x){
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'dich'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 rnorm(1L),
                 rnorm(1L, -2, .5),
                 rnorm(1L, 2, .5))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'nestlogit'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 rnorm(1L),
                 rnorm(1L, -2, .5),
                 rnorm(1L, 2, .5), 0, 
                 abs(rnorm(x@ncat-2L, (x@ncat-2L) / 2, sd = 1)), x@ncat,
                 rnorm(x@ncat-2L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'graded'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'rating'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE),
                 rnorm(1L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'gpcm'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=-.2, sdlog=.5), 0,
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'rsm'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5), 0, 
                 sort(rnorm(x@ncat-1L, sd = 2), decreasing=TRUE),
                 rnorm(1L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'nominal'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=-.2, sdlog=.5), 0, 
                 abs(rnorm(x@ncat-1L, (x@ncat-1L) / 2, sd = 1)), 0,
                 rnorm(x@ncat-1L))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'partcomp'),
    definition = function(x){
        par <- c(rlnorm(x@nfact, meanlog=0, sdlog=.5),
                 rlnorm(x@nfact, meanlog=.2, sdlog=.5),
                 rnorm(1L, -2, .5),
                 rnorm(1L, 2, .5))
        x@par[x@est] <- par[x@est]
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'GroupPars'),
    definition = function(x){
        x
    }
)

setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'RandomPars'),
    definition = function(x){
        x
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
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        P <- P.mirt(x@par, Theta=Theta, asMatrix=TRUE, ot=ot)        
        return(P)
    }
)

P.mirt <- function(par, Theta, asMatrix = FALSE, ot = 0)
{
    return(.Call("traceLinePts", par, Theta, asMatrix, ot, FALSE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nestlogit', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nestlogit(x@par, Theta=Theta, correct=x@correctcat, ncat=x@ncat))
    }
)

P.nestlogit <- function(par, Theta, correct, ncat)
{
    return(.Call("nestlogitTraceLinePts", par, Theta, correct, ncat, FALSE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'graded', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, useDesign = TRUE, ot=0){        
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.poly(x@par, Theta=Theta, itemexp=itemexp, ot=ot))
    }
)

P.poly <- function(par, Theta, itemexp = FALSE, ot = 0)
{
    return(.Call('gradedTraceLinePts', par, Theta, itemexp, ot, FALSE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rating', Theta = 'matrix'),
    definition = function(x, Theta, itemexp = TRUE, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.rating(x@par, Theta=Theta, itemexp=itemexp, ot=ot))
    }
)

P.rating <- function(par, Theta, itemexp = FALSE, ot = 0)
{
    return(.Call('gradedTraceLinePts', par, Theta, itemexp, ot, TRUE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'gpcm', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.gpcm(x@par, Theta=Theta, ot=ot))
    }
)

P.gpcm <- function(par, Theta, ot = 0)
{
    return(.Call("gpcmTraceLinePts", par, Theta, ot, FALSE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'rsm', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.rsm(x@par, Theta=Theta, ot=ot))
    }
)

P.rsm <- function(par, Theta, ot = 0)
{
    return(.Call("gpcmTraceLinePts", par, Theta, ot, TRUE))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'nominal', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        a <- x@par[1L:x@nfact]
        ak <- x@par[(x@nfact+1L):(x@nfact + x@ncat)]
        d <- x@par[(length(x@par) - x@ncat + 1L):length(x@par)]
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.nominal(a=a, ak=ak, d=d, Theta=Theta, ot=ot))
    }
)

#d[1] == 0, ak[1] == 0, ak[length(ak)] == length(ak) - 1
P.nominal <- function(a, ak, d, Theta, returnNum = FALSE, ot = 0){
    return(.Call("nominalTraceLinePts", a, ak, d, Theta, returnNum, ot))
}

setMethod(
    f = "ProbTrace",
    signature = signature(x = 'partcomp', Theta = 'matrix'),
    definition = function(x, Theta, useDesign = TRUE, ot=0){
        if(nrow(x@fixed.design) > 1L && useDesign)
            Theta <- cbind(x@fixed.design, Theta)
        return(P.comp(x@par, Theta=Theta, asMatrix=TRUE))
    }
)

P.comp <- function(par, Theta, asMatrix = FALSE, ot = 0)
{
    return(.Call('partcompTraceLinePts', par, Theta, asMatrix, ot, FALSE))
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
            totals <- .Call('denRowSums', fulldata, itemtrace0, itemtrace1, 
                            rep(0, nrow(fulldata)), rep(0, nrow(fulldata)))    
            total_0 <- totals[[1L]] 
            total_1 <- totals[[2L]]            
            total_0 <- tapply(total_0, x@mtch, sum) + log_den0
            total_1 <- tapply(total_1, x@mtch, sum) + log_den1
        } else {
            totals <- .Call('denRowSums', t(fulldata), t(itemtrace0), t(itemtrace1), 
                            rep(0, ncol(fulldata)), rep(0, ncol(fulldata)))
            tmp0 <- totals[[1L]]
            tmp1 <- totals[[2L]]
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
