# Methods
setMethod(
    f = "print",
    signature = signature(x = 'ConfirmatoryClass'),
    definition = function(x)
    {
        class(x) <- 'ExploratoryClass'
        print(x)
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object) {
        print(object)
    }
)

setMethod(
    f = "summary",
    signature = 'ConfirmatoryClass',
    definition = function(object, suppress = 0, digits = 3, verbose = TRUE, ...)
    {
        nfact <- ncol(object@F)
        itemnames <- colnames(object@data)
        F <- object@F
        h2 <- as.matrix(object@h2)
        colnames(h2) <- 'h2'
        rownames(F) <- itemnames
        SS <- apply(F^2,2,sum)
        gpars <- ExtractGroupPars(object@pars[[length(object@pars)]])
        Phi <- gpars$gcov
        Phi <- round(Phi, digits)
        colnames(Phi) <- rownames(Phi) <- paste('F',1:ncol(Phi), sep='')
        if(verbose){
            cat("\nFactor loadings metric: \n")
            print(cbind(F, h2),digits)
            cat("\nSS loadings: ",round(SS,digits), "\n")
            cat("\nFactor covariance: \n")
            print(Phi)
        }
        invisible(list(F=F, fcor=Phi))
    }
)

setMethod(
    f = "coef",
    signature = 'ConfirmatoryClass',
    definition = function(object, ...)
    {
        class(object) <- 'ExploratoryClass'
        coef(object,  ...)
    }
)

setMethod(
    f = "residuals",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, ...)
    {
        class(object) <- 'ExploratoryClass'
        residuals(object, ...)
    }
)

setMethod(
    f = "anova",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, object2, ...)
    {
        class(object) <- 'ExploratoryClass'
        anova(object, object2, ...)
    }
)

setMethod(
    f = "fitted",
    signature = signature(object = 'ConfirmatoryClass'),
    definition = function(object, ...){
        class(object) <- 'ExploratoryClass'
        fitted(object, ...)
    }
)

setMethod(
    f = "plot",
    signature = signature(x = 'ConfirmatoryClass', y = 'missing'),
    definition = function(x, y, ...)
    {
        class(x) <- 'ExploratoryClass'
        plot(x, ...)

    }
)

mirt2traditional <- function(x){
    cls <- class(x)
    par <- x@par
    if(cls != 'GroupPars')
        ncat <- x@ncat
    if(cls == 'dich'){
        par[2] <- -par[2]/par[1]
        names(par) <- c('a', 'b', 'g', 'u')
    } else if(cls == 'graded'){
        for(i in 2:ncat)
            par[i] <- -par[i]/par[1]
        names(par) <- c('a', paste0('b', 1:(length(par)-1)))
    } else if(cls == 'gpcm'){
        ds <- par[-1]/par[1]
        newd <- numeric(length(ds)-1)
        for(i in 1:length(newd))
            newd[i] <- -(ds[i+1] - ds[i])
        par <- c(par[1], newd)
        names(par) <- c('a', paste0('b', 1:length(newd)))
    } else if(cls == 'nominal'){
        as <- par[2:(ncat+1)] * par[1]
        as <- as - mean(as)
        ds <- par[(ncat+2):length(par)]
        ds <- ds - mean(ds)
        par <- c(as, ds)
        names(par) <- c(paste0('a', 1:ncat), paste0('c', 1:ncat))
    } else if(cls == 'nestlogit'){
        par1 <- par[1:4]
        par1[2] <- -par1[2]/par1[1]
        names(par1) <- c('a', 'b', 'g', 'u')
        par2 <- par[5:length(par)]
        as <- par2[1:(ncat-1)]
        as <- as - mean(as)
        ds <- par2[-c(1:(ncat-1))]
        ds <- ds - mean(ds)
        names(as) <- paste0('a', 1:(ncat-1))
        names(ds) <- paste0('c', 1:(ncat-1))
        par <- c(par1, as, ds)
    } else {
        names(par) <- names(x@est)
    }
    ret <- matrix(par, 1L, dimnames=list('par', names(par)))
    ret
}
