# ----------------------------------------------------------------

# continuous response model (requires original data to be extracted)

setClass("crm", contains = 'AllItemsClass',
         representation = representation(transdat='matrix'))

setMethod(
    f = "print",
    signature = signature(x = 'crm'),
    definition = function(x, ...){
        cat('Item object of class:', class(x))
    }
)

setMethod(
    f = "show",
    signature = signature(object = 'crm'),
    definition = function(object){
        print(object)
    }
)

#extract the slopes (should be a vector of length nfact)
setMethod(
    f = "ExtractLambdas",
    signature = signature(x = 'crm'),
    definition = function(x){
        x@par[seq_len(x@nfact)] #slopes
    }
)

#extract the intercepts
setMethod(
    f = "ExtractZetas",
    signature = signature(x = 'crm'),
    definition = function(x){
        x@par[2L] #intercepts
    }
)

# generating random starting values (only called when, e.g., mirt(..., GenRandomPars = TRUE))
setMethod(
    f = "GenRandomPars",
    signature = signature(x = 'crm'),
    definition = function(x){
        par <- c(rlnorm(1, .2, .2), rnorm(1), 1)
        x@par[x@est] <- par[x@est]
        x
    }
)

# how to set the null model to compute statistics like CFI and TLI (usually just fixing slopes to 0)
setMethod(
    f = "set_null_model",
    signature = signature(x = 'crm'),
    definition = function(x){
        x@par[x@nfact] <- 0
        x@est[x@nfact] <- FALSE
        x
    }
)

# probability trace line function. Must return a matrix with a trace line for each category
setMethod(
    f = "ProbTrace",
    signature = signature(x = 'crm', Theta = 'matrix'),
    definition = function(x, Theta){
        a <- x@par[1L:x@nfact]
        b <- x@par[length(x@par)-1L]
        alpha <- x@par[length(x@par)]
        Z <- x@transdat
        if(nrow(Theta) == 1L || nrow(Theta) == nrow(Z)){
            p <- a/(alpha*sqrt(2*pi))*exp(-(a*(Theta - b - Z/alpha ))^2 / 2)
        } else {
            stop('how did you get here?')
        }
        p <- ifelse(p < 1e-10, 1e-10, p) #numerical constraints to avoid log() problems
        p <- ifelse(p > 1 - 1e-10, 1 - 1e-10, p)
        p
    }
)

# complete-data derivative used in parameter estimation (here it is done numerically)
setMethod(
    f = "Deriv",
    signature = signature(x = 'crm', Theta = 'matrix'),
    definition = function(x, Theta, estHess = FALSE, offterm = numeric(1L)){
        grad <- rep(0, length(x@par))
        hess <- matrix(0, length(x@par), length(x@par))
        if(any(x@est)){
            a <- x@par[1L:x@nfact]
            b <- x@par[length(x@par)-1L]
            alpha <- x@par[length(x@par)]
            Z <- x@transdat
            # grad[1L] <- sum(1/a - a * (Theta - b - Z/alpha)^2)
            # grad[2L] <- sum(a^2 * (Theta - b - Z/alpha) )
            # grad[3L] <- sum(-a^2 * (Theta - b - Z/alpha) * (Z/alpha^2) - 1/alpha )
            grad[x@est] <- numerical_deriv(x@par[x@est], EML, obj=x, Theta=Theta)
            if(estHess){
                hess[x@est, x@est] <- numerical_deriv(x@par[x@est], EML, obj=x,
                                                      Theta=Theta, gradient=FALSE)
            }
        }
        ret <- list(grad=grad, hess=hess)
        # if(x@any.prior) ret <- DerivativePriors(x=x, grad=ret$grad, hess=ret$hess,
        #                                         fill.hess=FALSE)
        ret
    }
)

# derivative of the model wft to the Theta values (done numerically here)
setMethod(
    f = "DerivTheta",
    signature = signature(x = 'crm', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_DerivTheta(x, Theta) #replace with analytical derivatives
    }
)

# derivative of the probability trace line function wrt Theta (done numerically here)
setMethod(
    f = "dP",
    signature = signature(x = 'crm', Theta = 'matrix'),
    definition = function(x, Theta){
        numDeriv_dP(x, Theta) #replace with analytical derivatives
    }
)

