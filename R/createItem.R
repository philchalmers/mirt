#' Create a user defined item with correct generic functions
#' 
#' Initializes the proper S4 class and methods necessary for mirt functions to use in estimation. To use 
#' the defined objects pass to the \code{mirt(..., customItems = list())} command, and ensure that the classes
#' are properly labelled and unique in the list. 
#' 
#' Addionally, the \code{summary()} function will not return proper standardized loadings since the function
#' is not sure how to handle them (no slopes could be defined at all!). Instead loadings of .001 are filled in
#' as placeholders.
#' 
#' @aliases createItem
#' @param name a character indicating the item class name to be defined
#' @param par a named vector of the starting values for the parameters
#' @param est a logical vector indicating which parameters should be freely estimated by default
#' @param P the probability trace function for all categories (first column is category 1, second category two, etc). 
#' First input contains a vector of all the item parameters, while the second input must be a matrix called Theta.
#' Function also must return a \code{matrix} object of category probabilites
#' @param gr gradient function (vector of first derivatives) used in estimation. 
#' If not specified a numeric approximation will be used
#' @param hss hessian function (matrix of second derivatives) used in estimation. 
#' If not specified a numeric approximation will be used (required for the MH-RM algorithm only)
#' @param lbound optional vector indicating the lower bounds of the parameters. If not specified then
#' the bounds will be set to -Inf
#' @param ubound optional vector indicating the lower bounds of the parameters. If not specified then
#' the bounds will be set to Inf
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords createItem
#' @export createItem
#' @examples
#' 
#' \dontrun{
#' 
#' name <- 'old2PL'
#' par <- c(a = .5, b = -2)
#' est <- c(TRUE, TRUE)
#' P.old2PL <- function(par,Theta){
#'      a <- par[1]
#'      b <- par[2] 
#'      P1 <- 1 / (1 + exp(-1.702*a*(Theta - b)))
#'      cbind(1-P1, P1)
#' } 
#' lbound <- c(-Inf, -Inf)
#' ubound <- c(Inf, Inf)
#' 
#' x <- createItem(name, par=par, est=est, lbound=lbound, ubound=ubound, P=P.old2PL)
#' 
#' #So, let's estimate it!
#' dat <- expand.table(LSAT7) 
#' sv <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), pars = 'values')
#' tail(sv) #looks good
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = TRUE)
#' coef(mod)
#' mod2 <- confmirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = TRUE)
#' coef(mod2)
#' 
#' #nonlinear
#' name <- 'nonlin'
#' par <- c(a1 = .5, a2 = .1, d = 0)
#' est <- c(TRUE, TRUE, TRUE)
#' P.nonlin <- function(par,Theta){
#'      a1 <- par[1]
#'      a2 <- par[2] 
#'      d <- par[3]
#'      P1 <- 1 / (1 + exp(-1.702*(a1*Theta + a2*Theta^2 + d)))
#'      cbind(1-P1, P1)
#' } 
#' 
#' x2 <- createItem(name, par=par, est=est, P=P.nonlin)
#' 
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'nonlin'), customItems=list(nonlin=x2), verbose = TRUE)
#' coef(mod)
#' }
createItem <- function(name, par, est, P, gr=NULL, hss = NULL, lbound = NULL, ubound = NULL){    
    setClass(name, contains = 'AllItemsClass', representation = representation(P='function'))    
    setMethod("initialize",
              name,
              function(.Object, par, est, lbound, ubound, P) {                  
                  names(est) <- names(par)
                  .Object@par <- par
                  .Object@est <- est                      
                  .Object@P <- P                      
                  .Object@lbound <- if(!is.null(lbound)) lbound  else rep(-Inf, length(par)) 
                  .Object@ubound <- if(!is.null(ubound)) ubound  else rep(Inf, length(par))                   
                  .Object
              })
    setMethod(
        f = "print",
        signature = signature(x = name),
        definition = function(x, ...){
            cat('Item object of class:', class(x))
        }
    )
    setMethod(
        f = "show",
        signature = signature(object = name),
        definition = function(object){
            print(object)
        }
    )
    setMethod(
        f = "ExtractLambdas",
        signature = signature(x = name),
        definition = function(x){             
            a <- rep(.001, x@nfact)
            a        
        }
    )
    setMethod(
        f = "ProbTrace",
        signature = signature(x = name, Theta = 'matrix'),
        definition = function(x, Theta, fixed.design = NULL){              
            x@P(x@par, Theta=Theta)            
        }
    )
    setMethod(
        f = "LogLik",
        signature = signature(x = name, Theta = 'matrix'),
        definition = function(x, Theta, EM=FALSE, prior=NULL){          
            itemtrace <- ProbTrace(x=x, Theta=Theta)                
            Prior <- rep(1, nrow(itemtrace))
            if(EM) Prior <- prior
            LL <- (-1) * sum(x@rs * log(itemtrace) * Prior)
            LL <- LL.Priors(x=x, LL=LL)        
            return(LL)
        }
    )
    if(is.null(gr) && is.null(hss)){
        setMethod(
            f = "Deriv",
            signature = signature(x = name, Theta = 'matrix'),
            definition = function(x, Theta, EM = FALSE, BFACTOR = FALSE, prior = NULL){
                grad <- rep(0, length(x@par))
                hess <- matrix(0, length(x@par), length(x@par))
                Prior <- rep(1, nrow(x@rs))        
                if(BFACTOR) Prior <- prior
                if(EM){                
                    grad[x@est] <- numDeriv::grad(EML, x@par[x@est], obj=x, Theta=Theta, prior=Prior)
                    #hess[x@est, x@est] <- numDeriv::hessian(EML, x@par[x@est], obj=x, Theta=Theta, prior=Prior)     
                    return(list(grad = grad))            
                }        
                grad[x@est] <- numDeriv::grad(L, x@par[x@est], obj=x, Theta=Theta)
                hess[x@est, x@est] <- numDeriv::hessian(L, x@par[x@est], obj=x, Theta=Theta)
                return(list(grad=grad, hess=hess))
            }
        )       
    }    
    return(new(name, par=par, est=est, lbound=lbound, ubound=ubound, P=P))
}