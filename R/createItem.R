#' Create a user defined item with correct generic functions
#'
#' Initializes the proper S4 class and methods necessary for mirt functions to use in estimation. 
#' To use the defined objects pass to the \code{mirt(..., customItems = list())} command, and 
#' ensure that the classes are properly labeled and unique in the list.
#'
#' The \code{summary()} function will not return proper standardized loadings since the function
#' is not sure how to handle them (no slopes could be defined at all!). Instead loadings of .001 
#' are filled in as place-holders. 
#'
#' @aliases createItem
#' @param name a character indicating the item class name to be defined
#' @param par a named vector of the starting values for the parameters
#' @param est a logical vector indicating which parameters should be freely estimated by default
#' @param P the probability trace function for all categories (first column is category 1, second 
#'   category two, etc). First input contains a vector of all the item parameters, the second input
#'   must be a matrix called \code{Theta}, and the third input must be the number of categories 
#'   called \code{ncat}. Function also must return a \code{matrix} object of category probabilities
#' @param gr gradient function (vector of first derivatives) of the log-likelihood used in 
#'   estimation. The function must be of the form \code{gr(x, Theta)}, where \code{x} is the object 
#'   defined by \code{createItem()} and \code{Theta} is a matrix of latent trait parameters. 
#'   Tabulated (EM) or raw (MHRM) data are located in the \code{x@@dat} slot, and are used to form 
#'   the complete data log-likelihood. If not specified a numeric approximation will be used
#' @param hss Hessian function (matrix of second derivatives) of the log-likelihood used in 
#'   estimation. If not specified a numeric approximation will be used (required for the MH-RM 
#'   algorithm only). The input is idential to the \code{gr} argument
#' @param lbound optional vector indicating the lower bounds of the parameters. If not specified 
#'   then the bounds will be set to -Inf
#' @param ubound optional vector indicating the lower bounds of the parameters. If not specified 
#'   then the bounds will be set to Inf
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
#' P.old2PL <- function(par,Theta, ncat){
#'      a <- par[1]
#'      b <- par[2]
#'      P1 <- 1 / (1 + exp(-1*a*(Theta - b)))
#'      cbind(1-P1, P1)
#' }
#'
#' x <- createItem(name, par=par, est=est, P=P.old2PL)
#'
#' #So, let's estimate it!
#' dat <- expand.table(LSAT7)
#' sv <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), pars = 'values')
#' tail(sv) #looks good
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x))
#' coef(mod)
#' mod2 <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), method = 'MHRM')
#' coef(mod2)
#' 
#' ###non-linear
#' name <- 'nonlin'
#' par <- c(a1 = .5, a2 = .1, d = 0)
#' est <- c(TRUE, TRUE, TRUE)
#' P.nonlin <- function(par,Theta, ncat=2){
#'      a1 <- par[1]
#'      a2 <- par[2]
#'      d <- par[3]
#'      P1 <- 1 / (1 + exp(-1*(a1*Theta + a2*Theta^2 + d)))
#'      cbind(1-P1, P1)
#' }
#'
#' x2 <- createItem(name, par=par, est=est, P=P.nonlin)
#'
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'nonlin'), customItems=list(nonlin=x2))
#' coef(mod)
#'
#' ###nominal response model (Bock 1972 version)
#' Tnom.dev <- function(ncat) {
#'    T <- matrix(1/ncat, ncat, ncat - 1)
#'    diag(T[-1, ]) <-  diag(T[-1, ]) - 1
#'    return(T)
#' }
#'
#' name <- 'nom'
#' par <- c(alp=c(3,0,-3),gam=rep(.4,3))
#' est <- rep(TRUE, length(par))
#' P.nom <- function(par, Theta, ncat){
#'    alp <- par[1:(ncat-1)]
#'    gam <- par[ncat:length(par)]
#'    a <- Tnom.dev(ncat) %*% alp
#'    c <- Tnom.dev(ncat) %*% gam
#'    z <- matrix(0, nrow(Theta), ncat)
#'    for(i in 1:ncat)
#'        z[,i] <- a[i] * Theta + c[i]
#'    P <- exp(z) / rowSums(exp(z))
#'    P
#' }
#'
#' nom1 <- createItem(name, par=par, est=est, P=P.nom)
#' nommod <- mirt(Science, 1, 'nom1', customItems=list(nom1=nom1))
#' coef(nommod)
#' Tnom.dev(4) %*% coef(nommod)[[1]][1:3] #a
#' Tnom.dev(4) %*% coef(nommod)[[1]][4:6] #d
#'
#' }
createItem <- function(name, par, est, P, gr=NULL, hss = NULL, lbound = NULL, ubound = NULL){
    if(any(names(par) %in% c('g', 'u')) || any(names(est) %in% c('g', 'u')))
        stop('Parameter names cannot be \'g\' or \'u\', please change.')
    return(new('custom', name=name, par=par, est=est, lbound=lbound,
               ubound=ubound, P=P, gr=gr, hss=hss, userdata=NULL))
}
