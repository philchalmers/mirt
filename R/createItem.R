#' Create a user defined item with correct generic functions
#'
#' Initializes the proper S4 class and methods necessary for mirt functions to use in estimation. To use
#' the defined objects pass to the \code{mirt(..., customItems = list())} command, and ensure that the classes
#' are properly labelled and unique in the list.
#'
#' Additionally, the \code{summary()} function will not return proper standardized loadings since the function
#' is not sure how to handle them (no slopes could be defined at all!). Instead loadings of .001 are filled in
#' as place-holders.
#'
#' @aliases createItem
#' @param name a character indicating the item class name to be defined
#' @param par a named vector of the starting values for the parameters
#' @param est a logical vector indicating which parameters should be freely estimated by default
#' @param P the probability trace function for all categories (first column is category 1, second category two, etc).
#' First input contains a vector of all the item parameters, the second input must be a matrix called \code{Theta}, and
#' the third input must be the number of categories called \code{ncat}.
#' Function also must return a \code{matrix} object of category probabilities
#' @param gr gradient function (vector of first derivatives) used in estimation.
#' If not specified a numeric approximation will be used
#' @param hss Hessian function (matrix of second derivatives) used in estimation.
#' If not specified a numeric approximation will be used (required for the MH-RM algorithm only)
#' @param userdata an optional matrix of person level covariate data that can be used in estimation. This
#' matrix with be used in the probability function by passing \code{Theta = cbind(Theta, userdata)}. Note that
#' this only makes sense to use when the estimation uses the MH-RM engine since the number of rows in Theta
#' will be the same as the number of rows in the covariate data (similar to how \code{mixedmirt} works)
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
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), verbose = TRUE)
#' coef(mod)
#' mod2 <- mirt(dat, 1, c(rep('2PL',4), 'old2PL'), customItems=list(old2PL=x), 
#'    verbose = TRUE, method = 'MHRM')
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
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'nonlin'), customItems=list(nonlin=x2), verbose = TRUE)
#' coef(mod)
#'
#' ### covariate included data
#' name <- 'mycov'
#' par <- c(a1 = .5, a2 =.5, d = 0)
#' est <- c(TRUE, TRUE, TRUE)
#' P.mycov <- function(par,Theta, ncat){
#'      a1 <- par[1]
#'      a2 <- par[2]
#'      d <- par[3]
#'      #notice here that the covariate data is found in Theta,
#'      #    use browser() to jump in for debugging if needed
#'      P1 <- 1 / (1 + exp(-1*(a1 * Theta[,1] + a2*Theta[,2] + d)))
#'      cbind(1-P1, P1)
#' }
#'
#' covdata <- matrix(c(rep(0, 500), rep(1,500)), nrow=nrow(dat))
#' x3 <- createItem(name, par=par, est=est, P=P.mycov, userdata=covdata)
#' mod <- mirt(dat, 1, c(rep('2PL',4), 'mycov'), customItems=list(mycov=x3), method = 'MHRM')
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
createItem <- function(name, par, est, P, gr=NULL, hss = NULL, lbound = NULL, ubound = NULL, userdata = NULL){
    return(new('custom', name=name, par=par, est=est, lbound=lbound,
               ubound=ubound, P=P, gr=gr, hss=hss, userdata=userdata))
}
