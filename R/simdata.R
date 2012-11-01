#' Simulate response patterns 
#' 
#' Simulates response patterns for compensatory and noncompensatory MIRT models
#' from multivariate normally distributed factor (\eqn{\theta}) scores, or from
#' a user input matrix of \eqn{\theta}'s.
#' 
#' Returns a data matrix simulated from the parameters.
#' 
#' @param a a matrix of slope parameters. If slopes are to be constrained to
#' zero then use \code{NA}. \code{a} may also be a similar matrix specifying
#' factor loadings if \code{factor.loads = TRUE}
#' @param d a matrix of intercepts. The matrix should have as many columns as
#' the item with the largest number of categories, and filled empty locations
#' with \code{NA}
#' @param itemtype a character vector of length \code{nrow(a)} specifying the type of items to simulate. 
#' Can be \code{'dich', 'graded', 'gpcm','nominal'}, or \code{'partcomp'}, for 
#' dichotomous, graded, generalized 
#' partial credit, nominal, and partially compensatory models. Note that 
#' for the gpcm and nominal model there should be as many parameters as desired categories,
#' however to parameterized them for meaningful interpretation the first category intercept should 
#' equal 0 for both models
#' @param nominal a matrix of specific item category slopes for nominal models.
#' Should be the dimensions as the intercept specification with one less column, with \code{NA}
#' in locations where not applicable. Note that during estimation the first slope will be constrained
#' to 0 and the last will be constrained to the number of categories minus 1, 
#' so it is best to set these as the values for the first and last categories as well
#' @param N sample size
#' @param guess a vector of guessing parameters for each item; only applicable
#' for dichotomous items. Must be either a scalar value that will affect all of
#' the dichotomous items, or a vector with as many values as to be simulated
#' items
#' @param upper same as \code{guess}, but for upper bound parameters
#' @param D logistic scaling parameter, default 1.702
#' @param sigma a covariance matrix of the underlying distribution. Default is
#' the identity matrix
#' @param mu a mean vector of the underlying distribution. Default is a vector
#' of zeros
#' @param Theta a user specified matrix of the underlying ability parameters,
#' where \code{nrow(Theta) == N} and \code{ncol(Theta) == ncol(a)}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references Reckase, M. D. (2009). \emph{Multidimensional Item Response
#' Theory}. New York: Springer.
#' @keywords data
#' @export simdata
#' @examples
#' 
#' \dontrun{
#' ###Parameters from Reckase (2009), p. 153
#' a <- matrix(c(
#'  .7471, .0250, .1428,
#'  .4595, .0097, .0692,
#'  .8613, .0067, .4040,
#' 1.0141, .0080, .0470,
#'  .5521, .0204, .1482,
#' 1.3547, .0064, .5362,
#' 1.3761, .0861, .4676,
#'  .8525, .0383, .2574,
#' 1.0113, .0055, .2024,
#'  .9212, .0119, .3044,
#'  .0026, .0119, .8036,
#'  .0008, .1905,1.1945,
#'  .0575, .0853, .7077,
#'  .0182, .3307,2.1414,
#'  .0256, .0478, .8551,
#'  .0246, .1496, .9348,
#'  .0262, .2872,1.3561,
#'  .0038, .2229, .8993,
#'  .0039, .4720, .7318,
#'  .0068, .0949, .6416,
#'  .3073, .9704, .0031,
#'  .1819, .4980, .0020,
#'  .4115,1.1136, .2008,
#'  .1536,1.7251, .0345,
#'  .1530, .6688, .0020,
#'  .2890,1.2419, .0220,
#'  .1341,1.4882, .0050,
#'  .0524, .4754, .0012,
#'  .2139, .4612, .0063,
#'  .1761,1.1200, .0870),30,3,byrow=TRUE)
#' 
#' d <- matrix(c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
#'   .6431,.0122,.0912,.8082,-.1867,.4533,-1.8398,.4139,
#'   -.3004,-.1824,.5125,1.1342,.0230,.6172,-.1955,-.3668,
#'   -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329),ncol=1)
#' 
#' mu <- c(-.4, -.7, .1)
#' sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
#' itemtype <- rep('dich', nrow(a))
#' 
#' dataset1 <- simdata(a, d, 2000, itemtype)
#' dataset2 <- simdata(a, d, 2000, itemtype, mu = mu, sigma = sigma)
#' 
#' ###An example of a mixed item, bifactor loadings pattern with correlated specific factors
#' a <- matrix(c(
#' .8,.4,NA,
#' .4,.4,NA,
#' .7,.4,NA,
#' .8,NA,.4,
#' .4,NA,.4,
#' .7,NA,.4),ncol=3,byrow=TRUE)
#' 
#' d <- matrix(c(
#' -1.0,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 0.0,-1.0,1.5,  #the first 0 here is the recommended constraint for nominal
#' 0.0,1.0,-1, #the first 0 here is the recommended constraint for gpcm 
#' 2.0,0.0,NA),ncol=3,byrow=TRUE)
#' 
#' nominal <- matrix(NA, nrow(d), ncol(d))
#' nominal[4, ] <- c(0,1.2,2) #the first 0 and last (ncat - 1) = 2 values are the recommended constraints
#' 
#' sigma <- diag(3)
#' sigma[2,3] <- sigma[3,2] <- .25
#' items <- c('dich','dich','dich','nominal','gpcm','graded')
#' 
#' dataset <- simdata(a,d,1000,items,sigma=sigma,nominal=nominal)
#'
#' ####Unidimensional nonlinear factor pattern
#' theta <- rnorm(2000)
#' Theta <- cbind(theta,theta^2)
#'
#' a <- matrix(c(
#' .8,.4,
#' .4,.4,
#' .7,.4,
#' .8,NA,
#' .4,NA,
#' .7,NA),ncol=2,byrow=TRUE)
#' d <- matrix(rnorm(6))
#' itemtype <- rep('dich',6)
#' 
#' nonlindata <- simdata(a,d,2000,itemtype,Theta=Theta)
#'
#'    }
#' 
simdata <- function(a, d, N, itemtype, sigma = NULL, mu = NULL, guess = 0, 
	upper = 1, nominal = NULL, Theta = NULL, D = 1.702)
{     
	nfact <- ncol(a)
	nitems <- nrow(a)		
	K <- rep(0,nitems)	
	if(length(guess) == 1) guess <- rep(guess,nitems)	
	if(length(guess) != nitems) stop("Guessing parameter is incorrect")
	if(length(upper) == 1) upper <- rep(upper,nitems)    
	if(length(upper) != nitems) stop("Upper bound parameter is incorrect")
    for(i in 1:length(K)){
        K[i] <- length(na.omit(d[i, ])) + 1
        if(itemtype[i] =='partcomp') K[i] <- 2
        if(any(itemtype[i] == c('gpcm', 'nominal'))) K[i] <- K[i] - 1
    }
    guess[K > 2] <- upper[K > 2] <- NA	
	if(is.null(sigma)) sigma <- diag(nfact)
	if(is.null(mu)) mu <- rep(0,nfact)
	if(!is.null(Theta))
		if(ncol(Theta) != nfact || nrow(Theta) != N) 
			stop("The input Theta matrix does not have the correct dimensions")
	if(is.null(Theta)) Theta <- mvtnorm::rmvnorm(N,mu,sigma)     
    if(is.null(nominal)) nominal <- matrix(NA, nitems, 1)
	data <- matrix(0, N, nitems)	
    a[is.na(a)] <- 0    
	for(i in 1:nitems){
        par <- na.omit(c(a[i, ],nominal[i,],d[i,],guess[i],upper[i]))
        obj <- new(itemtype[i], par=par, nfact=nfact, D=D)
        if(any(itemtype[i] == c('gpcm','nominal'))) 
            obj@ncat <- K[i]
        P <- ProbTrace(obj, Theta)
		for (j in 1:N) 
            data[j,i] <- sample(1:ncol(P), 1, prob = P[j,])        
        if(any(itemtype[i] == c('dich', 'gpcm', 'partcomp'))) data[ ,i] <- data[ ,i] - 1 
	}
	colnames(data) <- paste("Item_", 1:nitems, sep="") 
	return(data)
}

