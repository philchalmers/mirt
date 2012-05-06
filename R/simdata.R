#' Simulate response patterns from an underlying multivariate normal
#' distribution
#' 
#' Simulates response patterns for compensetory and noncompensatory MIRT models
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
#' @param N sample size
#' @param guess a vector of guessing parameters for each item; only applicable
#' for dichotomous items. Must be either a scalar value that will affect all of
#' the dichotomous items, or a vector with as many values as to be simulated
#' items
#' @param sigma a covariance matrix of the underlying distribution. Default is
#' the identity matrix
#' @param mu a mean vector of the underlying distribution. Default is a vector
#' of zeros
#' @param partcomp a logical vector used to specify which items are to be
#' treated as partially compensatory items (see Reckase, 2009), and single
#' values are repeated for each item. Note that when simulating noncompensatory
#' data the slope parameters must equal the number of intercept parameters
#' @param factor.loads logical; are the slope parameters in \code{a} in the
#' factor loadings metric?
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
#' 
#' dataset1 <- simdata(a, d, 2000)
#' dataset2 <- simdata(a, d, 2000, mu = mu, sigma = sigma)
#' 
#' ###An example of a mixed item, bifactor loadings pattern with correlated specific factors
#' # can use factor loadings metric
#' a <- matrix(c(
#' .8,.4,NA,
#' .4,.4,NA,
#' .7,.4,NA,
#' .8,NA,.4,
#' .4,NA,.4,
#' .7,NA,.4),ncol=3,byrow=TRUE)
#' 
#' #first three items are dichotomous, next two have 4 categories, and the last has 3
#' d <- matrix(c(
#' -1.0,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 3.0,2.0,-0.5,
#' 2.5,1.0,-1,
#' 2.0,0.0,NA),ncol=3,byrow=TRUE)
#' 
#' sigma <- diag(3)
#' sigma[2,3] <- sigma[3,2] <- .25
#' 
#' dataset <- simdata(a,d,1000,sigma=sigma,factor.loads=TRUE)
#' 
#' ####Noncompensatory item example
#' a <- matrix(c(
#'   1,NA,
#' 1.5,NA,
#'  NA, 1,
#'  NA,1.6,
#' 1.5,.5,
#'  .7, 1), ncol=2,byrow=TRUE)
#' 
#' #notice that the partially compensatory items have the same number of intercepts as 
#' #factors influencing the item
#' d <- matrix(c(
#' -1.0,NA,
#'  1.5,NA,
#'  0.0,NA,
#'  3.0,NA,
#' 2.5,1.0,
#' 2.0, -1),ncol=2,byrow=TRUE)
#' 
#' compdata <- simdata(a,d,3000, partcomp = c(F,F,F,F,T,T))
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
#' 
#' nonlindata <- simdata(a,d,2000,Theta=Theta)
#'
#'    }
#' 
simdata <- function(a, d, N, sigma = NULL, mu = NULL, guess = 0, 
	partcomp = FALSE, factor.loads = FALSE, Theta = NULL)
{
	dist = 'normal'
	nfact <- ncol(a)
	nitems <- nrow(a)	
	if(factor.loads){
		cs <- sqrt(1 - rowSums(a^2, na.rm = TRUE))
		a <- a / cs
	}
	K <- rep(0,nitems)
	for(i in 1:nitems) K[i] <- sum(!is.na(d[i,]))		
	if(length(partcomp) == 1) partcomp <- rep(partcomp,nitems)
	if(length(partcomp) != nitems) stop("Logical partcomp vector is incorrect")  
	if(length(guess) == 1) guess <- rep(guess,nitems)	
	if(length(guess) != nitems) stop("Guessing parameter is incorrect")  
	guess[K > 1] <- 0
	if(is.null(sigma)) sigma <- diag(nfact)
	if(is.null(mu)) mu <- rep(0,nfact)
	if (!is.null(Theta))
		if(ncol(Theta) != nfact || nrow(Theta) != N) 
			stop("The input Theta matrix does not have the correct dimensions")
	if (dist == 'normal' && is.null(Theta)) Theta <- rmvnorm(N,mu,sigma)     
	data <- matrix(0,N,nitems)
	K[partcomp]	<- 1	
	for(i in 1:nitems){
		if(K[i] == 1){	
			slp <- a[i,!is.na(a[i,])]
			tht <- Theta[,!is.na(a[i,])]
			if(partcomp[i]){										
				P <- P.comp(slp, na.omit(d[i,]), tht, guess[i])				
			} else {
				if(length(slp) == 1) tht <- matrix(tht)			
				P <- P.mirt(slp, na.omit(d[i,1]), tht, guess[i])
			}	
			for (j in 1:N) data[j,i] <- sample(c(0,1), 1, prob = c((1 - P[j]), P[j]))		
		} else {			
			int <- d[i,!is.na(d[i,])]
			slp <- a[i,!is.na(a[i,])]
			tht <- Theta[,!is.na(a[i,])]
			if(length(slp) == 1) tht <- matrix(tht)
			P <- P.poly(slp, int, tht, itemexp = TRUE)	
			for (j in 1:N) data[j,i] <- sample(1:ncol(P), 1, prob = P[j,])				  
		}	  
	}
	colnames(data) <- paste("Item_", 1:nitems, sep="") 
	return(data)
}

