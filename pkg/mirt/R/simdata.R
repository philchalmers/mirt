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
