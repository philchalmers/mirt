simdata <- function(a, d, N, guess = 0,  sigma = NULL, mu = NULL,
  logconst = TRUE)
{  
  dist = 'normal'
  nfact <- ncol(a)
  nitems <- nrow(a)
  if(length(guess) == 1) guess <- rep(guess,nitems)
  if(length(guess) != nitems) stop("Guessing parameter is incorrect")  
  if(is.null(sigma)) sigma <- diag(nfact)
  if(is.null(mu)) mu <- rep(0,nfact)
  if(logconst) {
    a <- a/1.702
	d <- d/1.702
  }  
  if (dist == 'normal') Theta <- rmvnorm(N,mu,sigma)     
  data <- matrix(0,N,nitems)  
  for( i in 1:nitems){ 
    p <- P.mirt(a[i, ], d[i], Theta, guess[i])	
	for (j in 1:N) data[j,i] <- sample(c(0,1), 1, prob = c((1 - p[j]), p[j]))		
  }  
  colnames(data) <- paste("Item_", 1:nitems, sep="") 
  return(data)
}

