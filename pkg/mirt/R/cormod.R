cormod <- function(fulldata, K, guess, smooth = TRUE) 
{  
	fulldata <- as.matrix(fulldata) 
	nitems <- ncol(fulldata)         
	cormat <- cor(fulldata)      
	if (any(guess > 0)){
		for (i in 1:nitems){
			for (j in 1:nitems){
				if (i < j & K[i] == 2 & K[j] == 2 & guess[i]!= 0 ){         
					g1 <- guess[i]
					g2 <- guess[j]
					tabp <- tab <- table(fulldata[ ,i],fulldata[ ,j])/length(fulldata[ ,i])
					w1 <- (1 - g1)
					w2 <- (1 - g2)
					tabp[1,1] <- tab[1,1]/(w1*w2)
					tabp[1,2] <- (w2*tab[1,2] - g2*tab[1,1])/(w1*w2)
					tabp[2,1] <- (w1*tab[2,1] - g1*tab[1,1])/(w1*w2)
					tabp[2,2] <- 1 - tabp[1,1] - tabp[1,2] - tabp[2,1]
					tabp <- round(tabp*length(fulldata[ ,i]))		  
					cormat[i,j] <- cormat[j,i] <- 
						abs(phi(tabp,6))^(1/1.15)*sign(phi(tabp,6))          		  
				} 
				if(i < j & K[i] == 2 & K[j] > 2 & guess[i]!= 0) 
					cormat[i,j] <- cormat[j,i] <- abs(cormat[i,j])^(1/1.15) * sign(cormat[i,j])
			}	
		}      
	} 
	cormat <- abs(cormat)^(1/1.15) * sign(cormat)  
	if(smooth){  
		eig <- eigen(cormat)
		negvalues <- eig$values < 0 
		if (any(negvalues)) {
			negeig <- sum(abs(eig$value[eig$value < 0])) 
			eig$value[eig$value < 0] <- 0
			L <- nitems/sum(eig$value)*eig$value[!negvalues]
			V <- eig$vector[ ,!negvalues] 
			cormat <- V %*% diag(L) %*% t(V)    
		}      
	}	
	cormat
}  

