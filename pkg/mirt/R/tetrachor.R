tetrachor <- function(fulldata, guess = 0, smooth = TRUE) 
{  
  fulldata <- as.matrix(fulldata)
  if (!any(fulldata %in% c(0,1,NA))) stop("Data must contain only 0, 1, or NA.")
  nitems <- ncol(fulldata)
  if (length(guess) == 1) guess <- rep(guess,nitems)
    else if (length(guess) > nitems || length(guess) < nitems) 
      stop("The number of guessing parameters is incorrect.")  
  tmat <- diag(nitems)  
  fulldata[is.na(fulldata)] <- 0  
  cormat <- cor(fulldata)  
  datavec <- rep(0,(nitems*(nitems-1) /2 * 4))
  k = 0  
  if (sum(guess) > 0){
    for (i in 1:nitems){
     for (j in 1:nitems){
       if (i < j & !is.na(cormat[i,j])) 
       {         
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
          datavec[k*4 + 1:4] <- as.integer(tabp)
          k <- k + 1  		  
        }
      }
    }      
  } else {  
    for (i in 1:nitems){       
      for (j in 1:nitems){         
        if (i < j & !is.na(cormat[i,j])){ 
		  tabp <- table(fulldata[,i],fulldata[,j])
		  datavec[k*4 + 1:4] <- as.integer(tabp)
          k <- k + 1
		}
      }
	}		
  }  
  polycors <- .Call("tetrachoric", as.integer(datavec))
  k <- 1
  for(i in 1:nitems){
    for(j in 1:nitems){
	  if(i < j & !is.na(tmat[i,j])){
	    tmat[i,j] <- tmat[j,i] <- polycors[k]
		k <- k + 1
      }
	}
  }  
  if(smooth){  
    eig <- eigen(tmat)
    negvalues <- eig$values < 0 
    if (any(negvalues)) {
      negeig <- sum(abs(eig$value[eig$value < 0])) 
      eig$value[eig$value < 0] <- 0
      L <- nitems/sum(eig$value)*eig$value[!negvalues]
      V <- eig$vector[ ,!negvalues] 
      tmat <- V %*% diag(L) %*% t(V)    
    }      
  }	
  tmat  
}

