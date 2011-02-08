tetrachor <- function(fulldata, guess = 0, nowarn = TRUE) 
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
  if(nowarn) options(warn = -1)
  if (sum(guess) > 0){
    for (i in 1:nitems){
     for (j in 1:nitems){
       if (i < j & !is.na(cormat[i,j])) 
       {          
          d1 <- fulldata[ ,i]
          d2 <- fulldata[ ,j]          
          g1 <- guess[i]
          g2 <- guess[j]
          tabp <- tab <- table(d1,d2)/length(d1)
          w1 <- (1 - g1)
          w2 <- (1 - g2)
          tabp[1,1] <- tab[1,1]/(w1*w2)
          tabp[1,2] <- (w2*tab[1,2] - g2*tab[1,1])/(w1*w2)
          tabp[2,1] <- (w1*tab[2,1] - g1*tab[1,1])/(w1*w2)
          tabp[2,2] <- 1 - tabp[1,1] - tabp[1,2] - tabp[2,1]
          tabp <- tabp*length(d1)
          tmat[i,j] <- polychor(tabp)            
        }
      }
    }  
    tmat <- tmat + t(tmat) - diag(nitems)
  } else {  
    for (i in 1:nitems)       
      for (j in 1:nitems)         
        if (i < j & !is.na(cormat[i,j])) tmat[i,j] <- polychor(table(fulldata[,i],fulldata[,j]))
        tmat <- tmat + t(tmat) - diag(nitems)
  } 
  
  eig <- eigen(tmat)
  negvalues <- eig$values < 0 
  if (any(negvalues)) {
    negeig <- sum(abs(eig$value[eig$value < 0])) 
    eig$value[eig$value < 0] <- 0
    L <- nitems/sum(eig$value)*eig$value[!negvalues]
    V <- eig$vector[ ,!negvalues] 
    tmat <- V %*% diag(L) %*% t(V)    
  }    
  if(nowarn) options(warn = 0)
  tmat  
}
