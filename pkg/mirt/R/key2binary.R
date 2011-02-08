key2binary <- function (fulldata, key){    
  if (ncol(fulldata) != length(key)) stop("Key is not the correct length.\n") 
  colname <- colnames(fulldata)	
  X <- as.matrix(fulldata)
  colnames(X) <- colname                
  X <- t(t(X) == key) + 0   
  return(X)
}
