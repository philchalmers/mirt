expand.table <- function(tabdata) {
  if (sum(tabdata[,ncol(tabdata)]) <= nrow(tabdata)) 
    stop("Frequencies must be on the right of the data matrix.")
  itemnames <- colnames(tabdata[,1:(ncol(tabdata) - 1)])
  tabdata <- as.matrix(tabdata)  
  fulldata <- c()  
    for (i in 1:nrow(tabdata)) {
      for (j in 1:tabdata[i,ncol(tabdata)]){ 
        fulldata <- rbind(fulldata, tabdata[i,1:(ncol(tabdata) - 1)])
      }  
    }
  colnames(fulldata) <- itemnames  
  fulldata 
}
