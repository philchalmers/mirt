#' Score a test by converting response patterns to binary data
#'
#' The \code{key2binary} function will convert response pattern data to a
#' dichotomous format, given a response key.
#'
#' @param fulldata an object of class \code{data.frame}, \code{matrix}, or
#'   \code{table} with the response patterns
#' @param key a vector or matrix consisting of the 'correct' response to the items. Each
#'   value/row corresponds to each column in \code{fulldata}. If the input is a matrix, multiple
#'   scoring keys can be supplied for each item. NA values are used to indicate no scoring key (or
#'   in the case of a matrix input, no additional scoring keys)
#' @return Returns a numeric matrix with all the response patterns in
#'   dichotomous format
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export key2binary
#' @examples
#'
#' data(SAT12)
#' head(SAT12)
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
#'
#' dicho.SAT12 <- key2binary(SAT12, key)
#' head(dicho.SAT12)
#'
#' # multiple scoring keys
#' key2 <- cbind(c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5),
#'               c(2,3,NA,1,rep(NA, 28)))
#' dicho.SAT12 <- key2binary(SAT12, key2)
#'
#' # keys from raw character responses
#' resp <- as.data.frame(matrix(c(
#'   "B","B","D","D","E",
#'   "B","A","D","D","E",
#'   "B","A","D","C","E",
#'   "D","D","D","C","E",
#'   "B","C","A","D","A"), ncol=5, byrow=TRUE))
#'
#' key <- c("B", "D", "D", "C", "E")
#'
#' d01 <- key2binary(resp, key)
#' head(d01)
#'
key2binary <- function (fulldata, key){
    if(missing(fulldata)) missingMsg('fulldata')
    if(missing(key)) missingMsg('key')
    if(is.vector(key)) key <- matrix(key)
    if (ncol(fulldata) != nrow(key)) stop("Key is not the correct length.\n", call.=FALSE)
    colname <- colnames(fulldata)
    X <- matrix(0L, nrow(fulldata), ncol(fulldata))
    colnames(X) <- colname
    for(i in 1L:ncol(X)){
        if(all(is.na(key[i,]))) next
        X[,i] <- fulldata[,i] %in% key[i,] + 0L
    }
    return(X)
}
