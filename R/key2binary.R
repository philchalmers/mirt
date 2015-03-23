#' Convert response patterns to binary data
#'
#' The \code{key2binary} function will convert response pattern data to a
#' dichotomous format, given a response key.
#'
#'
#' @param fulldata an object of class \code{data.frame}, \code{matrix}, or
#'   \code{table} with the response patterns
#' @param key a vector consisting of the 'correct' response to the items. Each
#'   value corresponds to each column in \code{fulldata}
#' @return Returns a numeric matrix with all the response patterns in
#'   dichotomous format.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export key2binary
#' @examples
#'
#' \dontrun{
#' data(SAT12)
#' head(SAT12)
#' key <- c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5)
#'
#' dicho.SAT12 <- key2binary(SAT12,key)
#'     }
#'
key2binary <- function (fulldata, key){
    if(missing(fulldata)) missingMsg('fulldata')
    if(missing(key)) missingMsg('key')
    if (ncol(fulldata) != length(key)) stop("Key is not the correct length.\n")
    colname <- colnames(fulldata)
    X <- as.matrix(fulldata)
    colnames(X) <- colname
    X <- t(t(X) == key) + 0
    return(X)
}
