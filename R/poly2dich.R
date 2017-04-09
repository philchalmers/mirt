#' Change polytomous items to dichotomous item format
#'
#' Tranforms a matrix of items into a new matrix where the select polytomous items have been
#' converted into comperable dichotomous items with the same information.
#'
#' @param data an object of class \code{data.frame} or \code{matrix}
#' @param which.items a vector indicating which items should be transformed into the
#'   dichotomous form. Default uses all input items
#' @return Returns an integer matrix
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
#' @export poly2dich
#' @examples
#'
#' \dontrun{
#' data(Science)
#'
#' head(Science)
#' newScience <- poly2dich(Science)
#' head(newScience)
#'
#' newScience2 <- poly2dich(Science, which.items = 2)
#' head(newScience2)
#'
#' }
#'
poly2dich <- function(data, which.items = 1:ncol(data)) {
    if(missing(data)) missingMsg('data')
    stopifnot(is.data.frame(data) || is.matrix(data))
    stopifnot(length(which.items) > 0L)
    ret <- vector('list', ncol(data))
    nms <- colnames(data)
    for(i in seq_len(ncol(data))){
        if(i %in% which.items){
            old <- data[,i]
            u <- sort(na.omit(unique(old)))
            new <- matrix(0L, nrow(data), length(u))
            colnames(new) <- paste(nms[i], u, sep='.')
            new[is.na(old), ] <- NA
            for(j in seq_len(length(u))) new[u[j] == old, j] <- 1L
            ret[[i]] <- new
        } else {
            ret[[i]] <- data[ ,i, drop=FALSE]
        }

    }
    do.call(cbind, ret)
}

