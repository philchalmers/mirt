#' Expand summary table of patterns and frequencies
#'
#' The \code{expand.table} function expands a summary table of unique response
#' patterns to a full sized data-set. The response frequencies must be on the
#' rightmost column of the input data.
#'
#' @param tabdata An object of class \code{data.frame} or \code{matrix}
#'   with the unique response patterns and the number of frequencies
#'   in the rightmost column
#' @param sample logical; randomly switch the rows in the expanded table? This does not change the
#'   expanded data, only the row locations
#' @return Returns a numeric matrix with all the response patterns.
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords data
#' @export expand.table
#' @examples
#'
#' \dontrun{
#' data(LSAT7)
#' head(LSAT7)
#' LSAT7full <- expand.table(LSAT7)
#' head(LSAT7full)
#'
#' LSAT7full <- expand.table(LSAT7, sample = TRUE)
#' head(LSAT7full)
#'
#' }
#'
expand.table <- function(tabdata, sample = FALSE) {
    if(missing(tabdata)) missingMsg('tabdata')
    if (sum(tabdata[,ncol(tabdata)]) <= nrow(tabdata))
        stop("Frequencies must be on the right of the data matrix.", call.=FALSE)
    stopifnot(is.data.frame(tabdata) || is.matrix(tabdata))
    freq <- tabdata[,ncol(tabdata)]
    tabdata <- tabdata[,-ncol(tabdata), drop=FALSE]
    fulldata <- vector('list', nrow(tabdata))
    for (i in seq_len(nrow(tabdata)))
        fulldata[[i]] <- tabdata[rep(i, freq[i]), ]
    fulldata <- do.call(rbind, fulldata)
    if(sample) fulldata <- fulldata[sample(seq_len(nrow(fulldata))), ]
    rownames(fulldata) <- seq_len(nrow(fulldata))
    fulldata
}

