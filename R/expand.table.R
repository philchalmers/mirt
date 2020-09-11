#' Expand summary table of patterns and frequencies
#'
#' The \code{expand.table} function expands a summary table of unique response
#' patterns to a full sized data-set. By default the response frequencies are
#' assumed to be on rightmost column of the input data, though this can be modified.
#'
#' @param tabdata An object of class \code{data.frame} or \code{matrix}
#'   with the unique response patterns and the number of frequencies
#'   in the rightmost column (though see \code{freq} for details on how to omit this
#'   column)
#' @param freq either a character vector specifying the column in \code{tabdata}
#'   to be used as the frequency count indicator for each response pattern (defaults to
#'   the right-most column) or a integer vector of length \code{nrow(tabdata)} specifying
#'   the frequency counts. When using the latter approach the \code{tabdata} input should not
#'   include any information regarding the counts, and instead should only include the unique
#'   response patterns themselves
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
#' data(LSAT7)
#' head(LSAT7) # frequency in right-most column
#' LSAT7full <- expand.table(LSAT7)
#' head(LSAT7full)
#' dim(LSAT7full)
#'
#' # randomly switch rows in the expanded response table
#' LSAT7samp <- expand.table(LSAT7, sample = TRUE)
#' head(LSAT7samp)
#' colMeans(LSAT7full)
#' colMeans(LSAT7samp) #equal
#'
#' #--------
#'
#' \dontrun{
#' # Generate data from separate response pattern matrix and freq vector
#' # The following uses Table 2.1 from de Ayala (2009)
#' f <- c(691,2280,242,235,158,184,1685,1053,134,462,92,65,571,79,87,41,1682,702,
#'        370,63,626,412,166,52,28,15,2095,1219,500,187,40,3385)
#'
#' longpat <- scan()
#' 0 0 0 0 0
#' 1 0 0 0 0
#' 0 1 0 0 0
#' 0 0 1 0 0
#' 0 0 0 1 0
#' 0 0 0 0 1
#' 1 1 0 0 0
#' 1 0 1 0 0
#' 0 1 1 0 0
#' 1 0 0 1 0
#' 0 1 0 1 0
#' 0 0 1 1 0
#' 1 0 0 0 1
#' 0 1 0 0 1
#' 0 0 1 0 1
#' 0 0 0 1 1
#' 1 1 1 0 0
#' 1 1 0 1 0
#' 1 0 1 1 0
#' 0 1 1 1 0
#' 1 1 0 0 1
#' 1 0 1 0 1
#' 1 0 0 1 1
#' 0 1 1 0 1
#' 0 1 0 1 1
#' 0 0 1 1 1
#' 1 1 1 1 0
#' 1 1 1 0 1
#' 1 1 0 1 1
#' 1 0 1 1 1
#' 0 1 1 1 1
#' 1 1 1 1 1
#'
#'
#' pat <- matrix(longpat, byrow=TRUE, ncol=5)
#' colnames(pat) <- paste0('Item.', 1:5)
#' head(pat)
#'
#' table2.1 <- expand.table(pat, freq = f)
#' dim(table2.1)
#'
#' }
#'
#'
expand.table <- function(tabdata, freq = colnames(tabdata)[ncol(tabdata)],
                         sample = FALSE) {
    if(missing(tabdata)) missingMsg('tabdata')
    if(is.null(colnames(tabdata)) && is.null(freq))
        stop('Please either supply colnames to tabdata or provide a vector of counts in freq', call.=FALSE)
    stopifnot(is.data.frame(tabdata) || is.matrix(tabdata))
    tabdat <- as.matrix(tabdata)
    if(is.character(freq)){
        stopifnot(length(freq) == 1L)
        tmp <- tabdata[,freq]
        tabdata <- tabdata[, colnames(tabdata) != freq, drop=FALSE]
        freq <- tmp
    }
    stopifnot(length(freq) == nrow(tabdata))
    fulldata <- vector('list', nrow(tabdata))
    for (i in seq_len(nrow(tabdata)))
        fulldata[[i]] <- tabdata[rep(i, freq[i]), ]
    fulldata <- do.call(rbind, fulldata)
    if(sample) fulldata <- fulldata[sample(seq_len(nrow(fulldata))), ]
    rownames(fulldata) <- seq_len(nrow(fulldata))
    fulldata
}

