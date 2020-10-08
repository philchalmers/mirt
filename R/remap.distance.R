#' Remap item categories to have integer distances of 1
#'
#' The mirt package's estimation setup requires that all item responses have spaces
#' equal to 1 (e.g., a Likert scale scored from 1 through 5). In the event that categories
#' are missing the categories must be re-coded. This function is automatically called by
#' the package estimation functions (e.g., \code{\link{mirt}}), however for convince this
#' function has been extracted for users to better understand the remapping consequences.
#'
#' @param data the response data to remap as a data.frame or matrix
#' @param message logical; print message information pertaining to which items were remapped?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export
#' @examples
#'
#' # category 2 for item 1 missing
#' dat <- Science
#' dat[,1] <- ifelse(Science[,1] == 2, 1, Science[,1])
#' apply(dat, 2, table)
#'
#' # mirt() automatically remaps categories
#' mod <- mirt(dat, 1)
#' coef(mod, simplify=TRUE)
#'
#' # this is the transformed data used by mirt()
#' remap_dat <- remap.distance(dat)
#' apply(remap_dat, 2, table)
#'
#'
remap.distance <- function(data, message = TRUE) {
    if(missing(data)) missingMsg('dat')
    ind <- seq_len(ncol(data))
    nms <- colnames(data)
    if(is.null(nms)) nms <- paste0("Item ", ind)
    ret <- sapply(ind, function(i, data, nms, message){
        x <- data[,i]
        s <- sort(unique(x))
        se <- min(s, na.rm = TRUE):max(x, na.rm = TRUE)
        if(length(s) != length(se)){
            if(message)
                message(sprintf('\"%s\" re-mapped to ensure all categories have a distance of 1', nms[i]))
            for(i in 2L:length(s))
                x <- ifelse(x == s[i], se[i], x)
        }
        x
    }, data=data, nms=nms, message=message)
    rownames(ret) <- rownames(data)
    colnames(ret) <- colnames(data)
    ret
}

