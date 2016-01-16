#' Compute multidimensional difficulty index
#'
#' Returns a matrix containing the MDIFF values (Reckase, 2009).
#'
#' @aliases MDIFF
#' @param x an object of class 'SingleGroupClass'
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Reckase, M. D. (2009). Multidimensional Item Response Theory. Springer.
#'
#' @seealso \code{\link{extract.group}}, \code{\link{MDISC}}
#'
#' @keywords discrimination
#' @export MDIFF
#' @examples
#' \dontrun{
#'
#' mod <- mirt(Science, 2)
#' MDIFF(mod)
#'
#' }
MDIFF <- function(x){
    if(missing(x)) missingMsg('x')
    stopifnot(class(x) == 'SingleGroupClass')
    out <- vector('list', extract.mirt(x, 'nitems'))
    MD <- MDISC(x)
    for(i in 1L:length(out)){
        item <- extract.item(x, i)
        ds <- ExtractZetas(item)
        out[[i]] <- -ds / MD[i]
    }
    ret <- matrix(NA, length(out), max(sapply(out, length)))
    for(i in 1L:length(out))
        ret[i,1L:length(out[[i]])] <- out[[i]]
    rownames(ret) <- extract.mirt(x, 'itemnames')
    colnames(ret) <- paste0('MDIFF_', 1L:ncol(ret))
    ret
}
