#' Compute multidimensional difficulty index
#'
#' Returns a matrix containing the MDIFF values (Reckase, 2009). Only suppored for items of class
#' 'dich' and 'graded'.
#'
#' @aliases MDIFF
#' @param x an object of class 'SingleGroupClass'
#' @param which.items a vector indicating which items to select. If NULL is used
#'   (the default) then MDISC will be computed for all items
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
MDIFF <- function(x, which.items = NULL){
    if(missing(x)) missingMsg('x')
    stopifnot(class(x) == 'SingleGroupClass')
    J <- extract.mirt(x, 'nitems')
    if(is.null(which.items)) which.items <- 1L:J
    out <- vector('list', length(which.items))
    MD <- MDISC(x)
    for(i in 1L:length(out)){
        item <- extract.item(x, which.items[i])
        if(!(class(item) %in% c('dich', 'graded')))
            stop(sprintf('Item %i is not of class \"graded\" or \"dich\"', which.items[i]))
        ds <- ExtractZetas(item)
        out[[i]] <- -ds / MD[which.items[i]]
    }
    ret <- matrix(NA, length(out), max(sapply(out, length)))
    for(i in 1L:length(out))
        ret[i,1L:length(out[[i]])] <- out[[i]]
    rownames(ret) <- extract.mirt(x, 'itemnames')[which.items]
    colnames(ret) <- paste0('MDIFF_', 1L:ncol(ret))
    ret
}
