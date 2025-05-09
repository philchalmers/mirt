#' Compute multidimensional difficulty index
#'
#' Returns a matrix containing the MDIFF values (Reckase, 2009). Only supported for items of class
#' 'dich' and 'graded'.
#'
#' @aliases MDIFF
#' @param x an object of class 'SingleGroupClass', or an object of class 'MultipleGroupClass' if a suitable
#'   \code{group} input were supplied
#' @param which.items a vector indicating which items to select. If NULL is used
#'   (the default) then MDISC will be computed for all items
#' @param group group argument to pass to \code{\link{extract.group}} function. Required when the input object is
#'   a multiple-group model
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Reckase, M. D. (2009). Multidimensional Item Response Theory. Springer.
#'
#' @seealso \code{\link{extract.group}}, \code{\link{MDISC}}
#'
#' @keywords discrimination
#' @export MDIFF
#' @examples
#' \donttest{
#'
#' mod <- mirt(Science, 2)
#' MDIFF(mod)
#'
#' mod <- mirt(expand.table(LSAT7), 2)
#' MDIFF(mod)
#'
#' }
MDIFF <- function(x, which.items = NULL, group=NULL){
    if(missing(x)) missingMsg('x')
    if(is(x, 'MultipleGroupClass') && is.null(group))
        stop('Input must be a SingleGroupClass object or a MultipleGroupClass object with a suitable group input',
             call.=FALSE)
    if(!is.null(group) && is(x, 'MultipleGroupClass'))
        x <- extract.group(x=x, group=group)
    stopifnot(class(x) == 'SingleGroupClass')
    J <- extract.mirt(x, 'nitems')
    if(is.null(which.items)) which.items <- 1L:J
    out <- vector('list', length(which.items))
    MD <- MDISC(x)
    for(i in seq_len(length(which.items))){
        item <- extract.item(x, which.items[i])
        if(!(class(item) %in% c('dich', 'graded')))
            stop(sprintf('Item %i is not of class \"graded\" or \"dich\"', which.items[i]),
                 call.=FALSE)
        ds <- ExtractZetas(item)
        out[[i]] <- -ds / MD[which.items[i]]
    }
    ret <- matrix(NA, length(out), max(sapply(out, length)))
    for(i in seq_len(length(out)))
        ret[i,1L:length(out[[i]])] <- out[[i]]
    rownames(ret) <- extract.mirt(x, 'itemnames')[which.items]
    colnames(ret) <- paste0('MDIFF_', 1L:ncol(ret))
    ret
}
