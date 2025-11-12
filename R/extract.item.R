#' Extract an item object from \code{mirt} objects
#'
#' Extract the internal \code{mirt} objects from any estimated model.
#'
#' @aliases extract.item
#' @param x mirt model of class \code{'SingleGroupClass'}, \code{'MultipleGroupClass'}, or
#'   \code{'MixtureClass'}
#' @param item a number or character signifying which item to extract
#' @param group which group the item should be extracted from (applies to
#'   \code{'MultipleGroupClass'} and \code{'MixtureClass'} only). Can be a numeric
#'   value or the name of the group to be extracted
#' @param drop.zeros logical; drop slope values that are numerically close to zero to reduce
#'   dimensionality? Useful in objects returned from \code{\link{bfactor}} or other confirmatory
#'   models that contain several zero slopes
#' @keywords extract
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @seealso \code{\link{extract.group}}, \code{\link{extract.mirt}}
#' @export extract.item
#' @examples
#'
#' \donttest{
#' mod <- mirt(Science, 1)
#' extr.1 <- extract.item(mod, 1)
#' }
extract.item <- function(x, item, group = NULL, drop.zeros = FALSE){
    if(missing(x)) missingMsg('x')
    if(missing(item)) missingMsg('item')
    if(is(x, 'MixedClass'))
        stop('Lower-level functions do not support extracted items from MixedClass objects',
             call.=FALSE)
    inames <- extract.mirt(x, 'itemnames')
    ind <- 1L:length(inames)
    if(!is.numeric(item)) item <- ind[inames == item]
    if(is(x, 'MultipleGroupClass') || is(x, 'MixtureClass')){
        if(is.null(group)) stop('Which group are you trying to extract from?', call.=FALSE)
        if(is.character(group)){
            grp <- extract.mirt(x, 'groupNames')
            group <- which(group == grp)
        }
        ret <- extract.mirt(extract.mirt(x, 'pars')[[group]], 'pars')[[item]]
    } else {
        ret <- extract.mirt(x, 'pars')[[item]]
    }
    if(drop.zeros){
        zeros <- ret@par > -1e-8 & ret@par < 1e-8
        nfact <- extract.mirt(x, 'nfact')
        zeros[-c(1L:nfact)] <- FALSE
        ret@par <- ret@par[!zeros]
        ret@est <- ret@est[!zeros]
        ret@parnum <- ret@parnum[!zeros]
        ret@nfact <- sum(!zeros[c(1L:nfact)])
    }
    ret
}
