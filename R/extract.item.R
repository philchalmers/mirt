#' Extract an item object from mirt objects
#'
#' Extract the internal mirt objects from any estimated model.
#' 
#' @aliases extract.item
#' @param x mirt model of class 'ExploratoryClass', 'ConfirmatoryClass', or 'MultipleGroupClass'
#' @param item a number signifying which item to extract
#' @param group a number signifying which group the item should be extracted from (applies to 
#' 'MultipleGroupClass' only)
#' @keywords extract
#' @export extract.item
#' @examples 
#' 
#' \dontrun{
#' mod <- mirt(Science, 1)
#' extr.1 <- extract.item(mod, 1)
#' }
extract.item <- function(x, item, group = NULL){
    if(is(x, 'MultipleGroupClass')){
        if(is.null(group)) stop('Which group are you trying to extract from?')
        ret <- x@cmods[[group]]@pars[[item]]
    } else {
        ret <- x@pars[[item]]
    }
    ret
}
