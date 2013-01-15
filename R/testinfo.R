#' Function to calculate test information
#'
#' Given an estimated model compute the test information.
#' 
#' @aliases testinfo
#' @param x an estimated mirt object
#' @param Theta a matrix of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90 that jointly sum to 90. 
#' Only applicable when the input object is multidimensional
#' @param group a number signifying which group the item should be extracted from (applies to 
#' 'MultipleGroupClass' objects only)
#' @keywords information
#' @export testinfo
#' @examples 
#' 
#' \dontrun{
#' dat <- expand.table(deAyala)
#' mod <- mirt(dat, 1, '1PL')
#' 
#' Theta <- matrix(seq(-4,4,.01))
#' tinfo <- testinfo(mod, Theta)
#' plot(Theta, tinfo, type = 'l')
#' 
#' } 
testinfo <- function(x, Theta, degrees = NULL, group = NULL){
    if(is(x, 'MultipleGroupClass')) 
        J <- length(x@cmods[[1]]@pars) - 1
    else J <- length(x@pars) - 1
    info <- 0
    for(i in 1:J){
        item <- extract.item(x, i, group=group)
        info <- info + iteminfo(item, Theta=Theta, degrees=degrees)
    }
    return(info)
}
