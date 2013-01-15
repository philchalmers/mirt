#' Function to calculate test information
#'
#' Given an estimated model compute the test information.
#' 
#' @aliases testinfo
#' @param x an estimated mirt object
#' @param Theta a matrix of latent trait values
#' @param degrees a vector of angles in degrees that are between 0 and 90 that jointly sum to 90. 
#' Only applicable when the input object is multidimensional
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
testinfo <- function(x, Theta, degrees = NULL){    
    if(is(x, 'MultipleGroupClass')) 
        stop('Multiple group objects not supported. Create a 
             submodel using mirt() with fixed and non-estimated parameters.')
    J <- length(x@pars) - 1
    info <- 0
    for(i in 1:J){
        item <- extract.item(x, i)
        info <- info + iteminfo(item, Theta=Theta, degrees=degrees)
    }
    return(info)
}
