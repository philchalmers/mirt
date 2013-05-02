#' Function to calculate probability trace lines
#'
#' Given an internal mirt object extracted from an estimated model
#' compute the probability trace lines for all categories.
#'
#' @aliases probtrace
#' @param x an extracted internal mirt object containing item information
#' @param Theta a vector (unidimensional) or matrix (multidimensional) of latent trait values
#' @keywords tracelines
#' @export probtrace
#' @seealso
#' \code{\link{extract.item}}
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#' extr.2 <- extract.item(mod, 2)
#' Theta <- matrix(seq(-4,4, by = .1))
#' traceline <- probtrace(extr.2, Theta)
#'
#' head(data.frame(traceline, Theta=Theta))
#'
#' }
probtrace <- function(x, Theta){
    if(is(Theta, 'vector')) Theta <- as.matrix(Theta)
    if(!is.matrix(Theta)) stop('Theta input must be a matrix')
    P <- ProbTrace(x=x, Theta=Theta)
    cats <- 1:ncol(P)
    if(ncol(P) == 2) cats <- cats - 1
    colnames(P) <- paste0('P.', cats)
    P
}
