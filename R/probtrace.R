#' Function to calculate probability trace lines
#'
#' Given an internal mirt object extracted from an estimated model, or the
#' single-group estimated model itself, compute the probability trace
#' lines for all categories.
#'
#' @aliases probtrace
#' @param x either an extracted internal mirt object containing item information
#'   (see \code{\link{extract.item}}) or a model of class \code{SingleGroupClass}
#'   typically returned by the function \code{\link{mirt}} or \code{\link{bfactor}}
#' @param Theta a vector (unidimensional) or matrix (unidimensional/multidimensional) of
#'   latent trait values
#' @keywords tracelines
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export probtrace
#' @seealso
#' \code{\link{extract.item}}, \code{\link{extract.group}}
#' @examples
#'
#' mod <- mirt(Science, 1)
#'
#' # single item probabilty tracelines for Item 2
#' extr.2 <- extract.item(mod, 2)
#' Theta <- matrix(seq(-4,4, by = .1))
#' traceline <- probtrace(extr.2, Theta)
#' head(data.frame(traceline, Theta=Theta))
#'
#' # probability tracelines for all items in test
#' tracelines <- probtrace(mod, Theta)
#' head(tracelines)
#'
probtrace <- function(x, Theta){
    if(is(x, 'SingleGroupClass')){
        nitems <- extract.mirt(x, 'nitems')
        itemnames <- extract.mirt(x, 'itemnames')
        ret <- vector('list', nitems)
        ret <- lapply(seq_len(nitems), function(i){
            P <- probtrace(x=x@ParObjects$pars[[i]], Theta=Theta)
            colnames(P) <- paste0(itemnames[i], '.', colnames(P))
            P
        })
        return(do.call(cbind, ret))
    } else {
        if(missing(x)) missingMsg('x')
        if(missing(Theta)) missingMsg('Theta')
        if(is(Theta, 'vector')) Theta <- as.matrix(Theta)
        if(!is.matrix(Theta)) stop('Theta input must be a matrix', call.=FALSE)
        if(ncol(Theta) != x@nfact)
            stop('Theta does not have the correct number of dimensions', call.=FALSE)
        P <- ProbTrace(x=x, Theta=Theta)
        cats <- 1:ncol(P)
        if(ncol(P) == 2) cats <- cats - 1
        colnames(P) <- paste0('P.', cats)
        return(P)
    }
}
