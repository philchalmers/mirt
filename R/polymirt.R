#' Full-Information Item Factor Analysis for Mixed Data Formats
#'
#' \code{polymirt} fits an unconditional (exploratory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polytomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. The function is currently depreciated
#' and instead should be run by using the \code{\link{confmirt}} function.
#'
#' @aliases polymirt
#' @export polymirt
#' @param ... arguments to be passed to the \code{\link{confmirt}} estimation engine
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{confmirt}},
#' \code{\link{itemplot}}
polymirt <- function(...){
    .Deprecated(new = "confmirt")
    ret <- confmirt(...)
    ret
}

