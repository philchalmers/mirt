#' Full-Information Item Factor Analysis for Mixed Data Formats
#' 
#' \code{polymirt} fits an unconditional (exploratory) full-information
#' maximum-likelihood factor analysis model to dichotomous and polychotomous
#' data under the item response theory paradigm using Cai's (2010)
#' Metropolis-Hastings Robbins-Monro algorithm. If requested, lower and upper asymptote
#' parameters are estimated with a beta priors included automatically.
#'
#'
#' @param ... arguments to be passed to the \code{\link{confmirt}} estimation engine
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{expand.table}}, \code{\link{key2binary}}, \code{\link{confmirt}},
#' \code{\link{itemplot}}
polymirt <- function(...){         
    ret <- confmirt(...)    
    message('NOTE: polymirt() is now obsolete and will be removed completely in version 0.3.0,  
         use confmirt(data, nfact) for exploratory models instead.')    
    ret    
}

