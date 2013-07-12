#' Define a cluster object to be used in mirt functions
#'
#' This function defines a object that is placed in the users workspace 
#' (i.e., the\code{.GlobalEnv}). Relavent internal functions such as \code{calcLogLik},
#' \code{fscores}, etc, will utilize this object automatically to capitilze on parallel
#' processing architecture. The object defined is a call from \code{parallel::makeCluster()}.
#' @aliases mirtCluster
#' @param ncores number of cores to be used in the returned object which is 
#' passed to \code{parallel::makeCluster()}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords parallel
#' @export mirtCluster
#' @examples
#'
#' \dontrun{
#' #make 4 cores available for parallel computing
#' mirtCluster(4)
#' 
#' }
mirtCluster <- function(ncores){
    if(!is.null(globalenv()$MIRTCLUSTER))
        stop('MIRTCLUSTER object already defined')
    if(!require(parallel)) require(parallel)
    if(!is.numeric(ncores)) 
        stop('ncores must be numeric')
    .GlobalEnv$MIRTCLUSTER <- makeCluster(ncores)
    return(invisible())
}
