#' Define a parallel cluster object to be used in internal functions
#'
#' This function defines a object that is placed in a relevant internal environment defined in mirt.
#' Internal functions such as \code{calcLogLik}, \code{fscores}, etc, will utilize this object
#' automatically to capitalize on parallel
#' processing architecture. The object defined is a call from \code{parallel::makeCluster()}.
#' Note that if you are defining other parallel objects (for simulation desings, for example)
#' it is not recommended to define a mirtCluster.
#'
#' @aliases mirtCluster
#' @param spec input that is passed to \code{parallel::makeCluster()}. If no input is given the
#'   maximum number of available local cores will be used
#' @param ... additional arguments to pass to \code{parallel::makeCluster}
#' @param remove logical; remove previously defined \code{mirtCluster()}?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords parallel
#' @export mirtCluster
#' @examples
#'
#' \dontrun{
#'
#' #make 4 cores available for parallel computing
#' mirtCluster(4)
#'
#' #' #stop and remove cores
#' mirtCluster(remove = TRUE)
#'
#' }
mirtCluster <- function(spec, ..., remove = FALSE){
    if(requireNamespace("parallel", quietly = TRUE)){
        if(missing(spec))
            spec <- parallel::detectCores()
        if(remove){
            if(is.null(.mirtClusterEnv$MIRTCLUSTER)){
                message('There is no visible mirtCluster() definition')
                return(invisible())
            }
            parallel::stopCluster(.mirtClusterEnv$MIRTCLUSTER)
            .mirtClusterEnv$MIRTCLUSTER <- NULL
            .mirtClusterEnv$ncores <- 1L
            return(invisible())
        }
        if(!is.null(.mirtClusterEnv$MIRTCLUSTER)){
            message('mirtCluster() has already been defined')
            return(invisible())
        }
        .mirtClusterEnv$MIRTCLUSTER <- parallel::makeCluster(spec, ...)
        .mirtClusterEnv$ncores <- length(.mirtClusterEnv$MIRTCLUSTER)
    }
    return(invisible())
}
