#' Define a parallel cluster object to be used in internal functions
#'
#' This function defines a object that is placed in a relevant internal environment defined in mirt.
#' Internal functions such as \code{calcLogLik}, \code{fscores}, etc, will utilize this object
#' automatically to capitalize on parallel
#' processing architecture. The object defined is a call from \code{parallel::makeCluster()}.
#' Note that if you are defining other parallel objects (for simulation designs, for example)
#' it is not recommended to define a mirtCluster.
#'
#' @aliases mirtCluster
#' @param spec input that is passed to \code{parallel::makeCluster()}. If no input is given the
#'   maximum number of available local cores will be used. Setting this to NULL will skip a new definition (allows \code{omp_threads} to be used independently)
#' @param omp_threads number of OpenMP threads to use (currently applies to E-step computations only).
#'   Not used when argument input is missing
#' @param ... additional arguments to pass to \code{parallel::makeCluster}
#' @param remove logical; remove previously defined \code{mirtCluster()}?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords parallel
#' @export mirtCluster
#' @examples
#'
#' \dontrun{
#'
#' # use all available cores
#' mirtCluster()
#' mirtCluster(remove = TRUE)
#'
#' #make 4 cores available for parallel computing
#' mirtCluster(4)
#' mirtCluster(remove = TRUE)
#'
#' # create 3 core architecture in R, and 4 thread architecture with OpenMP
#' mirtCluster(spec = 3, omp_threads = 4)
#'
#' # leave previous multicore objects, but change omp_threads
#' mirtCluster(spec = NULL, omp_threads = 2)
#'
#' }
mirtCluster <- function(spec, omp_threads, remove = FALSE, ...){
    if(!missing(omp_threads)){
        stopifnot(is.numeric(omp_threads))
        .mirtClusterEnv$omp_threads <- omp_threads
    }
    if(!missing(spec) && is.null(spec))
        return(invisible(NULL))
    if(requireNamespace("parallel", quietly = TRUE)){
        if(missing(spec))
            spec <- parallel::detectCores()
        if(remove){
            if(is.null(.mirtClusterEnv$MIRTCLUSTER)){
                message('There is no visible mirtCluster() definition')
                return(invisible(NULL))
            }
            parallel::stopCluster(.mirtClusterEnv$MIRTCLUSTER)
            .mirtClusterEnv$MIRTCLUSTER <- NULL
            .mirtClusterEnv$ncores <- 1L
            .mirtClusterEnv$omp_threads <- 1L
            return(invisible(NULL))
        }
        if(!is.null(.mirtClusterEnv$MIRTCLUSTER)){
            message(
                sprintf('mirtCluster() previously defined for %i clusters',
                        .mirtClusterEnv$ncores))
            return(invisible(NULL))
        }
        .mirtClusterEnv$MIRTCLUSTER <- parallel::makeCluster(spec, ...)
        .mirtClusterEnv$ncores <- length(.mirtClusterEnv$MIRTCLUSTER)
        mySapply(1L:.mirtClusterEnv$ncores*2L, function(x) invisible(NULL))
    }
    return(invisible(NULL))
}
