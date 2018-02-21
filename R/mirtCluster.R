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
#' @import future
#' @import future.apply
#' @param spec input that is passed to \code{parallel::makeCluster()}. If no input is given the
#'   maximum number of available local cores will be used
#' @param ... additional arguments to pass to \code{parallel::makeCluster}
#' @param remove logical; remove previously defined \code{mirtCluster()}?
#' @param use_future logical; run asynchronous parallel works using \code{future::future()}
#'    and \code{future.apply::future_*apply} internally?
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
#' #make 4 cores available for parallel computing
#' mirtCluster(4)
#'
#' #make remote clusters for parallel computing (SSH setup required)
#' mirtCluster(c('n1', 'n2', 'n3'))
#'
#' #make remove clusters for asynchronous parallel computing using \code{future::future} API (SSH setup required)
#' mirtCluster(c('n1', 'n2', 'n3), use_future = TRUE)
#'
#' #' #stop and remove cores
#' mirtCluster(remove = TRUE)
#'
#' }
mirtCluster <- function(spec, ..., remove = FALSE, use_future = FALSE){
    if(requireNamespace("parallel", quietly = TRUE)){
        if(missing(spec))
            spec <- parallel::detectCores()
        if(remove){
            # remove parallel
            if(is.null(.mirtClusterEnv$MIRTCLUSTER)){
                message('There is no visible mirtCluster() definition')
                return(invisible())
            }
            parallel::stopCluster(.mirtClusterEnv$MIRTCLUSTER)
            .mirtClusterEnv$MIRTCLUSTER <- NULL
            .mirtClusterEnv$ncores <- 1L
            .mirtClusterEnv$future <- FALSE
            return(invisible())
        }
        if(!is.null(.mirtClusterEnv$MIRTCLUSTER)){
            # if already called, then return a warning
            message('mirtCluster() has already been defined')
            return(invisible())
        }

        if(use_future){
            .mirtClusterEnv$future <- TRUE
        }
        .mirtClusterEnv$MIRTCLUSTER <- parallel::makeCluster(spec, ...) # future also use parallel::makeCluster()
        .mirtClusterEnv$ncores <- length(.mirtClusterEnv$MIRTCLUSTER)

        if(.mirtClusterEnv$future){
            if(((requireNamespace('future', quietly = TRUE)))){
                if(is.numeric(.mirtClusterEnv$MIRTCLUSTER)){
                    future::plan(list(
                        future::tweak(future::multisession,
                                      workers = function() { max(1, round(.mirtClusterEnv$ncores)) }, gc = TRUE) # use ncores in local
                    ))
                } else {
                    # connect PSOCK (SSH)
                    future::plan(list(
                        future::tweak(future::cluster, workers = .mirtClusterEnv$MIRTCLUSTER, gc = TRUE), # connect remote first
                        future::tweak(future::multisession,
                                      workers = function() { max(1, round(0.7 * future::availableCores())) }, gc = TRUE) # then detect cores each remote
                    ))
                }
            }
        }
        mySapply(1L:.mirtClusterEnv$ncores*2L, function(x) invisible())
    }
    return(invisible())
}
