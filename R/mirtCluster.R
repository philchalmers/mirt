#' Define a parallel cluster object to be used in internal functions
#'
#' This function defines a object that is placed in the users workspace 
#' (i.e., the\code{.GlobalEnv}). Relavent internal functions such as \code{calcLogLik},
#' \code{fscores}, etc, will utilize this object automatically to capitilze on parallel
#' processing architecture. The object defined is a call from \code{parallel::makeCluster()} and 
#' defines an object called \code{MIRTCLUSTER}.
#' @aliases mirtCluster
#' @param ncores number of cores to be used in the returned object which is 
#' passed to \code{parallel::makeCluster()}. If no input is given the maximum number of available
#' cores will be used
#' @param remove logical; remove previously defined \code{mirtCluster()}?
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
#' #use all available cores
#' mirtCluster()
#' 
#' }
mirtCluster <- function(ncores, remove = FALSE){
    if(!require(parallel)) require(parallel)    
    if(remove){
        if(is.null(globalenv()$MIRTCLUSTER)){
            message('No MIRTCLUSTER object in workspace')
            return(invisible())
        }
        stopCluster(.GlobalEnv$MIRTCLUSTER)        
        .GlobalEnv$MIRTCLUSTER <- NULL
        return(invisible())
    }
    if(!is.null(globalenv()$MIRTCLUSTER)){
        message('MIRTCLUSTER object already defined')    
        return(invisible())
    }
    if(missing(ncores))
        ncores <- detectCores()
    if(!is.numeric(ncores)) 
        stop('ncores must be numeric')    
    .GlobalEnv$MIRTCLUSTER <- makeCluster(ncores)
    return(invisible())
}
