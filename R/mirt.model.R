#' Specify model loadings
#'
#' The \code{mirt.model} function scans/reads user input to specify the
#' confirmatory model.
#'
#' Factors are first named and then specify which numerical items they affect
#' (i.e., where the slope is not equal to 0), separated either by commas or by
#' - to indicate a range of items. Products between factors may be specified
#' by enclosing the left hand term within brackets. To finish the declaration of
#' a model simply enter a blank line with only a carriage return (i.e., the
#' 'enter' or 'return' key), or instead read in an input version of the model syntax.
#'
#' There is an optional keyword for specifying the correlation between relationships between factors
#' called \code{COV}, and non-linear factor products can be included by enclosing the product combination
#' on the left hand side of the declaration (e.g., \code{(F1*F1)} would create a quadratic factor for
#' \code{F1}).
#'
#' \describe{
#' \item{COV}{Specify the relationship between the latent factors.
#' Estimating a correlation between factors is declared by joining the two
#' factors with an asterisk (e.g., F1*F2).}
#' \item{CONSTRAIN}{A bracketed, comma seperate list specifying equality constrains between items. 
#' The input format is 
#' \code{CONSTRAIN = (items, ..., parameterName, Group), (items, ..., parameterName, Group)}. 
#' For single group analyses, the \code{Group} specification is not required. For example, in a single group 
#' 10-item dichotmous tests, using the default 2PL model, the first and last 5 item slopes 
#' can be constrained to be equal by using \code{CONSTRAIN = (1-5, a1), (6-10, a1)}, or some cobmination
#' such as \code{CONSTRAIN = (1-3,4,5,a1), (6,7,8-10,a1)}} 
#' }
#'
#' @param input input for writing out the model syntax. Can either be a string declaration of 
#' class character or the so-called Q-matrix or class \code{matrix} that specifies the model 
#' either with integer or logical values. If the Q-matrix method 
#' is chosen covariances terms can be sepcified with the \code{COV} input
#' @param file a input specifying an external file that declares the input.
#' @param COV a symmetric, logical matrix used to declare which covariance terms are estimated
#' @param quiet logical argument passed to \code{scan()} to suppress console read message
#' @param ... additional arguments for \code{scan()}
#' @return Returns a model specification object to be used in
#' \code{\link{mirt}}, \code{\link{multipleGroup}}, or \code{\link{mixedmirt}}.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export mirt.model
#' @examples
#'
#' \dontrun{
#'
#' model <- mirt.model()
#'   F1 = 1,2,3,4-10
#'   F2 = 10-20
#'   (F1*F2) = 1,2,3,4-10
#'   COV = F1*F2
#'
#'
#' #Or alternatively
#' s <- 'F1 = 1,2,3,4-10
#'       F2 = 10-20
#'       (F1*F2) = 1,2,3,4-10
#'       COV = F1*F2'
#' model <- mirt.model(s)
#' 
#' 
#' #Q-matrix specification
#' Q <- matrix(c(1,1,1,0,0,0,0,0,0,1,1,1), ncol=2, dimnames = list(NULL, c('Factor1', 'Factor2')))
#' COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2)
#' model <- mirt.model(Q, COV=COV) 
#' 
#' ## constrain various items slopes and all intercepts in single group model to be equal
#' s <- 'F = 1-10
#'       CONSTRAIN = (1-3, 5, 6, a1), (1-10, d)'
#' model <- mirt.model(s)
#'    
#' 
#'     }
mirt.model <- function(input = NULL, file = "", COV = NULL, quiet = TRUE, ...)
{
    if(is.matrix(input)){
        fnames <- colnames(input)
        if(is.null(fnames)) fnames <- paste0('F', 1:ncol(input))
        string <- c()
        vals <- 1:nrow(input)
        input <- matrix(as.logical(input), nrow(input), ncol(input)) 
        for(i in 1L:ncol(input)){
            tmp <- vals[input[,i]]
            if(length(tmp) > 1L){
                string <- c(string, paste(c(fnames[i], ' = ', paste0(tmp[-length(tmp)], ','), 
                                            tmp[length(tmp)], '\n'), collapse=''))
            } else {
                string <- c(string, paste(c(fnames[i], ' = ', tmp[length(tmp)], '\n'), collapse=''))
            }
        }        
        if(!is.null(COV)){
            tmp <- outer(fnames, fnames, FUN=function(x, y) paste0(x, paste0('*', y)))
            sel <- upper.tri(COV) & COV
            tmp <- tmp[sel]
            if(length(tmp) > 1L){
                string <- c(string, paste(c('COV = ', paste0(tmp[-length(tmp)], ','), 
                                                      tmp[length(tmp)], '\n'), collapse=''))
            } else {
                string <- c(string, paste(c('COV = ', tmp[length(tmp)], '\n'), collapse=''))
            }            
        }
        input <- paste(string, collapse='')        
    }
    if(!is.null(input)){        
        minput <- strsplit(input, '\\n')
        file <- tempfile()
        write.table(minput, file=file, row.names=FALSE, col.names=FALSE, quote=FALSE)
    }
    mod <- scan(file = file, what = list(type = "", pars = ""),
		sep = "=", strip.white = TRUE, comment.char = "#", fill = TRUE, quiet=quiet, ...)
	mod <- cbind(mod$type, mod$pars)
	colnames(mod) <- c("Type","Parameters")
	mod <- list(x = mod)
	class(mod) <- 'mirt.model'
	mod
}

