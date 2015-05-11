#' Specify model loadings
#'
#' The \code{mirt.model} function scans/reads user input to specify the
#' confirmatory model. Item locations must be used in the specifications if no
#' \code{itemnames} argument is supplied. This is called implicitly by estimation functions
#' when a string is passed to the \code{model} argument.
#'
#' Factors are first named and then specify which numerical items they affect
#' (i.e., where the slope is not equal to 0), separated either by commas or by
#' - to indicate a range of items. Products between factors may be specified
#' by enclosing the left hand term within brackets. To finish the declaration of
#' a model simply enter a blank line with only a carriage return (i.e., the
#' 'enter' or 'return' key), or instead read in an input version of the model syntax.
#'
#' There is an optional keyword for specifying the correlation between relationships between factors
#' called \code{COV}, and non-linear factor products can be included by enclosing the product
#' combination on the left hand side of the declaration (e.g., \code{(F1*F1)} would create a
#' quadratic factor for \code{F1}).
#'
#' \describe{
#'   \item{COV}{Specify the relationship between the latent factors.
#'   Estimating a correlation between factors is declared by joining the two
#'   factors with an asterisk (e.g., F1*F2), or with an asterisk between three or more factors
#'   to estimate all the possible correlations (e.g., F1*F2*F3)}
#'
#'  \item{MEAN}{A comma separated list specifying which latent factor means to freely estimate.
#'   E.g., \code{MEAN = F1, F2} will free the latent means for factors F1 and F2}
#'
#' \item{CONSTRAIN}{A bracketed, comma separated list specifying equality constrains between items.
#'   The input format is
#'   \code{CONSTRAIN = (items, ..., parameterName(s), OptionalGroup),
#'   (items, ..., parameterName, OptionalGroup)}.
#'   If \code{OptionalGroup} is omitted then the constraints are applied within all groups.
#'
#'   For example, in a single group 10-item dichotomous tests, using the default 2PL model,
#'   the first and last 5 item slopes (a1) can be constrained to be equal by using
#'   \code{CONSTRAIN = (1-5, a1), (6-10, a1)}, or some combination
#'   such as \code{CONSTRAIN = (1-3,4,5,a1), (6,7,8-10,a1)}.
#'
#'   When constraining parameters to be equal across items with different parameter names, a
#'   balanced bracketed vector must be supplied. E.g., setting the first slope for item 1 equal to
#'   the second slope in item 3 would be \code{CONSTRAIN = (1, 3, a1, a2)}
#'   }
#'
#'  \item{CONSTRAINB}{A bracketed, comma separate list specifying equality constrains between groups.
#'   The input format is \code{CONSTRAINB = (items, ..., parameterName),
#'   (items, ..., parameterName)}.
#'
#'   For example, in a two group 10-item dichotomous tests, using the default 2PL model, the first
#'   5 item slopes (a1) can be constrained to be equal across both groups by using
#'   \code{CONSTRAINB = (1-5, a1)}, or some combination such as \code{CONSTRAINB = (1-3,4,5,a1)}}
#'
#' \item{PRIOR}{A bracketed, comma separate list specifying prior parameter distributions.
#'   The input format is
#'   \code{PRIOR = (items, ..., parameterName, priorType, val1, val2, OptionalGroup),
#'   (items, ..., parameterName, priorType, val1, val2, OptionalGroup)}.
#'   If \code{OptionalGroup} is omitted then the priors are defined for all groups.
#'
#'   For example, in a single group 10-item dichotomous tests, using the default 2PL model,
#'   defining a normal prior of N(0,2) for the first 5 item intercepts (d) can be defined by
#'   \code{PRIOR = (1-5, d, norm, 0, 2)}}
#'
#' \item{LBOUND}{A bracketed, comma separate list specifying lower bounds for estimated
#'   parameters (used in optimizers such as \code{L-BFGS-B} and \code{nlminb}).
#'   The input format is \code{LBOUND = (items, ..., parameterName, value),
#'   (items, ..., parameterName, value)}.
#'
#'   For example, in a single group 10-item dichotomous tests, using the 3PL model and
#'   setting lower bounds for the 'g' parameters for the first 5 items to 0.2 is accomplished with
#'   \code{LBOUND = (1-5, g, 0.2)}}
#'
#' \item{UBOUND}{same as LBOUND, but specifying upper bounds in estimated parameters}
#'
#' \item{START}{A bracketed, comma separate list specifying the starting values for individual parameters.
#'   The input is of the form \code{(item_number, parname, value)}. For instance, setting the 10th item
#'   slope parameter (a1) to 1.0 is specified with \code{START = (1, a1, 1.0)}
#'
#'   For more hands on control of the starting values pass the argument \code{pars = 'values'} through
#'   whatever estimation function is being used}
#' }
#' @param input input for writing out the model syntax. Can either be a string declaration of
#'   class character or the so-called Q-matrix or class \code{matrix} that specifies the model
#'   either with integer or logical values. If the Q-matrix method
#'   is chosen covariances terms can be specified with the \code{COV} input
#' @param itemnames a character vector or factor indicating the item names. If a data.frame or
#'   matrix object is supplied the names will be extracted using \code{colnames(itemnames)}.
#'   Supplying this input allows the syntax to be specified with the raw item names rather than
#'   item locations
#' @param file a input specifying an external file that declares the input.
#' @param COV a symmetric, logical matrix used to declare which covariance terms are estimated
#' @param quiet logical argument passed to \code{scan()} to suppress console read message
#' @param ... additional arguments for \code{scan()}
#' @return Returns a model specification object to be used in
#'   \code{\link{mirt}}, \code{\link{bfactor}}, \code{\link{multipleGroup}}, or
#'   \code{\link{mixedmirt}}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com} and Alexander Robitzsch
#' @export mirt.model
#' @examples
#'
#' \dontrun{
#'
#' # interactively through the console (not run)
#' #model <- mirt.model()
#' #  F1 = 1,2,3,4-10
#' #  F2 = 10-20
#' #  (F1*F2) = 1,2,3,4-10
#' #  COV = F1*F2
#'
#'
#' #Or alternatively with a string input
#' s <- 'F1 = 1,2,3,4-10
#'       F2 = 10-20
#'       (F1*F2) = 1,2,3,4-10
#'       COV = F1*F2'
#' model <- mirt.model(s)
#'
#' # strings can also be passed to the estimation functions directly,
#' #   which silently calls mirt.model(). E.g., using the string above:
#' # mod <- mirt(data, s)
#'
#'
#' #Q-matrix specification
#' Q <- matrix(c(1,1,1,0,0,0,0,0,0,1,1,1), ncol=2, dimnames = list(NULL, c('Factor1', 'Factor2')))
#' COV <- matrix(c(FALSE, TRUE, TRUE, FALSE), 2)
#' model <- mirt.model(Q, COV=COV)
#'
#' ## constrain various items slopes and all intercepts in single group model to be equal,
#' #   and use a log-normal prior for all the slopes
#' s <- 'F = 1-10
#'       CONSTRAIN = (1-3, 5, 6, a1), (1-10, d)
#'       PRIOR = (1-10, a1, lnorm, .2, .2)'
#' model <- mirt.model(s)
#'
#'
#' ## constrain various items slopes and intercepts across groups for use in multipleGroup(),
#' #  and constrain first two slopes within 'group1' to be equal
#' s <- 'F = 1-10
#'       CONSTRAIN = (1-2, a1)
#'       CONSTRAINB = (1-3, 5, 6, a1), (1-10, d)'
#' model <- mirt.model(s)
#'
#'
#' ## specify model using raw item names
#' data(data.read, package = 'sirt')
#' dat <- data.read
#'
#' # syntax with variable names
#' mirtsyn2 <- "
#'        F1 = A1,B2,B3,C4
#'        F2 = A1-A4,C2,C4
#'        MEAN = F1
#'        COV = F1*F1, F1*F2
#'        CONSTRAIN=(A2-A4,a2),(A3,C2,d)
#'        PRIOR = (C3,A2-A4,a2,lnorm, .2, .2),(B3,d,norm,0,.0001)"
#' # create a mirt model
#' mirtmodel <- mirt.model(mirtsyn2, itemnames=dat)
#' # or equivelently:
#' # mirtmodel <- mirt.model(mirtsyn2, itemnames=colnames(dat))
#'
#' # mod <- mirt(dat , mirtmodel)
#'
#'     }
mirt.model <- function(input = NULL, itemnames = NULL, file = "", COV = NULL, quiet = TRUE, ...)
{
    # split_syn_string vectorized input
    split_syn_string_vec <- function( syn, vecstr ){
        for (vv in vecstr){
            syn <- split_syn_string( syn , vv  )
        }
        return(syn)
    }

    # cleans syntax in a vector from strings vv
    split_syn_string <- function( syn , vv ){
        syn <- as.list(syn )
        syn.vv <- grep( vv , syn )
        LL <- length(syn.vv)
        if (LL>0){
            for (ii in 1:LL){
                ll <- syn.vv[ii]
                syn.ll <- syn[[ll]]
                syn[[ll]] <- split_conc( syn.ll , vv )
            }
        }
        syn <- unlist(syn)
        return(syn)
    }

    # splits a string syn.ll and concatanates it with string vv
    split_conc <- function( syn.ll , vv ){
        g1 <- strsplit( syn.ll , vv , perl=FALSE )[[1]]
        Lg1 <- length(g1)
        vec <- NULL
        if (Lg1 == 1 ){ vec <- c( g1 , vv ) }
        if (Lg1 > 1 ){
            vec <- rep("" , Lg1 + (Lg1-1) )
            vec[ seq( 1 , 2*Lg1 , 2 ) ] <- g1
            vec[ seq( 2 , 2*Lg1 - 1 , 2 ) ] <- vv
            Ls1 <- nchar(syn.ll)
            if ( substring( syn.ll , Ls1 , Ls1 ) == vv ){
                vec <- c( vec , vv )
            }
        }
        return(vec)
    }

    if(!is.null(itemnames)){
        # the following block of code and above functions were largely contributed by Alexander

        # mirt syntax splitted
        inputsyntax <- input
        mirtsyn2 <- gsub( ";" , "\n" , input )
        msyn0 <- strsplit( mirtsyn2 , c(" ") )[[1]]
        syn <- msyn0
        if (is.matrix(itemnames) || is.data.frame(itemnames)){
            items <- colnames(itemnames)
        } else items <- as.character(itemnames)

        # admissible strings
        vecstr <- c( "\n" , "\\(" , "=" , "-" , "," , "\\)" )
        # process syntax
        for (vv in vecstr){
            syn <- split_syn_string( syn , vv  )
        }
        # postprocess syntax
        syn <- syn[ syn != "" ]
        syn[ syn == "\\(" ] <- "("
        syn[ syn == "\\)" ] <- ")"

        # substitute variables by numbers
        VV <- length(items)
        useditems <- NULL
        for (vv in 1:VV){
            ind <- which( syn == items[vv] )
            if (length(ind) > 0 ){
                syn[ ind ] <- vv
                useditems <- c( useditems , items[vv] )
            }
        }
        syn <- paste0( syn , collapse="")
        mirtmodel <- mirt.model(syn)
        attr(mirtmodel, 'item_syntax') <- inputsyntax
        return(mirtmodel)
    } else {
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
            if(!is.null(COV) && any(COV)){
                tmp <- outer(fnames, fnames, FUN=function(x, y) paste0(x, paste0('*', y)))
                sel <- upper.tri(COV,TRUE) & COV
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
            for(j in length(minput[[1L]]):1L){
                if(!grepl(pattern='=', minput[[1L]][j])){
                    minput[[1L]][j-1L] <- paste(minput[[1L]][j-1L], minput[[1L]][j])
                    minput[[1L]] <- minput[[1L]][-j]
                }
            }
            file <- tempfile()
            write.table(minput, file=file, row.names=FALSE, col.names=FALSE, quote=FALSE)
        }
        mod <- scan(file = file, what = list(type = "", pars = ""),
    		sep = "=", strip.white = TRUE, comment.char = "#", fill = TRUE, quiet=quiet, ...)
    	mod <- cbind(mod$type, mod$pars)
    	colnames(mod) <- c("Type","Parameters")
    	mod <- list(x = mod)
    	class(mod) <- 'mirt.model'
    	return(mod)
    }
}
