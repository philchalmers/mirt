#' Specify model loadings
#' 
#' The \code{specifyModel} function scans user input to specify the
#' confirmatory model.
#' 
#' Factors are first named and then specify which numerical items they affect
#' (i.e., where the slope is not equal to 0), separated either by commas or by
#' - to indicate a range of items. Products between factors may be specified 
#' by enclosing the left hand term within brackets. To finish the declaration of
#' a model simply enter a blank line with only a carriage return (i.e., the 
#' 'enter' or 'return' key).
#'
#' There is an optional keyword for specifying the correlation between relationships between factors 
#' called \code{COV}, and nonlinear factor products can be included by enclosing the product combination
#' on the left hand side of the declaration (e.g., \code{(F1*F1)} would create a quadratic factor for 
#' \code{F1}). 
#' 
#' \describe{ 
#' \item{COV}{Specify the relationship between the latent factors.
#' Estimating a correlation between factors is declared by joining the two
#' factors with an asterisk (e.g., F1*F2), while fixing this value to a constant
#' is performed by using the keyword 'eq' followed by some numeric value (e.g.,
#' \code{F1*F2 eq .5} fixes the correlation between F1 and F2 to .5).}
#' }
#' 
#' @param file a string specifying an external file that declares the input.
#' @param ... additional arguments for \code{scan()}
#' @return Returns a model specification object to be used in
#' \code{\link{confmirt}}.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export specifyModel
#' @examples
#' 
#' \dontrun{
#' 
#' model <- specifyModel()
#'   F1 = 1,2,3,4-10
#'   F2 = 10-20
#'   (F1*F2) = 1,2,3,4-10
#'   COV = F1*F2
#'     
#'     }
#' 
specifyModel <- function(file = "", ...)
{
    mod <- scan(file = file, what = list(type = "", pars = ""), 
		sep = "=", strip.white = TRUE, comment.char = "#", fill = TRUE, ...)
	mod <- cbind(mod$type, mod$pars)
	colnames(mod) <- c("Type","Parameters")	
	mod <- list(x = mod)
	class(mod) <- 'confmirt.model'
	mod
}

