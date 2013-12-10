#' Displays item surface and information plots
#'
#' \code{itemplot} displays various item based IRT plots, with special options for plotting items
#' that contain several 0 slope parameters. Supports up to three dimensional models.
#'
#'
#' @aliases itemplot
#' @param object a computed model object of class \code{ExploratoryClass}, \code{ConfirmatoryClass}, or
#'   \code{MultipleGroupClass}. Input may also be a \code{list} for comparing similar item types (e.g., 1PL vs 2PL)
#' @param item a single numeric value, or the item name, indicating which item to plot
#' @param type plot type to use, information (\code{'info'}), standard errors (\code{'SE'}), information and
#'   standard errors (\code{'infoSE'}), item trace lines (\code{'trace'}), relative
#'   efficiency lines (\code{'RE'}), expected score \code{'score'}, or information contours \code{('infocontour')}
#'   (not for \code{MultipleGroupClass} objects)
#' @param degrees the degrees argument to be used if there are exactly two factors. See \code{\link{iteminfo}}
#'   for more detail
#' @param CE logical; plot confidence envelope?
#' @param CEalpha area remaining in the tail for confidence envelope. Default gives 95\% confidence region
#' @param CEdraws draws number of draws to use for confidence envelope
#' @param rot a list of rotation coordinates to be used for 3 dimensional plots
#' @param drop.zeros logical; drop slope values that are numerically close to zero to reduce dimensionality?
#'   Useful in objects returned from \code{\link{bfactor}} or other confirmatory models that contain several
#'   zero slopes
#' @param shiny logical; run interactive display for item plots using the \code{shiny} interface.
#'   This primarily is an instructive tool for demonstrating how item response curves
#'   behave when adjusting their parameters
#' @param ... additional arguments to be passed to \code{lattice} and \code{coef()}
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords plot
#' @export itemplot
#' @examples
#' \dontrun{
#'
#' data(LSAT7)
#' fulldata <- expand.table(LSAT7)
#' mod1 <- mirt(fulldata,1,SE=TRUE)
#' mod2 <- mirt(fulldata,1, itemtype = 'Rasch')
#' mod3 <- mirt(fulldata,2)
#'
#' itemplot(mod1, 2)
#' itemplot(mod1, 2, CE = TRUE)
#' itemplot(mod1, 2, type = 'info')
#' itemplot(mod1, 2, type = 'info', CE = TRUE)
#'
#' mods <- list(twoPL = mod1, onePL = mod2)
#' itemplot(mods, 1, type = 'RE')
#'
#' #multidimensional info
#' itemplot(mod3, 3, type = 'info')
#'
#' #polytomous items
#' pmod <- mirt(Science, 1, SE=TRUE, SE.type = 'MHRM')
#' itemplot(pmod, 3)
#' itemplot(pmod, 3, CE = TRUE)
#'
#' #interactive shiny applet
#' itemplot(shiny = TRUE)
#'     }
#'
itemplot <- function(object, item, type = 'trace', degrees = 45, CE = FALSE, CEalpha = .05,
                     CEdraws = 1000, drop.zeros = FALSE, rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                     shiny = FALSE, ...){
    if(shiny){
        require(shiny)
        runApp(shinyItemplot())
    }
    if(is.list(object)) inames <- colnames(object[[1]]@data)
    else inames <- colnames(object@data)
    ind <- 1:length(inames)
    if(!is.numeric(item)) item <- ind[inames == item]
    rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
    ret <- itemplot.internal(object=object, item=item, type=type, degrees=degrees, CE=CE,
                             CEalpha=CEalpha, CEdraws=CEdraws, drop.zeros=drop.zeros, rot=rot, ...)
    if(is.null(ret)) return(invisible(ret))
    else return(ret)
}
