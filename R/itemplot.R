#' Displays item surface and information plots
#'
#' \code{itemplot} displays various item based IRT plots, with special options for plotting items
#' that contain several 0 slope parameters. Supports up to three dimensional models.
#'
#'
#' @aliases itemplot
#' @param object a computed model object of class \code{SingleGroupClass}
#'   or \code{MultipleGroupClass}. Input may also be a \code{list} for comparing similar item types
#'   (e.g., 1PL vs 2PL)
#' @param item a single numeric value, or the item name, indicating which item to plot
#' @param type plot type to use, information (\code{'info'}), standard errors (\code{'SE'}),
#'   item trace lines (\code{'trace'}), information and standard errors (\code{'infoSE'}) or
#'   information and trace lines (\code{'infotrace'}), relative efficiency lines (\code{'RE'}),
#'   expected score \code{'score'}, or information and trace line contours (\code{'infocontour'} and
#'   \code{'tracecontour'}; not supported for \code{MultipleGroupClass} objects)
#' @param degrees the degrees argument to be used if there are two or three factors.
#'   See \code{\link{iteminfo}} for more detail. A new vector will be required for three dimensional
#'   models to override the default
#' @param CE logical; plot confidence envelope?
#' @param CEalpha area remaining in the tail for confidence envelope. Default gives 95\% confidence
#'   region
#' @param CEdraws draws number of draws to use for confidence envelope
#' @param rot a list of rotation coordinates to be used for 3 dimensional plots
#' @param drop.zeros logical; drop slope values that are numerically close to zero to reduce
#'   dimensionality? Useful in objects returned from \code{\link{bfactor}} or other confirmatory
#'   models that contain several zero slopes
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param npts number of quadrature points to be used for plotting features.
#'   Larger values make plots look smoother
#' @param shiny logical; run interactive display for item plots using the \code{shiny} interface.
#'   This primarily is an instructive tool for demonstrating how item response curves
#'   behave when adjusting their parameters
#' @param auto.key plotting argument passed to \code{\link{lattice}}
#' @param par.strip.text plotting argument passed to \code{\link{lattice}}
#' @param par.settings plotting argument passed to \code{\link{lattice}}
#' @param ... additional arguments to be passed to \code{\link{lattice}} and \code{coef()}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
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
#' #multidimensional
#' itemplot(mod3, 4, type = 'info')
#' itemplot(mod3, 4, type = 'infocontour')
#' itemplot(mod3, 4, type = 'tracecontour')
#'
#' #polytomous items
#' pmod <- mirt(Science, 1, SE=TRUE)
#' itemplot(pmod, 3)
#' itemplot(pmod, 3, CE = TRUE)
#' itemplot(pmod, 3, type = 'score')
#' itemplot(pmod, 3, type = 'infotrace')
#'
#' # use the directlabels package to put labels on tracelines
#' library(directlabels)
#' plt <- itemplot(pmod, 3)
#' direct.label(plt, 'top.points')
#'
#' # change colour theme of plots
#' bwtheme <- standard.theme("pdf", color=FALSE)
#' plot(pmod, type='trace', par.settings=bwtheme)
#' itemplot(pmod, 1, type = 'trace', par.settings=bwtheme)
#'
#' itemplot(pmod, 1, type = 'infoSE')
#' update(trellis.last.object(), par.settings = bwtheme)
#'
#' # uncomment to run interactive shiny applet
#' # itemplot(shiny = TRUE)
#'     }
#'
itemplot <- function(object, item, type = 'trace', degrees = 45, CE = FALSE, CEalpha = .05,
                     CEdraws = 1000, drop.zeros = FALSE, theta_lim = c(-6,6), shiny = FALSE,
                     rot = list(xaxis = -70, yaxis = 30, zaxis = 10),
                     par.strip.text = list(cex = 0.7), npts = 200,
                     par.settings = list(strip.background = list(col = '#9ECAE1'),
                                         strip.border = list(col = "black")),
                     auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...){
    if(shiny){
        if(requireNamespace("shiny", quietly = TRUE)){
            shiny::runApp(shinyItemplot(), ...)
        } else {
            stop('shiny package is not available. Please install.', call.=FALSE)
        }
    }
    dots <- list(...)
    if(missing(object)) missingMsg('object')
    if(missing(item)) missingMsg('item')
    if(is(object, 'MixtureClass')) class(object) <- "MultipleGroupClass"
    if(is(object, 'DiscreteClass'))
        stop('Discrete latent structures not yet supported', call.=FALSE)
    if(is.list(object)) inames <- colnames(object[[1]]@Data$data)
    else inames <- colnames(object@Data$data)
    ind <- 1:length(inames)
    if(!is.numeric(item)) item <- ind[inames == item]
    rot <- list(x = rot[[1]], y = rot[[2]], z = rot[[3]])
    if(is.list(object)){
        if(object[[1]]@Model$nfact == 1L) degrees <- 0
    } else if(object@Model$nfact == 1L) degrees <- 0
    rotate <- if(is.null(dots$rotate)) 'none' else dots$rotate
    ret <- itemplot.internal(object=object, item=item, type=type, degrees=degrees, CE=CE,
                             CEalpha=CEalpha, CEdraws=CEdraws, drop.zeros=drop.zeros, rot=rot,
                             theta_lim=theta_lim, par.strip.text=par.strip.text,
                             par.settings=par.settings, auto.key=auto.key, npts=npts, ...)
    if(!is.list(object) && object@Options$exploratory)
        ret$main <- paste0(ret$main, ' (rotate = \'', rotate, '\')')
    return(ret)
}
