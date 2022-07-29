#' Function to calculate the area under a selection of information curves
#'
#' Compute the area of a test or item information function over a definite integral range.
#'
#' @aliases areainfo
#' @param x an object of class 'SingleGroupClass', or an object of class 'MultipleGroupClass' if a suitable
#'   \code{group} input were supplied
#' @param theta_lim range of integration to be computed
#' @param which.items an integer vector indicating which items to include in the expected information function.
#'   Default uses all possible items
#' @param group group argument to pass to \code{\link{extract.group}} function. Required when the input object is
#'   a multiple-group model
#' @param ... additional arguments passed to \code{\link{integrate}}
#'
#' @return a \code{data.frame} with the lower and upper integration range, the information area
#'   within the range (Info), the information area over the range -10 to 10 (Total.Info), proportion
#'   of total information given the integration range (Info.Proportion), and the number of items included (nitems)
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords information area
#' @export areainfo
#' @examples
#'
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#'
#' areainfo(mod, c(-2,0), which.items = 1) #item 1
#' \dontrun{
#' areainfo(mod, c(-2,0), which.items = 1:3) #items 1 to 3
#' areainfo(mod, c(-2,0)) # all items (total test information)
#'
#' # plot the area
#' area <- areainfo(mod, c(-2,0))
#' Theta <- matrix(seq(-3,3, length.out=1000))
#' info <- testinfo(mod, Theta)
#' plot(info ~ Theta, type = 'l')
#'
#' pick <- Theta >= -2 & Theta <=0
#' polygon(c(-2, Theta[pick], 0), c(0, info[pick], 0), col='lightblue')
#' text(x = 2, y = 0.5, labels = paste("Total Information:", round(area$TotalInfo, 3),
#'            "\n\nInformation in (-2, 0):", round(area$Info, 3),
#'            paste("(", round(100 * area$Proportion, 2), "%)", sep = "")), cex = 1.2)
#'
#' }
areainfo <- function(x, theta_lim, which.items = 1:extract.mirt(x, 'nitems'), group = NULL, ...){
    if(is(x, 'MultipleGroupClass') && is.null(group))
        stop('Input must be a SingleGroupClass object or a MultipleGroupClass object with a suitable group input',
             call.=FALSE)
    if(!is.null(group) && is(x, 'MultipleGroupClass'))
        x <- extract.group(x=x, group=group)
    f <- function(theta, x, which.items)
        testinfo(x=x, Theta=matrix(theta), which.items=which.items)
    if(missing(x)) missingMsg('x')
    if(missing(theta_lim)) missingMsg('theta_lim')
    stopifnot(is(x, 'SingleGroupClass'))
    stopifnot(length(theta_lim) == 2L)
    stopifnot(extract.mirt(x, 'nfact') == 1L)

    nitems <- ifelse(is.null(which.items), extract.mirt(x, 'nitems'), length(which.items))
    ii <- integrate(f, lower = theta_lim[1L], upper = theta_lim[2L],
                    x=x, which.items=which.items, ...)
    iT <- integrate(f, lower = -Inf, upper = Inf, x=x, which.items=which.items, ...)
    ret <- data.frame(LowerBound=min(theta_lim), UpperBound=max(theta_lim),
                      Info=ii$value, TotalInfo=iT$value, Proportion=ii$value/iT$value,
                      nitems=nitems)
    rownames(ret) <- ''
    ret
}
