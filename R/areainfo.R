#' Function to calculate the area under a selection of information curves
#'
#' Compute the area within test or item information over a definite integral range.
#'
#' @aliases areainfo
#' @param x an estimated mirt object
#' @param theta_lim range of integration to be computed
#' @param which.items an integer vector indicating which items to include in the expected information function.
#'   Default uses all possible items
#'
#' @return a \code{data.frame} with the lower and upper integration range, the information area
#'   within the range (Info), the information area over the range -10 to 10 (Total.Info), proportion
#'   of total information given the integration range (Info.Proportion), and the number of items included (nitems)
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords information area
#' @export areainfo
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' mod <- mirt(dat, 1)
#'
#' areainfo(mod, c(-2,0), which.items = 1) #item 1
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
areainfo <- function(x, theta_lim, which.items = NULL, ...){
    f <- function(theta, mod, which.items)
        testinfo(x=mod, Theta=matrix(theta), which.items=which.items)
    if(missing(x)) missingMsg('x')
    if(missing(theta_lim)) missingMsg('theta_lim')
    stopifnot(length(theta_lim) == 2L)
    stopifnot(extract.mirt(x, 'nfact') == 1L)

    nitems <- ifelse(is.null(which.items), extract.mirt(x, 'nitems'), length(which.items))
    ii <- integrate(f, lower = theta_lim[1L], upper = theta_lim[2L],
                    mod=mod, which.items=which.items, ...)
    iT <- integrate(f, lower = -10, upper = 10, mod=mod, which.items=which.items, ...)
    ret <- data.frame(LowerBound=min(theta_lim), UpperBound=max(theta_lim),
                      Info=ii$value, TotalInfo=iT$value, Proportion=ii$value/iT$value,
                      nitems=nitems)
    rownames(ret) <- ''
    ret
}
