#' Extract a group from a multiple group mirt object
#'
#' Extract a single group from an object defined by \code{\link{multipleGroup}}.
#'
#' @aliases extract.group
#' @param x mirt model of class 'MultipleGroupClass'
#' @param group a number signifying which group should be extracted
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords extract
#' @seealso \code{\link{extract.item}}, \code{\link{extract.mirt}}
#' @export extract.group
#' @examples
#'
#' \dontrun{
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#' models <- 'F1 = 1-15'
#'
#' mod_configural <- multipleGroup(dat, models, group = group)
#' group.1 <- extract.group(mod_configural, 1) #extract first group
#' summary(group.1)
#' plot(group.1)
#' }
extract.group <- function(x, group){
    if(missing(x)) missingMsg('x')
    if(missing(group)) missingMsg('group')
    if(!is(x, 'MultipleGroupClass'))
        stop('Model was not estimated with multipleGroup()', call.=FALSE)
    if(missing(group)) stop('Must specify group number', call.=FALSE)
    vals <- mod2values(x)
    vals <- vals[vals$group == x@Data$groupNames[group], ]
    dat <- extract.mirt(x, 'data')
    nfact <- extract.mirt(x, 'nfact')
    K <- extract.mirt(x, 'K')
    groupvec <- extract.mirt(x, 'group')
    groupNames <- extract.mirt(x, 'groupNames')
    sv <- mirt(dat[groupvec == groupNames[group], ], nfact,
                pars = 'values', technical = list(customK = K))
    sv$value <- vals$value
    mod <- mirt(dat[groupvec == groupNames[group], ], nfact,
                pars = sv, technical = list(customK = K, warn=FALSE), TOL = NaN, quadpts=1L)
    return(mod)
}
