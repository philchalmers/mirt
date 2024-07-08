#' Reverse score one or more items from a response matrix
#'
#' Reverse score specific items given empirical range or specific scoring
#' range.
#'
#' @param data an object of class \code{data.frame}, \code{matrix}, or
#'   \code{table} with the response patterns
#' @param which names of items in \code{data} that should be rescored
#' @param range (optional) a named \code{list} to specify the low and high
#'   score ranges. Specified names must match the names found in
#'   \code{data}, and each element of this list should contain only two values.
#'   If items specified in \code{which} are omitted
#'   then the empirical min/max information will be used instead
#'
#' @return returns the original \code{data} object with the specified
#'   items reverse scored replacing the original scoring scheme
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @export
#' @examples
#'
#' a <- rlnorm(20)
#' a[c(1,5,10)] <- -a[c(1,5,10)]
#' diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
#' diffs <- -(diffs - rowMeans(diffs))
#' d <- diffs + rnorm(20)
#' dat <- simdata(a,d,itemtype='graded', N=300)
#' head(dat)
#'
#' \dontrun{
#' # fitted model has negative slopes due to flipped scoring
#' mod <- mirt(dat)
#' coef(mod, simplify=TRUE)$items
#' plot(mod, type = 'itemscore')
#' }
#'
#' # reverse the scoring for items 1, 5, and 10 only using empirical min/max
#' revdat <- reverse.score(dat, c('Item_1', 'Item_5', 'Item_10'))
#' head(revdat)
#'
#' head(dat[,c(1,5,10)]) # compare
#' head(revdat[,c(1,5,10)])
#'
#' \dontrun{
#' # slopes all positive now
#' mod2 <- mirt(revdat)
#' coef(mod2, simplify=TRUE)$items
#' plot(mod2, type = 'itemscore')
#' }
#'
#' # use different empirical scoring information due to options not used
#'   # 0 score not observed for item 1, though should have been rescored to a 4
#' dat[dat[,1] == 0, 1] <- 1
#' table(dat[,1])
#'
#' # specify theoretical scoring values in the range list
#' revdat2 <- reverse.score(dat, c('Item_1', 'Item_5', 'Item_10'),
#'                               range = list(Item_1 = c(0,4)))
#' head(revdat2[,1:3])
#' table(revdat2[,1])
#'
#'
reverse.score <- function (data, which, range = NULL){
    if(missing(data)) missingMsg('data')
    if(length(colnames(data)) != ncol(data))
        stop('data must contain suitable column names')
    subdat <- data[,which, drop=FALSE]
    min <- apply(subdat, 2L, min, na.rm=TRUE)
    max <- apply(subdat, 2L, max, na.rm=TRUE)
    if(!is.null(range)){
        stopifnot(length(names(range))>0)
        stopifnot("range contains names that do not match data" =
                      all(names(range) %in% colnames(subdat)))
        for(i in length(range)){
            pick <- which(names(range)[i] == colnames(subdat))
            min[pick] <- min(range[[i]])
            max[pick] <- max(range[[i]])
        }
    }
    subdat <- t((max - min) - t(subdat))
    data[,which] <- subdat
    data
}
