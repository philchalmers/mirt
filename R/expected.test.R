#' Function to calculate expected test score
#'
#' Given an estimated model compute the expected test score. Returns the expected values in the 
#' same form as the data used to estimate the model.
#'
#' @aliases expected.test
#' @param x an estimated mirt object
#' @param Theta a matrix of latent trait values
#' @param group a number signifying which group the item should be extracted from (applies to
#'   'MultipleGroupClass' objects only)
#' @keywords expected score
#' @seealso \code{\link{extract.item}}, \code{\link{expected.item}}
#' @export expected.test
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(deAyala)
#' model <- mirt.model('F = 1-5
#'                     CONSTRAIN = (1-5, a1)')
#' mod <- mirt(dat, model)
#'
#' Theta <- matrix(seq(-6,6,.01))
#' tscore <- expected.test(mod, Theta)
#' tail(cbind(Theta, tscore))
#' 
#' }
expected.test <- function(x, Theta, group = NULL){
    J <- ncol(x@Data$data)
    score <- numeric(nrow(Theta))
    mins <- apply(x@Data$data, 2L, min, na.rm=TRUE)
    for(i in 1L:J){
        item <- extract.item(x, i, group=group)
        score <- score + expected.item(item, Theta=Theta, min=0L) + mins[i]
    }
    return(score)
}