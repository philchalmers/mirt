#' Wald statistics for mirt models
#'
#' Compute a Wald test given an \code{L} vector or matrix of numeric contrasts. Requires that the
#' model information matrix be computed (by passing \code{SE = TRUE} when estimating the model). Use
#' \code{wald(model)} to observe how the information matrix columns are named, especially if
#' the estimated model contains constrained parameters (e.g., 1PL).
#'
#'
#' @aliases wald
#' @param L a coefficient matrix with dimensions \code{nconstrasts x npars.estimated},
#'   or a character vector giving the hypothesis in symbolic form
#'   (syntax format borrowed from the \code{car} package; see \code{Details} below).
#'   Omitting this value will return the column names of the
#'   information matrix used to identify the (potentially constrained) parameters
#' @param object estimated object from \code{mirt}, \code{bfactor},
#'   \code{multipleGroup}, \code{mixedmirt}, or \code{mdirt}
#' @param C a constant vector of population parameters to be compared along side L, where
#'   \code{length(C) == row(L)}. By default a vector of 0's is constructed. Note that when using
#'   the syntax input for \code{L} this argument is omitted
#'
#' The following description is borrowed from \code{car} package documentation pertaining to the character vector
#' input to the argument \code{L}: "The hypothesis matrix can be supplied as a numeric matrix (or vector), the rows of which
#' specify linear combinations of the model
#' coefficients, which are tested equal to the corresponding entries in the right-hand-side vector, which defaults to a vector of zeroes.
#'
#' Alternatively, the hypothesis can be specified symbolically as a character vector with one or more elements, each of which gives either
#' a linear combination of coefficients, or a linear equation in the coefficients (i.e., with both a left and right side separated by an
#' equals sign). Components of a linear expression or linear equation can consist of numeric constants, or numeric constants multiplying
#' coefficient names (in which case the number precedes the coefficient, and may be separated from it by spaces or an asterisk);
#' constants of 1 or -1 may be omitted. Spaces are always optional. Components are separated by plus or minus signs. Newlines or tabs
#' in hypotheses will be treated as spaces. See the examples below."
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords wald
#' @export wald
#' @examples
#' \dontrun{
#'
#' #View parnumber index
#' data(LSAT7)
#' data <- expand.table(LSAT7)
#' mod <- mirt(data, 1, SE = TRUE)
#' coef(mod)
#'
#' # see how the information matrix relates to estimated parameters, and how it lines up
#' #   with the parameter index
#' (infonames <- wald(mod))
#' index <- mod2values(mod)
#' index[index$est, ]
#'
#' #second item slope equal to 0?
#' L <- matrix(0, 1, 10)
#' L[1,3] <- 1
#' wald(mod, L)
#'
#' # same as above using character syntax input
#' infonames
#' wald(mod, "a1.5 = 0")
#'
#' #simultaneously test equal factor slopes for item 1 and 2, and 4 and 5
#' L <- matrix(0, 2, 10)
#' L[1,1] <- L[2, 7] <- 1
#' L[1,3] <- L[2, 9] <- -1
#' L
#' wald(mod, L)
#'
#' # Again, using more efficient syntax
#' infonames
#' wald(mod, c("a1.1 = a1.5", "a1.13 = a1.17"))
#'
#' #log-Liklihood tests (requires estimating a new model)
#' cmodel <- 'theta = 1-5
#'            CONSTRAIN = (1,2, a1), (4,5, a1)'
#' mod2 <- mirt(data, cmodel)
#' #or, equivalently
#' #mod2 <- mirt(data, 1, constrain = list(c(1,5), c(13,17)))
#' anova(mod2, mod)
#'
#' #####
#' # test equality of means in multi-group model:
#' #    H0: (mu1 - mu2) = (mu3 - mu4)
#'
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 500
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .5)
#' dataset3 <- simdata(a, d, N, itemtype, mu = -1)
#' dataset4 <- simdata(a, d, N, itemtype, mu = -.5)
#' dat <- rbind(dataset1, dataset2, dataset3, dataset4)
#' group <- factor(rep(paste0('D', 1:4), each=N))
#' levels(group)
#' models <- 'F1 = 1-15'
#'
#' # 3 means estimated
#' mod_free <- multipleGroup(dat, models, group = group, SE=TRUE,
#'                           invariance=c('slopes', 'intercepts', 'free_var','free_means'))
#' wald(mod_free) # obtain parameter names
#' # View(mod2values(mod_free))
#'
#' # reference group mean = 0 by default
#' wald(mod_free, c("0 - MEAN_1.123 = MEAN_1.185 - MEAN_1.247"))
#'
#'
#' }
wald <- function(object, L, C = 0){
    if(missing(object)) missingMsg('object')
    covB <- extract.mirt(object, 'vcov')
    if(all(dim(covB) == c(1,1)))
        if(covB[1,1] == 0L || is.na(covB[1,1]))
            stop('No information matrix has been calculated for the model', call.=FALSE)
    Names <- colnames(covB)
    B <- extract.mirt(object, 'parvec')
    names(B) <- colnames(covB)
    if(missing(L)){
        ret <- round(B, 3)
        names(ret) <- Names
        return(ret)
    }
    if(is.character(L)){
        tmp <- makeHypothesis(names(B), L)
        C <- tmp[, NCOL(tmp)]
        L <- tmp[, -NCOL(tmp), drop = FALSE]
        # rownames(L) <- L
    }
    if(!is.matrix(L)){
        stop('L must be a matrix', call.=FALSE)
    } else if(ncol(L) != length(Names)){
        stop('L does not have an appropriate number of columns', call.=FALSE)
    }
    if(!is.vector(C))
        stop('C must be a vector of constant population parameters', call.=FALSE)
    if(length(C) == 1L){
        C <- numeric(ncol(L))
    } else if(length(C) != ncol(L)){
        stop('length(C) must be the same as ncol(L)', call.=FALSE)
    }
    W <- t(L %*% (B - C)) %*% solve(L %*% covB %*% t(L)) %*% (L %*% (B - C))
    W <- ifelse(W < 0, 0, W)
    ret <- list(W=W, df = nrow(L))
    p <- 1 - pchisq(ret$W, ret$df)
    ret$p <- p
    as.data.frame(ret)
}
