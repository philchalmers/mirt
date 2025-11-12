#' Function to calculate the marginal moments for items and bundles fit via IRT models
#'
#' Given an estimated model and a prior density function, compute the marginal moments for either
#' and item or a bundle of items. Function returns the first found moments implied by the model
#' and select density function (MEAN, VAR, SKEW, and KURT). Currently limited to unidimensional IRT models.
#'
#' @param mod an object of class \code{'SingleGroupClass'} or \code{'MultipleGroupClass'}
#' @param which.items vector indicating which items to use in the computation of the expected values.
#'   Default (\code{NULL}) uses all available items
#' @param bundle logical; given \code{which.items}, should the composite of the response functions be used
#'  as a collective bundle? If \code{TRUE} then \code{\link{expected.test}} will be used, otherwise
#'  \code{\link{expected.item}} will be used
#' @param group optional indicator to return only specific group information for multiple group models.
#'   Default compute moments for each group, returning a \code{list}
#' @param density a density function to use for integration. Default assumes the latent traits are from a
#'   normal (Gaussian) distribution. Function definition must be of the form \code{function(quadrature, mean, sd)}
#'   as the values of the mean/variance are extracted and passed from the supplied model
#' @param Theta_lim range of integration grid to use when forming expectations
#' @param quadpts number of discrete quadrature to use in the computations
#' @param ... additional arguments passed to the density function
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' @export
#' @seealso \code{\link{expected.item}}, \code{\link{expected.test}}
#' @examples
#'
#' \donttest{
#'
#' # single group
#' dat <- expand.table(deAyala)
#' mod <- mirt(dat)
#' TS <- rowSums(dat)
#'
#' # expected moments of total scores given model
#' marginal_moments(mod)
#' c(mean=mean(TS), var=var(TS))  # first two moments of data
#'
#' # same, however focusing on individual items
#' marginal_moments(mod, bundle = FALSE)
#' cbind(mean=colMeans(dat), var=apply(dat, 2, var)) # first two moments of data
#'
#' ############################################
#' ## same as above, however with multiple group model
#'
#' set.seed(1234)
#' group <- sample(c('G1', 'G2'), nrow(dat), replace=TRUE)
#' modMG <- multipleGroup(dat, group=group,
#'  invariance=c(colnames(dat), 'free_mean', 'free_var'))
#' coef(modMG, simplify=TRUE)
#'
#' # expected moments of total scores given model
#' marginal_moments(modMG)
#' marginal_moments(modMG, group = 'G1') # specific group only
#'
#' # same, however focusing on individual items
#' marginal_moments(modMG, bundle = FALSE)
#'
#'
#' }
marginal_moments <- function(mod, which.items = NULL, group = NULL, bundle = TRUE,
                             density = NULL, Theta_lim = c(-6,6), quadpts = 121, ...){
    stopifnot(extract.mirt(mod, 'nfact') == 1L)
    if(!is.null(group))
        mod <- extract.group(mod, group=group)
    stopifnot(extract.mirt(mod, 'nfact') == 1)        ## TODO: can make into a grid
    if(is(mod, 'MultipleGroupClass') && is.null(group)){
        grp <- extract.mirt(mod, 'groupNames')
        out <- lapply(grp, \(g)
              marginal_moments(mod=mod, which.items=which.items, group=g,
                               bundle=bundle, density=density, Theta_lim=Theta_lim,
                               quadpts=quadpts, ...))
        names(out) <- grp
        return(out)
    }
    stopifnot(is(mod, 'SingleGroupClass'))
    if(is.null(which.items))
        which.items <- 1:extract.mirt(mod, 'nitems')
    Theta <- matrix(seq(Theta_lim[1], Theta_lim[2], length.out=quadpts))
    if(is.null(density))
        density <- function(theta, mean, sd)
            dnorm(theta, mean=mean, sd=sd)
    cfs <- coef(mod, simplify=TRUE)
    den <- density(Theta, mean=cfs$means[1], sd=sqrt(cfs$cov[1,1]))
    den <- den / sum(den)
    if(bundle){
        Eitem <- expected.test(mod, Theta=Theta, which.items=which.items)
        E <- sum(Eitem * den)
        VAR <- sum((Eitem - E)^2 * den)
        SKEW <- sum((Eitem - E)^3 * den) / VAR^(3/2)
        KURT <- sum((Eitem - E)^4 * den) / VAR^2
        ret <- data.frame(MEAN=E, VAR=VAR, SKEW=SKEW, KURT=KURT)
    } else {
        ret <- as.data.frame(matrix(NA, ncol=4, nrow=length(which.items)))
        rownames(ret) <- extract.mirt(mod, 'itemnames')[which.items]
        for(i in 1:length(which.items)){
            pick <- which.items[i]
            ii <- extract.item(mod, item=pick)
            Eitem <- expected.item(ii, Theta=Theta)
            E <- sum(Eitem * den)
            VAR <- sum((Eitem - E)^2 * den)
            SKEW <- sum((Eitem - E)^3 * den) / VAR^(3/2)
            KURT <- sum((Eitem - E)^4 * den) / VAR^2
            ret[i, ] <- c(E, VAR, SKEW, KURT)
        }
    }
    ret <- as.mirt_df(ret)
    colnames(ret) <- c('MEAN', 'VAR', 'SKEW', 'KURT')
    ret
}
