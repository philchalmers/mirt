#' Simulate response patterns
#'
#' Simulates response patterns for compensatory and noncompensatory MIRT models
#' from multivariate normally distributed factor (\eqn{\theta}) scores, or from
#' a user input matrix of \eqn{\theta}'s.
#'
#' Returns a data matrix simulated from the parameters, or a list containing the data,
#' item objects, and Theta matrix.
#'
#' @param a a matrix/vector of slope parameters. If slopes are to be constrained to
#'   zero then use \code{NA} or simply set them equal to 0
#' @param d a matrix/vector of intercepts. The matrix should have as many columns as
#'   the item with the largest number of categories, and filled empty locations
#'   with \code{NA}. When a vector is used the test is assumed to consist only of dichotomous items
#'   (because only one intercept per item is provided). When \code{itemtype = 'lca'} intercepts will not
#'   be used
#' @param itemtype a character vector of length \code{nrow(a)} (or 1, if all the item types are
#'   the same) specifying the type of items to simulate. Inputs can either be the same as
#'   the inputs found in the \code{itemtype} argument in \code{\link{mirt}} or the
#'   internal clases defined by the package. Typical \code{itemtype} inputs that
#'   are passed to \code{\link{mirt}} are used then these will be converted into
#'   the respective internal classes automatically.
#'
#'   If the internal class of the object is specified instead, the inputs can
#'   be \code{'dich', 'graded', 'gpcm','nominal', 'nestlogit', 'partcomp'}, or \code{'lca'}, for
#'   dichotomous, graded, generalized partial credit, nominal, nested logit, partially compensatory,
#'   and latent class analysis model. Note that for the gpcm, nominal, and nested logit models there should
#'   be as many parameters as desired categories, however to parametrized them for meaningful
#'   interpretation the first category intercept should
#'   equal 0 for these models (second column for \code{'nestlogit'}, since first column is for the
#'   correct item traceline). For nested logit models the 'correct' category is always the lowest
#'   category (i.e., == 1). It may be helpful to use \code{\link{mod2values}} on data-sets that
#'   have already been estimated to understand the itemtypes more intimately
#' @param nominal a matrix of specific item category slopes for nominal models.
#'   Should be the dimensions as the intercept specification with one less column, with \code{NA}
#'   in locations where not applicable. Note that during estimation the first slope will be
#'   constrained to 0 and the last will be constrained to the number of categories minus 1,
#'   so it is best to set these as the values for the first and last categories as well
#' @param N sample size
#' @param guess a vector of guessing parameters for each item; only applicable
#'   for dichotomous items. Must be either a scalar value that will affect all of
#'   the dichotomous items, or a vector with as many values as to be simulated items
#' @param upper same as \code{guess}, but for upper bound parameters
#' @param gpcm_mats a list of matricies specifying the scoring scheme for generalized partial
#'   credit models (see \code{\link{mirt}} for details)
#' @param sigma a covariance matrix of the underlying distribution. Default is
#'   the identity matrix. Used when \code{Theta} is not supplied
#' @param mu a mean vector of the underlying distribution. Default is a vector
#'   of zeros. Used when \code{Theta} is not supplied
#' @param Theta a user specified matrix of the underlying ability parameters,
#'   where \code{nrow(Theta) == N} and \code{ncol(Theta) == ncol(a)}. When this is supplied the
#'   \code{N} input is not required
#' @param returnList logical; return a list containing the data, item objects defined
#'   by \code{mirt} containing the population parameters and item structure, and the
#'   latent trait matrix \code{Theta}? Default is FALSE
#' @param model a single group object, typically returned by functions such as \code{\link{mirt}} or
#'   \code{\link{bfactor}}. Supplying this will render all other parameter elements (excluding the
#'   \code{Theta}, \code{N}, \code{mu}, and \code{sigma} inputs) redundent (unless explicitly provided)
#' @param which.items an integer vector used to indicate which items to simulate when a
#'   \code{model} input is included. Default simulates all items
#' @param mins an integer vector (or single value to be used for each item) indicating what
#'   the lowest category should be. If \code{model} is supplied then this will be extracted from
#'   \code{slot(mod, 'Data')$mins}, otherwise the default is 0
#' @param lca_cats a vector indicating how many categories each lca item should have. If not supplied
#'   then it is assumed that 2 categories should be generated for each item
#' @param prob.list an optional list containing matrix/data.frames of probabilities values for
#'   each category to be simulated. This is useful when creating customized probability functions
#'   to be sampled from
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Reckase, M. D. (2009). \emph{Multidimensional Item Response Theory}. New York: Springer.
#' @keywords data
#' @export simdata
#' @examples
#'
#' \dontrun{
#' ### Parameters from Reckase (2009), p. 153
#'
#' set.seed(1234)
#'
#' a <- matrix(c(
#'  .7471, .0250, .1428,
#'  .4595, .0097, .0692,
#'  .8613, .0067, .4040,
#' 1.0141, .0080, .0470,
#'  .5521, .0204, .1482,
#' 1.3547, .0064, .5362,
#' 1.3761, .0861, .4676,
#'  .8525, .0383, .2574,
#' 1.0113, .0055, .2024,
#'  .9212, .0119, .3044,
#'  .0026, .0119, .8036,
#'  .0008, .1905,1.1945,
#'  .0575, .0853, .7077,
#'  .0182, .3307,2.1414,
#'  .0256, .0478, .8551,
#'  .0246, .1496, .9348,
#'  .0262, .2872,1.3561,
#'  .0038, .2229, .8993,
#'  .0039, .4720, .7318,
#'  .0068, .0949, .6416,
#'  .3073, .9704, .0031,
#'  .1819, .4980, .0020,
#'  .4115,1.1136, .2008,
#'  .1536,1.7251, .0345,
#'  .1530, .6688, .0020,
#'  .2890,1.2419, .0220,
#'  .1341,1.4882, .0050,
#'  .0524, .4754, .0012,
#'  .2139, .4612, .0063,
#'  .1761,1.1200, .0870),30,3,byrow=TRUE)*1.702
#'
#' d <- matrix(c(.1826,-.1924,-.4656,-.4336,-.4428,-.5845,-1.0403,
#'   .6431,.0122,.0912,.8082,-.1867,.4533,-1.8398,.4139,
#'   -.3004,-.1824,.5125,1.1342,.0230,.6172,-.1955,-.3668,
#'   -1.7590,-.2434,.4925,-.3410,.2896,.006,.0329),ncol=1)*1.702
#'
#' mu <- c(-.4, -.7, .1)
#' sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
#'
#' dataset1 <- simdata(a, d, 2000, itemtype = '2PL')
#' dataset2 <- simdata(a, d, 2000, itemtype = '2PL', mu = mu, sigma = sigma)
#'
#' #mod <- mirt(dataset1, 3, method = 'MHRM')
#' #coef(mod)
#'
#' ### Unidimensional graded response model with 5 categories each
#'
#' a <- matrix(rlnorm(20,.2,.3))
#'
#' # for the graded model, ensure that there is enough space between the intercepts,
#' # otherwise closer categories will not be selected often (minimum distance of 0.3 here)
#' diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
#' diffs <- -(diffs - rowMeans(diffs))
#' d <- diffs + rnorm(20)
#'
#' dat <- simdata(a, d, 500, itemtype = 'graded')
#' # mod <- mirt(dat, 1)
#'
#' ### An example of a mixed item, bifactor loadings pattern with correlated specific factors
#'
#' a <- matrix(c(
#' .8,.4,NA,
#' .4,.4,NA,
#' .7,.4,NA,
#' .8,NA,.4,
#' .4,NA,.4,
#' .7,NA,.4),ncol=3,byrow=TRUE)
#'
#' d <- matrix(c(
#' -1.0,NA,NA,
#'  1.5,NA,NA,
#'  0.0,NA,NA,
#' 0.0,-1.0,1.5,  #the first 0 here is the recommended constraint for nominal
#' 0.0,1.0,-1, #the first 0 here is the recommended constraint for gpcm
#' 2.0,0.0,NA),ncol=3,byrow=TRUE)
#'
#' nominal <- matrix(NA, nrow(d), ncol(d))
#' #the first 0 and last (ncat - 1) = 2 values are the recommended constraints
#' nominal[4, ] <- c(0,1.2,2)
#'
#' sigma <- diag(3)
#' sigma[2,3] <- sigma[3,2] <- .25
#' items <- c('2PL','2PL','2PL','nominal','gpcm','graded')
#'
#' dataset <- simdata(a,d,2000,items,sigma=sigma,nominal=nominal)
#'
#' #mod <- bfactor(dataset, c(1,1,1,2,2,2), itemtype=c(rep('2PL', 3), 'nominal', 'gpcm','graded'))
#' #coef(mod)
#'
#' #### Convert standardized factor loadings to slopes
#'
#' F2a <- function(F, D=1.702){
#'     h2 <- rowSums(F^2)
#'     a <- (F / sqrt(1 - h2)) * D
#'     a
#' }
#'
#' (F <- matrix(c(rep(.7, 5), rep(.5,5))))
#' (a <- F2a(F))
#' d <- rnorm(10)
#'
#' dat <- simdata(a, d, 5000, itemtype = '2PL')
#' mod <- mirt(dat, 1)
#' coef(mod, simplify=TRUE)$items
#' summary(mod)
#'
#' mod2 <- mirt(dat, 'F1 = 1-10
#'                    CONSTRAIN = (1-5, a1), (6-10, a1)')
#' summary(mod2)
#' anova(mod, mod2)
#'
#' #### Unidimensional nonlinear factor pattern
#'
#' theta <- rnorm(2000)
#' Theta <- cbind(theta,theta^2)
#'
#' a <- matrix(c(
#' .8,.4,
#' .4,.4,
#' .7,.4,
#' .8,NA,
#' .4,NA,
#' .7,NA),ncol=2,byrow=TRUE)
#' d <- matrix(rnorm(6))
#' itemtype <- rep('2PL',6)
#'
#' nonlindata <- simdata(a=a, d=d, itemtype=itemtype, Theta=Theta)
#'
#' #model <- '
#' #F1 = 1-6
#' #(F1 * F1) = 1-3'
#' #mod <- mirt(nonlindata, model)
#' #coef(mod)
#'
#' #### 2PLNRM model for item 4 (with 4 categories), 2PL otherwise
#'
#' a <- matrix(rlnorm(4,0,.2))
#'
#' #first column of item 4 is the intercept for the correct category of 2PL model,
#' #    otherwise nominal model configuration
#' d <- matrix(c(
#' -1.0,NA,NA,NA,
#'  1.5,NA,NA,NA,
#'  0.0,NA,NA,NA,
#'  1, 0.0,-0.5,0.5),ncol=4,byrow=TRUE)
#'
#' nominal <- matrix(NA, nrow(d), ncol(d))
#' nominal[4, ] <- c(NA,0,.5,.6)
#'
#' items <- c(rep('2PL',3),'nestlogit')
#'
#' dataset <- simdata(a,d,2000,items,nominal=nominal)
#'
#' #mod <- mirt(dataset, 1, itemtype = c('2PL', '2PL', '2PL', '2PLNRM'), key=c(NA,NA,NA,1))
#' #coef(mod)
#' #itemplot(mod,4)
#'
#' #return list of simulation parameters
#' listobj <- simdata(a,d,2000,items,nominal=nominal, returnList=TRUE)
#' str(listobj)
#'
#' # generate dataset from converged model
#' mod <- mirt(Science, 1, itemtype = c(rep('gpcm', 3), 'nominal'))
#' sim <- simdata(model=mod, N=1000)
#' head(sim)
#'
#' Theta <- matrix(rnorm(100))
#' sim <- simdata(model=mod, Theta=Theta)
#' head(sim)
#'
#' # alternatively, define a suitable object with functions from the mirtCAT package
#' # help(generate.mirt_object)
#' library(mirtCAT)
#'
#' nitems <- 50
#' a1 <- rlnorm(nitems, .2,.2)
#' d <- rnorm(nitems)
#' g <- rbeta(nitems, 20, 80)
#' pars <- data.frame(a1=a1, d=d, g=g)
#' head(pars)
#'
#' obj <- generate.mirt_object(pars, '3PL')
#' dat <- simdata(N=200, model=obj)
#'
#' ######
#' # prob.list example
#'
#' # custom probabilty function that returns a matrix
#' fun <- function(a, b, theta){
#'     P <- 1 / (1 + exp(-a * (theta-b)))
#'     cbind(1-P, P)
#' }
#'
#' set.seed(1)
#' theta <- matrix(rnorm(100))
#' prob.list <- list()
#' nitems <- 5
#' a <- rlnorm(nitems, .2, .2); b <- rnorm(nitems, 0, 1/2)
#' for(i in 1:nitems) prob.list[[i]] <- fun(a[i], b[i], theta)
#' str(prob.list)
#'
#' dat <- simdata(prob.list=prob.list)
#' head(dat)
#'
#' }
simdata <- function(a, d, N, itemtype, sigma = NULL, mu = NULL, guess = 0,
	upper = 1, nominal = NULL, Theta = NULL, gpcm_mats = list(), returnList = FALSE,
	model = NULL, which.items = NULL, mins = 0, lca_cats = NULL, prob.list = NULL)
{
    if(!is.null(prob.list)){
        if(!all(sapply(prob.list, function(x) is.matrix(x) || is.data.frame(x))))
            stop('Elements of prob.list must be either a matrix or data.frame')
        prob.list <- lapply(prob.list, as.matrix)
        if(!all(sapply(prob.list, nrow) == nrow(prob.list[[1L]])))
            stop('prob.list elements have unequal rows')
        K <- sapply(prob.list, ncol)
        nitems <- length(K)
        if(any(K == 1L)) stop('prob.list elements should have more than 1 column')
        if(length(mins) == 1L) mins <- rep(mins, nitems)
        stopifnot(length(mins) == nitems)
        data <- matrix(NA, nrow(prob.list[[1L]]), nitems)
        for(i in 1L:nitems) data[,i] <- respSample(prob.list[[i]])
        data <- (t(t(data) + mins))
        colnames(data) <- paste("Item_", 1L:nitems, sep="")
        return(data)
    }
    if(!is.null(model)){
        stopifnot(is(model, 'SingleGroupClass'))
        nitems <- extract.mirt(model, 'nitems')
        if(is.null(which.items)) which.items <- 1L:nitems
        nfact <- extract.mirt(model, 'nfact')
        cfs <- coef(model, simplify=TRUE)
        if(is.null(sigma)) sigma <- cfs$cov
        if(is.null(mu)) mu <- cfs$means
        if(is.null(Theta)){
            if(missing(N)) N <- nrow(extract.mirt(model, 'data'))
            Theta <- mirt_rmvnorm(N,mu,sigma,check=TRUE)
        } else N <- nrow(Theta)
        data <- matrix(0, N, nitems)
        colnames(data) <- extract.mirt(model, 'itemnames')
        for(i in which.items){
            obj <- extract.item(model, i)
            P <- ProbTrace(obj, Theta)
            data[,i] <- respSample(P)
        }
        ret <- t(t(data) + model@Data$mins)
        return(ret[,which.items, drop=FALSE])
    }
    if(missing(N) && is.null(Theta)) missingMsg('N or Theta')
    if(missing(a)) missingMsg('a')
    if(missing(itemtype)) missingMsg('itemtype')
    if(missing(d) && !all(itemtype == 'lca')) missingMsg('d')
    if(is.vector(a)) a <- matrix(a)
    if(missing(d)) d <- matrix(1, nrow(a))
    if(is.vector(d)) d <- matrix(d)
    if(any(itemtype == 'nominal') && is.null(nominal))
        stop('nominal itemtypes require a \'nominal\' matrix input of scoring coefs (the ak values)',
             call.=FALSE)
	nfact <- ncol(a)
	nitems <- nrow(a)
	if(length(mins) == 1L) mins <- rep(mins, nitems)
	stopifnot(length(mins) == nitems)
	K <- rep(0L,nitems)
	if(length(guess) == 1L) guess <- rep(guess,nitems)
	if(length(guess) != nitems) stop("Guessing parameter is incorrect", call.=FALSE)
	if(length(upper) == 1L) upper <- rep(upper,nitems)
	if(length(upper) != nitems) stop("Upper bound parameter is incorrect", call.=FALSE)
    if(length(itemtype) == 1L) itemtype <- rep(itemtype, nitems)
    if(length(gpcm_mats)){
        stopifnot(length(gpcm_mats) == nitems)
        use_gpcm_mats <- sapply(gpcm_mats, is.matrix)
    } else use_gpcm_mats <- rep(FALSE, nitems)
    if(any(itemtype %in% Valid_iteminputs())){
        if(any(itemtype %in% c('grsm', 'grsmIRT')))
            stop('Please rewrite rating scale models as gpcm', call.=FALSE)
        if(any(itemtype %in% c('Rasch')))
            stop('Rasch itemtype is ambiguous, please specifiy either gpcm or dich/2PL class', call.=FALSE)
        itemtype <- ifelse(itemtype %in% c('2PL', '3PL', '3PLu', '4PL'), 'dich', itemtype)
        itemtype <- ifelse(itemtype %in% c('PC2PL', 'PC3PL'), 'partcomp', itemtype)
        itemtype <- ifelse(itemtype %in% c("2PLNRM", "3PLNRM", "3PLuNRM", "4PLNRM"), 'nestlogit', itemtype)
    }
	for(i in 1L:length(K)){
	    K[i] <- length(na.omit(d[i, ])) + 1L
	    if(itemtype[i] =='partcomp') K[i] <- 2L
	    if(any(itemtype[i] == c('gpcm', 'nominal', 'nestlogit'))) K[i] <- K[i] - 1L
	}
	if(!is.null(lca_cats)) K[itemtype == 'lca'] <- lca_cats[itemtype == 'lca']
	K <- as.integer(K)
    if(any(guess > 1 | guess < 0)) stop('guess input must be between 0 and 1', call.=FALSE)
    if(any(upper > 1 | upper < 0)) stop('upper input must be between 0 and 1', call.=FALSE)
    guess <- logit(guess)
    upper <- logit(upper)
    oldguess <- guess
    oldupper <- upper
    guess[K > 2L] <- upper[K > 2L] <- NA
    guess[itemtype == 'nestlogit'] <- oldguess[itemtype == 'nestlogit']
    upper[itemtype == 'nestlogit'] <- oldupper[itemtype == 'nestlogit']
	if(is.null(sigma)) sigma <- diag(nfact)
	if(is.null(mu)) mu <- rep(0,nfact)
	if(!is.null(Theta)){
		if(ncol(Theta) != nfact)
			stop("The input Theta matrix does not have the correct dimensions", call.=FALSE)
	    N <- nrow(Theta)
	}
	if(is.null(Theta)) Theta <- mirt_rmvnorm(N,mu,sigma,check=TRUE)
    if(is.null(nominal)) nominal <- matrix(NA, nitems, max(K))
	data <- matrix(0, N, nitems)
    a[is.na(a)] <- 0
    itemobjects <- vector('list', nitems)
	for(i in 1L:nitems){
	    if(itemtype[i] == 'nestlogit'){
	        par <- na.omit(c(a[i, ],d[i,1], guess[i], upper[i], nominal[i,-1L],d[i,-1L]))
	        obj <- new(itemtype[i], par=par, nfact=nfact, correctcat=1L)
	    } else {
            if(itemtype[i] == 'gpcm'){
                if(!use_gpcm_mats[i]){
                    par <- na.omit(c(a[i, ],0:(K[i]-1), d[i,]))
                } else {
                    stopifnot(nrow(gpcm_mats[[i]]) == K[i])
                    stopifnot(ncol(gpcm_mats[[i]]) == nfact)
                    par <- na.omit(c(a[i, ],as.vector(gpcm_mats[[i]]), d[i,]))
                }
            } else if(itemtype[i] == 'ideal'){
                if(K[i] > 2) stop('ideal point models for dichotomous items only', call.=FALSE)
                if(d[i,1] > 0) stop('ideal point intercepts must be negative', call.=FALSE)
                par <- na.omit(c(a[i, ],d[i,]))
            } else if(itemtype[i] == 'lca'){
                par <- na.omit(a[i, ])
            } else {
                par <- na.omit(c(a[i, ],nominal[i,],d[i,],guess[i],upper[i]))
            }
            obj <- new(itemtype[i], par=par, nfact=nfact, ncat=K[i])
            if(itemtype[i] %in% c('gpcm', 'nominal')) obj@mat <- FALSE
            if(use_gpcm_mats[i]) obj@mat <- TRUE
	    }
        if(any(itemtype[i] == c('gpcm','nominal', 'nestlogit')))
            obj@ncat <- K[i]
        P <- ProbTrace(obj, Theta)
        data[,i] <- respSample(P)
        itemobjects[[i]] <- obj
	}
    data <- (t(t(data) + mins))
	colnames(data) <- paste("Item_", 1L:nitems, sep="")
    if(returnList){
        return(list(itemobjects=itemobjects, data=data, Theta=Theta))
    } else {
	    return(data)
    }
}

