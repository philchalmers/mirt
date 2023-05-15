#' Differential Response Functioning statistics
#'
#' Function performs various omnibus differential item (DIF), bundle (DBF), and test (DTF)
#' functioning procedures on an object
#' estimated with \code{multipleGroup()}. The compensatory and non-compensatory statistics provided
#' are described in Chalmers (2018), which generally can be interpreted as IRT generalizations
#' of the SIBTEST and CSIBTEST statistics. For hypothesis tests, these measures
#' require the ACOV matrix to be computed in the
#' fitted multiple-group model (otherwise, sets of plausible draws from the posterior are explicitly
#' required).
#'
#' The effect sizes estimates by the DRF function are
#' \deqn{sDRF = \int [S(C|\bm{\Psi}^{(R)},\theta) S(C|\bm{\Psi}^{(F)},\theta)] f(\theta)d\theta,}
#' \deqn{uDRF = \int |S(C|\bm{\Psi}^{(R)},\theta) S(C|\bm{\Psi}^{(F)},\theta)| f(\theta)d\theta,}
#' and
#' \deqn{dDRF = \sqrt{\int [S(C|\bm{\Psi}^{(R)},\theta) S(C|\bm{\Psi}^{(F)},\theta)]^2 f(\theta)d\theta}}
#' where \eqn{S(.)} are the scoring equations used to evaluate the model-implied
#' difference between the focal and reference group. The \eqn{f(\theta)}
#' terms can either be estimated from the posterior via an empirical
#' histogram approach (default), or can use the best
#' fitting prior distribution that is obtain post-convergence (default is a Guassian
#' distribution). Note that, in comparison to Chalmers (2018), the focal group is
#' the leftmost scoring function while the reference group is the rightmost
#' scoring function. This is largely to keep consistent with similar effect
#' size statistics, such as SIBTEST, DFIT, Wainer's measures of impact, etc,
#' which in general can be seen as special-case estimators of this family.
#'
#' @aliases DRF
#' @param mod a multipleGroup object which estimated only 2 groups
#' @param draws a number indicating how many draws to take to form a suitable multiple imputation
#'   or bootstrap estimate of the expected test scores (100 or more). If \code{boot = FALSE},
#'   requires an estimated parameter information matrix. Returns a list containing the
#'   bootstrap/imputation distribution and null hypothesis test for the sDRF statistics
#' @param focal_items a numeric vector indicating which items to include in the DRF tests. The
#'   default uses all of the items (note that including anchors in the focal items has no effect
#'   because they are exactly equal across groups). Selecting fewer items will result in tests of
#'   'differential bundle functioning'
#' @param param_set an N x p matrix of parameter values drawn from the posterior (e.g., using the
#'   parametric sampling approach, bootstrap, of MCMC). If supplied, then these will be used to compute
#'   the DRF measures. Can be much more efficient to pre-compute these values if DIF, DBF, or DTF are
#'   being evaluated within the same model (especially when using the bootstrap method).
#'   See \code{\link{draw_parameters}}
#' @param CI range of confidence interval when using draws input
#' @param npts number of points to use for plotting. Default is 1000
#' @param quadpts number of quadrature nodes to use when constructing DRF statistics. Default is extracted from
#'   the input model object
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{quadpts} and \code{npts}
#' @param Theta_nodes an optional matrix of Theta values to be evaluated in the draws for the
#'   sDRF statistics. However, these values are not averaged across, and instead give the bootstrap
#'   confidence intervals at the respective Theta nodes. Useful when following up a large
#'   sDRF or uDRF statistic, for example, to determine where the difference between the test curves are large
#'   (while still accounting for sampling variability). Returns a matrix with observed
#'   variability
#' @param plot logical; plot the 'sDRF' functions for the evaluated sDBF or sDTF values across the
#'    integration grid or, if \code{DIF = TRUE}, the selected items as a faceted plot of individual items?
#'    If plausible parameter sets were obtained/supplied then imputed confidence intervals will be included
# @param type character indicating the type of test scoring function use. Can be 'score' or 'info'
#' @param DIF logical; return a list of item-level imputation properties using the DRF statistics?
#'   These can generally be used as a DIF detection method and as a graphical display for
#'   understanding DIF within each item
#' @param best_fitting logical; use the best fitting parametric distribution (Gaussian by default)
#'  that was used at the time of model estimation? This will result in much fast computations, however
#'  the results are more dependent upon the underlying modelling assumptions. Default is FALSE, which
#'  uses the empirical histogram approach
#' @param p.adjust string to be passed to the \code{\link{p.adjust}} function to adjust p-values.
#'   Adjustments are located in the \code{adj_pvals} element in the returned list. Only applicable when
#'   \code{DIF = TRUE}
#' @param den.type character specifying how the density of the latent traits is computed.
#'   Default is \code{'marginal'} to include the proportional information from both groups,
#'   \code{'focal'} for just the focal group, and \code{'reference'} for the reference group
#' @param groups2test when more than 2 groups are being investigated which two groups
#'   should be used in the effect size comparisons?
#' @param pairwise logical; perform pairwise computations when the applying to multi-group settings
#' @param simplify logical; attempt to simplify the output rather than returning larger lists?
#' @param auto.key plotting argument passed to \code{\link{lattice}}
#' @param par.strip.text plotting argument passed to \code{\link{lattice}}
#' @param par.settings plotting argument passed to \code{\link{lattice}}
#' @param ... additional arguments to be passed to \code{lattice}
#' @param verbose logical; include additional information in the console?
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R. P. (2018). Model-Based Measures for Detecting and Quantifying Response Bias.
#'   \emph{Psychometrika, 83}(3), 696-732. \doi{10.1007/s11336-018-9626-9}
#' @seealso \code{\link{multipleGroup}}, \code{\link{DIF}}
#' @keywords differential response functioning
#' @export DRF
#' @examples
#' \dontrun{
#'
#' set.seed(1234)
#' n <- 30
#' N <- 500
#'
#' # only first 5 items as anchors
#' model <- 'F = 1-30
#'           CONSTRAINB = (1-5, a1), (1-5, d)'
#'
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('Group_1', N), rep('Group_2', N))
#'
#' ## -------------
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                      invariance=c('free_means', 'free_var'))
#' plot(mod)
#' plot(mod, which.items = 6:10) #DBF
#' plot(mod, type = 'itemscore')
#' plot(mod, type = 'itemscore', which.items = 10:15)
#'
#' # empirical histogram approach
#' DRF(mod)
#' DRF(mod, focal_items = 6:10) #DBF
#' DRF(mod, DIF=TRUE)
#' DRF(mod, DIF=TRUE, focal_items = 10:15)
#'
#' # Best-fitting Gaussian distributions
#' DRF(mod, best_fitting=TRUE)
#' DRF(mod, focal_items = 6:10, best_fitting=TRUE) #DBF
#' DRF(mod, DIF=TRUE, best_fitting=TRUE)
#' DRF(mod, DIF=TRUE, focal_items = 10:15, best_fitting=TRUE)
#'
#' DRF(mod, plot = TRUE)
#' DRF(mod, focal_items = 6:10, plot = TRUE) #DBF
#' DRF(mod, DIF=TRUE, plot = TRUE)
#' DRF(mod, DIF=TRUE, focal_items = 10:15, plot = TRUE)
#'
#' if(interactive()) mirtCluster()
#' DRF(mod, draws = 500)
#' DRF(mod, draws = 500, best_fitting=TRUE)
#' DRF(mod, draws = 500, plot=TRUE)
#'
#' # pre-draw parameter set to save computations
#' #  (more useful when using non-parametric bootstrap)
#' param_set <- draw_parameters(mod, draws = 500)
#' DRF(mod, focal_items = 6, param_set=param_set) #DIF test
#' DRF(mod, DIF=TRUE, param_set=param_set) #DIF test
#' DRF(mod, focal_items = 6:10, param_set=param_set) #DBF test
#' DRF(mod, param_set=param_set) #DTF test
#'
#' DRF(mod, focal_items = 6:10, draws=500) #DBF test
#' DRF(mod, focal_items = 10:15, draws=500) #DBF test
#'
#' DIFs <- DRF(mod, draws = 500, DIF=TRUE)
#' print(DIFs)
#' DRF(mod, draws = 500, DIF=TRUE, plot=TRUE)
#'
#' DIFs <- DRF(mod, draws = 500, DIF=TRUE, focal_items = 6:10)
#' print(DIFs)
#' DRF(mod, draws = 500, DIF=TRUE, focal_items = 6:10, plot = TRUE)
#'
#' DRF(mod, DIF=TRUE, focal_items = 6)
#' DRF(mod, draws=500, DIF=TRUE, focal_items = 6)
#'
#' # evaluate specific values for sDRF
#' Theta_nodes <- matrix(seq(-6,6,length.out = 100))
#'
#' sDTF <- DRF(mod, Theta_nodes=Theta_nodes)
#' head(sDTF)
#' sDTF <- DRF(mod, Theta_nodes=Theta_nodes, draws=200)
#' head(sDTF)
#'
#' # sDIF (isolate single item)
#' sDIF <- DRF(mod, Theta_nodes=Theta_nodes, focal_items=6)
#' head(sDIF)
#' sDIF <- DRF(mod, Theta_nodes=Theta_nodes, focal_items = 6, draws=200)
#' head(sDIF)
#'
#' ## -------------
#' ## random slopes and intercepts for 15 items, and latent mean difference
#' ##    (no systematic DTF should exist, but DIF will be present)
#' set.seed(1234)
#' dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 0, .25)),
#'                 d + c(numeric(15), rnorm(n-15, 0, .5)), N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod1 <- multipleGroup(dat, 1, group=group)
#' plot(mod1)
#' DRF(mod1) #does not account for group differences! Need anchors
#'
#' mod2 <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod2)
#'
#' #significant DIF in multiple items....
#' # DIF(mod2, which.par=c('a1', 'd'), items2test=16:30)
#' DRF(mod2)
#' DRF(mod2, draws=500) #non-sig DTF due to item cancellation
#'
#' ## -------------
#' ## systematic differing slopes and intercepts (clear DTF)
#' set.seed(1234)
#' dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)),
#'                 d + c(numeric(15), rnorm(n-15, 1, .5)),
#'                 N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod3 <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod3) #visable DTF happening
#'
#' # DIF(mod3, c('a1', 'd'), items2test=16:30)
#' DRF(mod3) #unsigned bias. Signed bias (group 2 scores higher on average)
#' DRF(mod3, draws=500)
#' DRF(mod3, draws=500, plot=TRUE) #multiple DRF areas along Theta
#'
#' # plot the DIF
#' DRF(mod3, draws=500, DIF=TRUE, plot=TRUE)
#'
#' # evaluate specific values for sDRF
#' Theta_nodes <- matrix(seq(-6,6,length.out = 100))
#' sDTF <- DRF(mod3, Theta_nodes=Theta_nodes, draws=200)
#' head(sDTF)
#'
#' # DIF
#' sDIF <- DRF(mod3, Theta_nodes=Theta_nodes, focal_items = 30, draws=200)
#' car::some(sDIF)
#'
#' ## ----------------------------------------------------------------
#' ### multidimensional DTF
#'
#' set.seed(1234)
#' n <- 50
#' N <- 1000
#'
#' # only first 5 items as anchors within each dimension
#' model <- 'F1 = 1-25
#'           F2 = 26-50
#'           COV = F1*F2
#'           CONSTRAINB = (1-5, a1), (1-5, 26-30, d), (26-30, a2)'
#'
#' a <- matrix(c(rep(1, 25), numeric(50), rep(1, 25)), n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('Group_1', N), rep('Group_2', N))
#' Cov <- matrix(c(1, .5, .5, 1.5), 2)
#' Mean <- c(0, 0.5)
#'
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich', sigma = cov2cor(Cov))
#' dat2 <- simdata(a, d, N, itemtype = 'dich', sigma = Cov, mu = Mean)
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                      invariance=c('free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#' plot(mod, degrees = c(45,45))
#' DRF(mod)
#'
#' # some intercepts slightly higher in Group 2
#' d2 <- d
#' d2[c(10:15, 31:35)] <- d2[c(10:15, 31:35)] + 1
#' dat1 <- simdata(a, d, N, itemtype = 'dich', sigma = cov2cor(Cov))
#' dat2 <- simdata(a, d2, N, itemtype = 'dich', sigma = Cov, mu = Mean)
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                      invariance=c('free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#' plot(mod, degrees = c(45,45))
#'
#' DRF(mod)
#' DRF(mod, draws = 500)
#'
# ####
# # bifactor model
#
# #simulate data
# a <- matrix(c(rlnorm(30, .2, .2), rlnorm(15, .2, .2),
#             numeric(30), rlnorm(15, .2, .2)), 30)
# d <- matrix(rnorm(30))
#
# sigma <- sigma2 <- diag(3)
# sigma2[1,1] <- 2
# d2 <- d + matrix(c(numeric(5), rnorm(25, 0, .1)))
# dat1 <- simdata(a, d, 1000, itemtype='dich', sigma=sigma)
# dat2 <- simdata(a, d2, 1000, itemtype='dich', sigma=sigma2, mu=c(0.5, 0, 0))
# dat <- rbind(dat1, dat2)
# group <- rep(c('group1', 'group2'), each = 1000)
#
# specific <- c(rep(1,15),rep(2,15))
# nms <- colnames(dat)
# simmod <- bfactor(dat, specific, group=group, SE=TRUE,
#                   invariance = c('free_means', 'free_var', nms[1:5]))
# coef(simmod, simplify=TRUE)
#
# DRF(simmod)
# DRF(simmod, draws = 500)
#
#' }
DRF <- function(mod, draws = NULL, focal_items = 1L:extract.mirt(mod, 'nitems'),
                param_set = NULL, den.type = 'marginal', best_fitting=FALSE,
                CI = .95, npts = 1000,
                quadpts = NULL, theta_lim=c(-6,6), Theta_nodes = NULL,
                plot = FALSE, DIF = FALSE, groups2test = 'all',
                pairwise = FALSE, simplify = TRUE, p.adjust = 'none',
                par.strip.text = list(cex = 0.7),
                par.settings = list(strip.background = list(col = '#9ECAE1'),
                                 strip.border = list(col = "black")),
                auto.key = list(space = 'right', points=FALSE, lines=TRUE),
                verbose = TRUE, ...){

    compute_ps <- function(x, xs, X2=FALSE){
        if(X2){
            if(ncol(xs) > 2L){
                se <- apply(xs, 2L, sd)
                se <- ifelse(se == 0, 1, se)
                tmp <- (as.vector(x) / se)^2
                ts2 <- tmp[1L:(length(tmp) / 2L)] + tmp[(length(tmp)/2L + 1L):length(tmp)]
                df <- 2L - abs(rowSums(x == 0))
            } else {
                se <- apply(xs, 2L, sd)
                se <- ifelse(se == 0, 1, se)
                ts2 <- sum((x / se)^2)
                df <- 2L - abs(sum(x == 0))
            }
            p <- unname(sapply(1L:length(df), function(i) pchisq(ts2[i], df[i], lower.tail = FALSE)))
            ret <- cbind(X2 = unname(ts2), df=df, p = p)
            ret[p==1, ] <- NA
        } else {
            if(length(x) > 1L){
                ts <- x / apply(xs, 2L, sd)
            } else {
                ts <- x / sd(xs)
            }
            p <- unname(pnorm(abs(ts), lower.tail = FALSE) * 2)
            ret <- cbind(X2 = unname(ts^2), df=ifelse(p==1, NA, 1L), p = p)
            ret[is.nan(p), ] <- NA
        }
        ret
    }

    fn <- function(x, omod, Theta, max_score, Theta_nodes = NULL, best_fitting,
                   plot, DIF, focal_items, details, whc_grp,
                   signs=NULL, rs=NULL, den.type){
        mod <- omod
        if(!is.null(Theta_nodes)){
            T1 <- expected.test(mod, Theta_nodes, group=whc_grp[1L], mins=FALSE, individual=DIF,
                                which.items=focal_items)
            T2 <- expected.test(mod, Theta_nodes, group=whc_grp[2L], mins=FALSE, individual=DIF,
                                which.items=focal_items)
            ret <- T2 - T1
            if(!DIF) ret <- c("sDRF." = ret)
            return(ret)
        }
        calc_DRFs(mod=mod, Theta=Theta, plot=plot, max_score=max_score, DIF=DIF, den.type=den.type,
                  focal_items=focal_items, details=details, signs=signs, rs=rs,
                  best_fitting=best_fitting, whc_grp=whc_grp)
    }
    fn2 <- function(ind, pars, MGmod, param_set, rslist, whc_grp, ...){
        pars <- reloadPars(longpars=param_set[ind,],
                           pars=pars, ngroups=length(pars), J=length(pars[[1L]])-1L)
        rs <- rslist[[ind]]
        MGmod@ParObjects$pars[[whc_grp[1L]]]@ParObjects$pars <- pars[[whc_grp[1L]]]
        MGmod@ParObjects$pars[[whc_grp[2L]]]@ParObjects$pars <- pars[[whc_grp[2L]]]
        fn(NA, omod=MGmod, rs=rs, whc_grp=whc_grp, ...)
    }

    if(missing(mod)) missingMsg('mod')
    if(pairwise){
        if(length(groups2test) == 1L && groups2test == 'all')
            groups2test <- extract.mirt(mod, 'groupNames')
        ngroups <- length(groups2test)
        compare <- vector('list', ngroups*(ngroups-1L)/2L)
        ind <- 1L
        for(i in 1L:ngroups){
            for(j in 1L:ngroups){
                if(i < j){
                    compare[[ind]] <- DRF(mod=mod, draws=draws, focal_items=focal_items,
                                          groups2test = groups2test[c(i,j)],
                                          param_set=param_set, den.type=den.type, best_fitting=best_fitting,
                                          CI=CI, npts=npts, quadpts=quadpts, theta_lim=theta_lim,
                                          Theta_nodes=Theta_nodes, plot = FALSE, DIF=DIF, pairwise=FALSE,
                                          p.adjust=p.adjust, par.strip.text=par.strip.text,
                                          par.settings=par.settings, auto.key=auto.key,
                                          verbose=verbose, ...)
                    ind <- ind + 1L
                }
            }
        }
        if(simplify && is(compare[[1L]], 'mirt_df')){
            compare <- do.call(rbind, compare)
            rownames(compare) <- NULL
        }
        if(DIF && !is.null(draws) && !is.null(compare[[1L]]$sDIF)){
            nms <- sapply(compare, function(x) x$sDIF$groups[1L])
            names(compare) <- nms
        }
        return(compare)
    }
    stopifnot(length(p.adjust) == 1L)
    groupNames <- extract.mirt(mod, 'groupNames')
    if(length(groups2test) == 1L && groups2test == 'all')
        groups2test <- groupNames
    stopifnot(den.type %in% c('marginal', 'focal', 'reference'))
    stopifnot(is.logical(plot))
    if(DIF && !is.null(Theta_nodes))
        stop('DIF must be FALSE when using Theta_nodes', call.=FALSE)
    type <- 'score'
    if(!is(mod, 'MultipleGroupClass'))
        stop('mod input was not estimated by multipleGroup()', call.=FALSE)
    if(length(groups2test) > 2L)
        stop('DRF only supports two group models at a time', call.=FALSE)
    if(!any(sapply(mod@ParObjects$pars, function(x) x@ParObjects$pars[[length(x@ParObjects$pars)]]@est)))
        message('No hyper-parameters were estimated in the DIF model. For effective
                \tDRF testing freeing the focal group hyper-parameters is recommended.')
    if(!is.null(Theta_nodes)){
        if(!is.matrix(Theta_nodes))
            stop('Theta_nodes must be a matrix', call.=FALSE)
        if(ncol(Theta_nodes) != mod@Model$nfact)
            stop('Theta_nodes input does not have the correct number of factors', call.=FALSE)
        colnames(Theta_nodes) <- if(ncol(Theta_nodes) > 1L)
            paste0('Theta.', 1L:ncol(Theta_nodes)) else 'Theta'
    }
    if(extract.mirt(mod, 'nfact') != 1L && plot)
        stop('plot arguments only supported for unidimensional models')
    if(length(type) > 1L && (plot || !is.null(Theta_nodes)))
        stop('Multiple type arguments cannot be combined with plot or Theta_nodes arguments')
    m2v <- mod2values(mod)
    is_logit <- m2v$name %in% c('g', 'u')
    longpars <- do.call(c, lapply(1L:length(groupNames), function(ind)
        do.call(c, lapply(mod@ParObjects$pars[[ind]]@ParObjects$pars, function(x) x@par))))
    whc_grp <- which(groupNames %in% groups2test)
    if(!is.null(param_set)){
        draws <- nrow(param_set)
        param_set[,is_logit] <- logit(param_set[,is_logit])
        impute <- TRUE
    }
    if(is.null(draws)){
        draws <- 1L
        impute <- FALSE
    } else if(is.null(param_set)){
        if(length(mod@vcov) == 1L)
            stop('Stop an information matrix must be computed', call.=FALSE)
        if(!mod@OptimInfo$secondordertest)
            stop('ACOV matrix is not positive definite')
        impute <- TRUE
        covB <- mod@vcov
        names <- colnames(covB)
        pars <- lapply(1L:length(groupNames), function(ind)
            mod@ParObjects$pars[[ind]]@ParObjects$pars)
        param_set <- draw_parameters(mod, draws=draws, ...)
        param_set[,is_logit] <- logit(param_set[,is_logit])
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=length(groupNames),
                           J=length(pars[[1L]])-1L)
    }
    shortpars <- mod@Internals$shortpars
    names <- names(shortpars)
    if(is.null(quadpts)) quadpts <- mod@Options$quadpts
    if(is.null(Theta_nodes)){
        theta <- matrix(seq(theta_lim[1L], theta_lim[2L], length.out=quadpts))
        Theta <- thetaComb(theta, mod@Model$nfact)
    } else Theta <- Theta_nodes
    max_score <- sum(mod@Data$mins + mod@Data$K - 1L)

    large <- multipleGroup(extract.mirt(mod, 'data'), 1, group=extract.mirt(mod, 'group'),
                           large='return')
    details <- list(data = extract.mirt(mod, 'data'),
                    model = extract.mirt(mod, 'model'),
                    group = extract.mirt(mod, 'group'),
                    constrain = extract.mirt(mod, 'constrain'),
                    itemtype = extract.mirt(mod, 'itemtype'),
                    customItems = extract.mirt(mod, 'customItems'),
                    customGroup = extract.mirt(mod, 'customGroup'),
                    technical = list(storeEtable=TRUE, theta_lim=theta_lim, omp=FALSE),
                    quadpts=quadpts, large=large, TOL = NaN)
    if(plot) Theta_nodes <- matrix(seq(theta_lim[1L], theta_lim[2L], length.out=1000))
    oCM <- lapply(1L, fn, omod=mod, Theta_nodes=Theta_nodes, best_fitting=best_fitting,
                  max_score=max_score, Theta=Theta, plot=plot, den.type=den.type, whc_grp=whc_grp,
                  DIF=DIF, focal_items=focal_items, details=details)[[1L]]
    signs <- attr(oCM, 'signs')
    if(plot && !impute) return(plot.DRF(Theta_nodes, oCM, DIF=DIF,
                                 itemnames = extract.mirt(mod, 'itemnames')[focal_items], ...))
    if(!is.null(Theta_nodes) && !impute)
        return(as.mirt_df(data.frame(Theta=Theta_nodes, sDRF=oCM)))
    groups_tested <- paste0(groups2test,collapse=',')
    if(impute){
        pars <- lapply(1L:length(groupNames), function(ind)
            mod@ParObjects$pars[[ind]]@ParObjects$pars)
        .mirtClusterEnv$param_set <- param_set
        try(with(details, multipleGroup(data=data, model=model, group=group, itemtype=itemtype, large=large,
                                    quadpts=quadpts, TOL=TOL,
                                    customItems = extract.mirt(mod, 'customItems'),
                                    customGroup = extract.mirt(mod, 'customGroup'),
                                    pars=mod2values(mod), technical=technical)), TRUE)
        rslist <- .mirtClusterEnv$rslist
        on.exit({.mirtClusterEnv$rslist <- .mirtClusterEnv$param_set <- NULL
            reloadPars(longpars=longpars, pars=pars, ngroups=length(groupNames),
                       J=length(pars[[1L]])-1L)
        }, add=TRUE)
        list_scores <- myLapply(1L:nrow(param_set), fn2, whc_grp=whc_grp,
                                pars=pars, MGmod=mod, param_set=param_set, best_fitting=best_fitting,
                                max_score=max_score, Theta=Theta, rslist=rslist,
                                Theta_nodes=Theta_nodes, plot=plot, details=details, progress=verbose,
                                DIF=DIF, focal_items=focal_items, signs=signs, den.type=den.type)
        scores <- do.call(rbind, list_scores)
        pars <- lapply(1L:length(groupNames), function(ind)
            mod@ParObjects$pars[[ind]]@ParObjects$pars)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=length(groupNames),
                           J=length(pars[[1L]])-1L)
        if(plot) return(plot.DRF(Theta_nodes, oCM, CIs=scores, DIF=DIF, CI=CI,
                                 itemnames = extract.mirt(mod, 'itemnames')[focal_items], ...))
        CIs <- apply(scores, 2, bs_range, CI=CI)
        if(is.null(Theta_nodes))
            CIs <- CIs[,1L:(ncol(CIs) * 3/5)]
        rownames(CIs) <- c(paste0('CI_', round((1-CI)/2, 3L)*100),
                           paste0('CI_', round(CI + (1-CI)/2, 3L)*100))
        if(!is.null(Theta_nodes))
            return(as.mirt_df(data.frame(Theta=Theta_nodes, sDRF=oCM, t(CIs))))
        if(DIF){
            oCM <- matrix(oCM, length(focal_items))
            t1 <- compute_ps(oCM[,1L], scores[,1L:length(focal_items), drop=FALSE])
            t2 <- compute_ps(oCM[,3L:4L], scores[,1L:(length(focal_items)*2L) +
                                                     length(focal_items)*2L, drop=FALSE],
                             X2=TRUE)
            ret <- list(sDIF = as.mirt_df(data.frame(
                                groups=groups_tested,
                                item=extract.mirt(mod, 'itemnames')[focal_items],
                                sDIF = oCM[,1L], t(CIs[,1L:length(focal_items)]), t1,
                                row.names=NULL)),
                uDIF = as.mirt_df(data.frame(
                            groups=groups_tested,
                            item=extract.mirt(mod, 'itemnames')[focal_items],
                            uDIF = oCM[,2L],
                            t(CIs[,1L:length(focal_items) + length(focal_items)]),
                            t2, row.names=NULL)),
                dDIF = as.mirt_df(data.frame(
                            groups=groups_tested,
                            item=extract.mirt(mod, 'itemnames')[focal_items],
                            dDIF=oCM[,3L], t(CIs[,1L:length(focal_items) + length(focal_items)*2]),
                            row.names=NULL)))
            if(p.adjust != 'none'){
                ret$sDIF$adj_pvals <- p.adjust(ret$sDIF$p, method=p.adjust)
                ret$uDIF$adj_pvals <- p.adjust(ret$uDIF$p, method=p.adjust)
            }
        } else {
            t1 <- compute_ps(oCM[1L], scores[,1L])
            t2 <- compute_ps(oCM[3L:4L], scores[,3L:4L], X2=TRUE)
            tests <- rbind(t1, t2, NA)
            ret <- data.frame(groups=groups_tested,
                              n_focal_items=length(focal_items),
                              stat = oCM[1L:3L], t(CIs), tests, check.names = FALSE)
            ret <- as.mirt_df(ret)
        }
    } else {
        # no imputations
        if(DIF){
            ret <- data.frame(groups=groups_tested,
                              item=extract.mirt(mod, 'itemnames')[focal_items],
                              matrix(oCM, length(oCM)/5L), row.names=NULL)
            ret <- ret[,-c(6L:7L)]
            colnames(ret) <- c('groups', 'item', 'sDIF', 'uDIF', 'dDIF')
            ret <- as.mirt_df(ret)
        } else {
            ret <- data.frame(groups=groups_tested,
                              n_focal_items=length(focal_items),
                              sDRF=oCM[1L], uDRF=oCM[2L], dDRF=oCM[3L],
                              row.names=NULL)
            ret <- as.mirt_df(ret)
        }
    }
    ret
}

calc_DRFs <- function(mod, Theta, DIF, plot, max_score, focal_items, details, den.type,
                      best_fitting, whc_grp, rs=NULL, signs=NULL){
    if(DIF){
        T1 <- expected.test(mod, Theta, group=whc_grp[1L], mins=FALSE, individual = TRUE,
                            which.items=focal_items)
        T2 <- expected.test(mod, Theta, group=whc_grp[2L], mins=FALSE, individual = TRUE,
                            which.items=focal_items)
    } else {
        T1 <- matrix(expected.test(mod, Theta, group=whc_grp[1L],
                            mins=FALSE, which.items=focal_items))
        T2 <- matrix(expected.test(mod, Theta, group=whc_grp[2L],
                            mins=FALSE, which.items=focal_items))
    }
    if(plot) return(c(T1, T2))
    if(best_fitting){
        pars <- extract.mirt(mod, 'pars')
        nitems <- extract.mirt(mod, 'nitems')
        den1 <- pars[[whc_grp[1L]]]@ParObjects$pars[[nitems + 1L]]@den(
            pars[[whc_grp[1L]]]@ParObjects$pars[[nitems + 1L]], Theta)
        den1 <- den1 / sum(den1)
        den2 <- pars[[whc_grp[2L]]]@ParObjects$pars[[nitems + 1L]]@den(
            pars[[whc_grp[2L]]]@ParObjects$pars[[nitems + 1L]], Theta)
        den2 <- den2 / sum(den2)
        p <- if(den.type == 'marginal'){
            Ns <- table(extract.mirt(mod, 'group'))
            (Ns[whc_grp[1L]] * den1 + Ns[whc_grp[2L]] * den2) / sum(Ns[whc_grp])
        } else if(den.type == 'focal') den2
          else den1
    } else {
        if(is.null(rs)){
            mod2 <- with(details, multipleGroup(data=data, model=model, group=group,
                                                itemtype=itemtype, large=large,
                                                constrain=constrain, quadpts=quadpts, TOL=TOL,
                                                customItems = customItems,
                                                customGroup = customGroup,
                                                pars=mod2values(mod), technical=technical))
            r1 <- rowSums(mod2@Internals$Etable[[whc_grp[1L]]]$r1)
            r2 <- rowSums(mod2@Internals$Etable[[whc_grp[2L]]]$r1)
        } else {
            r1 <- rs[,whc_grp[1L]]
            r2 <- rs[,whc_grp[2L]]
        }
        p <- if(den.type == 'marginal')
            (r1 + r2) / sum(r1 + r2)
        else if(den.type == 'focal')
            r2 / sum(r2)
        else r1/ sum(r1)
    }
    D <- T2 - T1
    uDRF <- colSums(abs(D) * p)
    sDRF <- colSums(D * p)
    dDRF <- sqrt(colSums(D^2 * p))
    attach_signs <- FALSE
    ret <- if(is.null(signs)){
        signs <- D < 0
        attach_signs <- TRUE
    }
    uDRF_L <- colSums(D * p * signs)
    uDRF_U <- colSums(D * p * !signs)
    ret <- c(sDRF=sDRF, uDRF=uDRF, dDRF=dDRF,
             uDRF_L=uDRF_L, uDRF_U=uDRF_U)
    if(attach_signs) attr(ret, 'signs') <- signs
    ret
}

#' Draw plausible parameter instantiations from a given model
#'
#' Draws plausible parameters from a model using parametric sampling (if the information matrix
#' was computed) or via bootstrap sampling. Primarily for use with the \code{\link{DRF}} function.
#'
#' @param mod estimated single or multiple-group model
#' @param draws number of draws to obtain
#' @param method type of plausible values to obtain. Can be 'parametric', for the parametric sampling
#'   scheme which uses the estimated information matrix, or 'boostrap' to obtain values from the \code{\link{boot}}
#'   function. Default is 'parametric'
#' @param redraws number of redraws to perform when the given parameteric sample does not satisfy the
#'   upper and lower parameter bounds. If a valid set cannot be found within this number of draws then
#'   an error will be thrown
#' @param verbose logical; include additional information in the console?
#' @param ... additional arguments to be passed
#' @return returns a draws x p matrix of plausible parameters, where each row correspeonds to a single
#'   set
#'
#' @export
#' @examples
#'
#' \dontrun{
#' set.seed(1234)
#' n <- 40
#' N <- 500
#'
#' # only first 5 items as anchors
#' model <- 'F = 1-40
#'           CONSTRAINB = (1-5, a1), (1-5, d)'
#'
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('Group_1', N), rep('Group_2', N))
#'
#' ## -------------
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                      invariance=c('free_means', 'free_var'))
#'
#' param_set <- draw_parameters(mod, 100)
#' head(param_set)
#' }
#'
draw_parameters <- function(mod, draws, method = c('parametric', 'boostrap'),
                            redraws = 20, verbose = FALSE, ...){
    fn_param <- function(ind, shortpars, longpars, lbound, ubound, est,
                         pre.ev, constrain, imputenums, MGmod, redraws, pars){
        count <- 0L
        while(TRUE){
            count <- count + 1L
            if(count == redraws)
                stop('Invalid parameter set drawn too often', call.=FALSE)
            shift <- mirt_rmvnorm(1L, mean=shortpars, pre.ev=pre.ev)
            longpars[imputenums] <- shift[1L,]
            for(i in seq_len(length(constrain)))
                longpars[constrain[[i]][-1L]] <- longpars[constrain[[i]][1L]]
            if(any((longpars < lbound | longpars > ubound) & est)) next
            if(any(MGmod@Model$itemtype %in% c('graded', 'grsm'))){
                pars <- reloadPars(longpars=longpars, pars=pars, ngroups=2L, J=length(pars[[1L]])-1L)
                pick <- c(MGmod@Model$itemtype %in% c('graded', 'grsm'), FALSE)
                if(!all(sapply(pars[[1L]][pick], CheckIntercepts) &
                        sapply(pars[[2L]][pick], CheckIntercepts))) next
            }
            return(longpars)
        }
    }

    method <- match.arg(method)
    shortpars <- mod@Internals$shortpars
    m2v <- mod2values(mod)
    longpars <- m2v$value
    logits <- m2v$name %in% c('g', 'u')
    lbound <- m2v$lbound
    ubound <- m2v$ubound
    longpars[logits] <- logit(longpars[logits])
    lbound[logits] <- logit(lbound[logits])
    ubound[logits] <- logit(ubound[logits])
    est <- m2v$est
    constrain <- mod@Model$constrain
    ngroups <- extract.mirt(mod, 'ngroups')
    pars <- if(ngroups > 1L){
        lapply(1:ngroups, function(i) mod@ParObjects$pars[[i]]@ParObjects$pars)
    } else {
        list(mod@ParObjects$pars)
    }

    if(method == 'parametric'){
        if(!mod@OptimInfo$secondordertest)
            stop('ACOV matrix is not positive definite')
        on.exit(reloadPars(longpars=longpars, pars=pars,
                           ngroups=ngroups, J=extract.mirt(mod, 'nitems')), add=TRUE)
        covB <- vcov(mod)
        names <- colnames(covB)
        imputenums <- sapply(strsplit(names, '\\.'), function(x) as.integer(x[2L]))
        pre.ev <- eigen(covB)
        ret <- myLapply(1L:draws, fn_param, progress=verbose,
                        shortpars=shortpars, longpars=longpars, lbound=lbound,
                        ubound=ubound, pre.ev=pre.ev, constrain=constrain, est=est,
                        imputenums=imputenums, MGmod=mod, redraws=redraws, pars=pars)
        if(verbose) cat("\n")
        ret <- do.call(rbind, ret)
        if(any(logits))
            ret[,logits] <- antilogit(ret[,logits])
        return(ret)
    } else stop('bootstrap not supported yet') #TODO

}

plot.DRF <- function(Theta, DV, itemnames, CIs = NULL, DIF = FALSE, CI,
                     main = 'Signed DIF', par.strip.text = list(cex = 0.7),
                     par.settings = list(strip.background = list(col = '#9ECAE1'),
                                         strip.border = list(col = "black")),
                     auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...){

    panel.bands <- function(x, y, upper, lower, fill, col,
                            subscripts, ..., font, fontface){
                                upper <- upper[subscripts]
                                lower <- lower[subscripts]
                                panel.polygon(c(x, rev(x)), c(upper, rev(lower)),
                                              col = fill, border = FALSE, ...)
    }

    nquad <- nrow(Theta)
    ID <- 1L:nquad
    nfact <- ncol(Theta)
    stopifnot(nfact < 2L)
    if(nfact == 1) colnames(Theta) <- 'Theta'
    else colnames(Theta) <- c('Theta.1', 'Theta.2')
    if(is.null(CIs)){
        if(DIF){
            item <- factor(rep(itemnames, each=length(ID)), levels = itemnames)
            dat <- data.frame(ID, Theta, as.vector(DV), item=item)
            return(lattice::xyplot(DV ~ Theta | item, dat, main = main, col = 'black',
                                   panel = function(x, y, ...){
                                       panel.xyplot(x, y, type='l', lty=1,...)
                                       panel.abline(h = 0, col = 'red')
                                   }, ylab = expression(sDIF[theta]), xlab = expression(theta),
                                   par.settings=par.settings,
                                   auto.key=auto.key, par.strip.text=par.strip.text, ...))
        } else {
            dat <- data.frame(ID, Theta, DV)
            main <-  'Signed DRF'
            return(lattice::xyplot(DV ~ Theta, dat, main = main, col = 'black',
                                   panel = function(x, y, ...){
                                       panel.xyplot(x, y, type='l', lty=1,...)
                                       panel.abline(h = 0, col = 'red')
                                   }, ylab = expression(sDRF[theta]), xlab = expression(theta),
                                   par.settings=par.settings,
                                   auto.key=auto.key, par.strip.text=par.strip.text, ...))
        }
    } else {
        if(DIF){
            item <- factor(rep(itemnames, each=length(ID)), levels = itemnames)
            draws <- nrow(CIs) / nquad
            tmp <- lapply(1L:nquad, function(x){
                pick <- seq(x, draws*nquad, by = nquad)
                t(apply(CIs[pick, , drop=FALSE], 2L, bs_range, CI=CI))
            })
            tmp2 <- do.call(rbind, tmp)
            lower <- upper <- c()
            for(i in seq_len(ncol(CIs))){
                pick <- seq(i, nrow(tmp2), by=ncol(CIs))
                lower <- c(lower, tmp2[pick,1L])
                upper <- c(upper, tmp2[pick,2L])
            }
            dat <- data.frame(ID, Theta, as.vector(DV), item=item, lower=lower, upper=upper)
            return(lattice::xyplot(DV ~ Theta | item, dat, main = main, lower=dat$lower, upper=dat$upper,
                                   fill = 'darkgrey', alpha = .2,
                                   panel = function(x, y, lower, upper, fill, alpha, subscripts, ...){
                                       panel.xyplot(c(min(x), max(x)), c(0,0), col = 'red', type = 'l')
                                       panel.xyplot(x, upper[subscripts], col='black',
                                                    lty=2, type='l', alpha = .4, ...)
                                       panel.xyplot(rev(x), rev(lower[subscripts]), col='black',
                                                    lty=2, type='l', alpha = .4, ...)
                                       panel.polygon(c(x, rev(x)), c(upper[subscripts], rev(lower[subscripts])),
                                                     col = fill, border = FALSE, alpha=alpha, ...)
                                       panel.xyplot(x, y, type='l', lty=1, col = 'black', ...)
                                   }, ylab = expression(sDIF[theta]), xlab = expression(theta),
                                   ylim = c(min(dat$lower), max(dat$upper)),
                                   par.settings=par.settings,
                                   auto.key=auto.key, par.strip.text=par.strip.text, ...))

        } else {

            tmp <- apply(t(CIs), 1L, bs_range, CI=CI)
            dat <- data.frame(ID, Theta, as.vector(DV), lower=tmp[1L,], upper=tmp[2L,])
            main <-  'Signed DRF'
            return(lattice::xyplot(DV ~ Theta, dat, main = main, lower=dat$lower, upper=dat$upper,
                                   fill = 'darkgrey', alpha = .2,
                                   panel = function(x, y, lower, upper, fill, alpha, subscripts, ...){
                                       panel.xyplot(c(min(x), max(x)), c(0,0), col = 'red', type = 'l')
                                       panel.xyplot(x, upper[subscripts], col='black',
                                                    lty=2, type='l', alpha = .4, ...)
                                       panel.xyplot(rev(x), rev(lower[subscripts]), col='black',
                                                    lty=2, type='l', alpha = .4, ...)
                                       panel.polygon(c(x, rev(x)), c(upper[subscripts], rev(lower[subscripts])),
                                                     col = fill, border = FALSE, alpha=alpha, ...)
                                       panel.xyplot(x, y, type='l', lty=1, col = 'black', ...)
                                   }, ylab = expression(sDIF[theta]), xlab = expression(theta),
                                   ylim = c(min(dat$lower), max(dat$upper)),
                                   par.settings=par.settings,
                                   auto.key=auto.key, par.strip.text=par.strip.text, ...))

        }
    }

}
