#' Differential test functioning statistics
#'
#' Function performs various omnibus differential test functioning procedures on an object
#' estimated with \code{multipleGroup()}. If the latent means/covariances are suspected to differ
#' then the input object should contain a set of 'anchor' items to ensure that only differential
#' test features are being detected rather than group differences. Returns signed (average area
#' above and below) and unsigned (total area) statistics, with descriptives such as the percent
#' average bias between group total scores for each statistic. If a grid of Theta values is passed,
#' these can be evaluated as well to determine specific DTF location effects.  For best results,
#' the baseline model should contain a set of 'anchor' items and have freely estimated
#' hyper-parameters in the focal groups. See \code{\link{DIF}} for details.
#'
#' @aliases DTF
#' @param mod a multipleGroup object which estimated only 2 groups
#' @param MI a number indicating how many draws to take to form a suitable multiple imputation
#'   for the expected test scores (100 or more). Requires an estimated parameter
#'   information matrix. Returns a list containing the bootstrap distribution and null hypothesis
#'   test for the sDTF statistic
#' @param CI range of confidence interval when using MI input
#' @param npts number of points to use in the integration. Default is 1000
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param Theta_nodes an optional matrix of Theta values to be evaluated in the MI draws for the
#'   sDTF statistic. However, these values are not averaged across, and instead give the bootstrap
#'   confidence intervals at the respective Theta nodes. Useful when following up a large
#'   uDTF/sDTF statistic to determine where the difference between the test curves are large
#'   (while still accounting for sampling variability). Returns a matrix with observed
#'   variability
#' @param plot logical; plot the test score functions with imputed confidence envelopes?
#' @param auto.key logical; automatically generate key in lattice plot?
#' @param ... additional arguments to be passed to lattice
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R. P., Counsell, A., and Flora, D. B. (in press). It might not
#'   make a big DIF: Improved Differential Test Functioning statistics that account for
#'   sampling variability. \emph{Educational and Psychological Measurement}.
#' @seealso \code{\link{multipleGroup}}, \code{\link{DIF}}
#' @keywords DTF
#' @export DTF
#' @examples
#' \dontrun{
#' set.seed(1234)
#' n <- 30
#' N <- 500
#'
#' # only first 5 items as anchors
#' model <- mirt.model('F = 1-30
#'                     CONSTRAINB = (1-5, a1), (1-5, d)')
#'
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('G1', N), rep('G2', N))
#'
#' ## -------------
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = 'dich')
#' dat2 <- simdata(a, d, N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type='crossprod',
#'                      invariance=c('free_means', 'free_var'))
#' plot(mod)
#'
#' DTF(mod)
#' mirtCluster()
#' DTF(mod, MI = 1000) #95% C.I. for sDTF containing 0. uDTF is very small
#'
#' ## -------------
#' ## random slopes and intercepts for 15 items, and latent mean difference
#' ##    (no systematic DTF should exist, but DIF will be present)
#' dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), sign(rnorm(n-15))*runif(n-15, .25, .5)),
#'                 d + c(numeric(15), sign(rnorm(n-15))*runif(n-15, .5, 1)), N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod1 <- multipleGroup(dat, 1, group=group)
#' plot(mod1) #does not account for group differences! Need anchors
#'
#' mod2 <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type = 'crossprod',
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod2)
#'
#' #significant DIF in multiple items....
#' DIF(mod2, which.par=c('a1', 'd'), items2test=16:30)
#' DTF(mod2)
#' DTF(mod2, MI=1000)
#'
#' ## -------------
#' ## systematic differing slopes and intercepts (clear DTF)
#' dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, .5)),
#'                 N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod3 <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type='crossprod',
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod3) #visable DTF happening
#'
#' DIF(mod3, c('a1', 'd'), items2test=16:30)
#' DTF(mod3) #unsigned bias. Signed bias indicates group 2 scores generally higher on average
#' DTF(mod3, MI=1000)
#' DTF(mod3, MI=1000, plot=TRUE)
#'
#' }
DTF <- function(mod, MI = NULL, CI = .95, npts = 1000, theta_lim=c(-6,6), Theta_nodes = NULL,
                plot = FALSE, auto.key = TRUE, ...){

    fn <- function(x, omod, impute, covBs, imputenums, Theta, max_score, Theta_nodes = NULL,
                   plot){
        mod <- omod
        if(impute){
            for(g in 1L:2L){
                while(TRUE){
                    tmp <- try(imputePars(pars=omod@pars[[g]]@pars, covB=covBs[[g]],
                                          imputenums=imputenums[[g]], constrain=omod@constrain),
                               silent=TRUE)
                    if(!is(tmp, 'try-error')) break
                }
                mod@pars[[g]]@pars <- tmp
            }
        }
        if(!is.null(Theta_nodes)){
            T1 <- expected.test(mod, Theta_nodes, group=1L, mins=FALSE)
            T2 <- expected.test(mod, Theta_nodes, group=2L, mins=FALSE)
            D <- T1 - T2
            ret <- c("sDTF." = D)
            return(ret)
        }
        T1 <- expected.test(mod, Theta, group=1L)
        T2 <- expected.test(mod, Theta, group=2L)
        if(plot) return(c(T1, T2))
        D <- T1 - T2
        uDTF <- mean(abs(D))
        uDTF_percent <- uDTF/max_score * 100
        sDTF <- mean(D)
        sDTF_percent <- sDTF/max_score * 100
        ret <- c("signed.DTF" = sDTF, "signed.DTF(%)" = sDTF_percent,
                    "unsigned.DTF" = uDTF, "unsigned.DTF(%)" = uDTF_percent)
        ret
    }
    if(missing(mod)) missingMsg('mod')
    if(class(mod) != 'MultipleGroupClass')
        stop('mod input was not estimated by multipleGroup()')
    if(length(mod@pars) != 2L)
        stop('DTF only supports two group models at a time')
    if(!any(sapply(mod@pars, function(x) x@pars[[length(x@pars)]]@est)))
        message('No hyper-parameters were estimated in the DIF model. For effective
                \tDTF testing, freeing the focal group hyper-parameters is recommend.')
    if(!is.null(Theta_nodes)){
        if(!is.matrix(Theta_nodes))
            stop('Theta_nodes must be a matrix')
        if(ncol(Theta_nodes) != mod@nfact)
            stop('Theta_nodes input does not have the correct number of factors')
        colnames(Theta_nodes) <- if(ncol(Theta_nodes) > 1)
            paste0('Theta.', 1:ncol(Theta_nodes)) else 'Theta'
    }
    if(plot){
        if(is.null(MI))
            stop('Must specificy number of imputations to generate plot')
        Theta_nodes <- NULL
    }

    J <- length(mod@K)
    if(is.null(MI)){
        MI <- 1L
        impute <- FALSE
    } else {
        if(length(mod@information) == 1L)
            stop('Stop an information matrix must be computed')
        info <- mod@information
        is_na <- is.na(diag(info))
        info <- info[!is_na, !is_na]
        if(is(try(chol(info), silent=TRUE), 'try-error')){
            stop('Proper information matrix must be precomputed in model')
        } else {
            impute <- TRUE
            list_scores <- vector('list', MI)
            mod <- assignInformationMG(mod)
            covBs <- try(lapply(mod@pars, function(x){
                info <- x@information
                is_na <- is.na(diag(info))
                return(solve(info[!is_na, !is_na]))
                }), silent=TRUE)
            if(is(covBs, 'try-error'))
                stop('Could not compute inverse of information matrix')
            imputenums <- vector('list', 2L)
            for(g in 1L:2L){
                names <- colnames(covBs[[g]])
                tmp <- lapply(names, function(x, split){
                    as.numeric(strsplit(x, split=split)[[1L]][-1L])
                }, split='\\.')
                imputenums[[g]] <- do.call(c, tmp)
            }
        }
    }

    theta <- matrix(seq(theta_lim[1L], theta_lim[2L], length.out=npts))
    if(mod@nfact != 1L)
        stop('DTF only supports unidimensional tests for now.')
    Theta <- thetaComb(theta, mod@nfact)
    max_score <- sum(mod@Data$mins + mod@Data$K - 1L)
    list_scores <- myLapply(1L, fn, omod=mod, impute=FALSE, covBs=NULL, Theta_nodes=Theta_nodes,
                            imputenums=NULL, max_score=max_score, Theta=Theta, plot=plot)
    if(impute){

        bs_range <- function(x, CI){
            ss <- sort(x)
            N <- length(ss)
            ret <- c(upper = ss[ceiling(N * (1 - (1-CI)/2))],
                     middle = median(x),
                     lower = ss[floor(N * (1-CI)/2)])
            ret
        }

        oCM <- list_scores[[1L]]
        list_scores <- myLapply(1L:MI, fn, omod=mod, impute=TRUE, covBs=covBs, max_score=max_score,
                                imputenums=imputenums, Theta=Theta, Theta_nodes=Theta_nodes,
                                plot=plot)
        scores <- do.call(rbind, list_scores)
        if(plot){
            panel.bands <- function(x, y, upper, lower, fill, col,
                                       subscripts, ..., font, fontface){
                upper <- upper[subscripts]
                lower <- lower[subscripts]
                panel.polygon(c(x, rev(x)), c(upper, rev(lower)), col = fill, border = FALSE,
                              ...)
            }
            #for some reason the plot is drawn backwords....use rev() for now
            group <- rev(factor(rep(mod@Data$groupNames, each=nrow(Theta))))
            CIs <- apply(scores, 2L, bs_range, CI=CI)
            CIs <- CIs[-2L, ]
            df <- data.frame(Theta=rbind(Theta, Theta), group, TS=oCM, t(CIs))
            return(xyplot(TS ~ Theta, data=df, groups=group, auto.key=auto.key,
                   upper=df$upper, lower=df$lower, col=c('red', 'blue'),
                   fill=c('red', 'blue'), alpha=0.2,
                   panel = function(x, y, alpha, ...){
                       panel.superpose(x, y, panel.groups = panel.bands, type='l', alpha=alpha, ...)
                       panel.xyplot(x, y, type='l', lty=1,...)
                   },
                   xlab = expression(theta), ylab = expression(T(theta)),
                   main = 'Expected Total Score', ...))
        }
        if(!is.null(Theta_nodes)){
            CIs <- apply(scores, 2L, bs_range, CI=CI)
            rownames(CIs) <- rownames(CIs) <-
                c(paste0('CI_', round(CI + (1-CI)/2, 3L)*100), paste0('CI_', 50),
                  paste0('CI_', round((1-CI)/2, 3L)*100))
            return(cbind(Theta_nodes, "sDTF(Theta)"=oCM, t(CIs)))
        }
        t_sDTF <- oCM['signed.DTF'] / sd(scores[,'signed.DTF'])
        p_sDTF <- pt(abs(t_sDTF), df=MI-1L, lower.tail=FALSE) * 2
        CIs <- apply(scores, 2L, bs_range, CI=CI)
        if(!is.matrix(CIs)) stop('Too few MI draws were specified')
        tests <- c("P(sDTF = 0)" = as.numeric(p_sDTF))
        rownames(CIs) <- rownames(CIs) <-
            c(paste0('CI_', round(CI + (1-CI)/2, 3L)*100), paste0('CI_', 50),
              paste0('CI_', round((1-CI)/2, 3L)*100))
        ret <- list(observed=oCM, CIs=CIs, tests=tests)
    } else ret <- list_scores[[1L]]
    ret
}
