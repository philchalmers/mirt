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
#' @param draws a number indicating how many draws to take to form a suitable multiple imputation
#'   estimate of the expected test scores (usually 100 or more). Returns a list containing the
#'   imputation distribution and null hypothesis test for the sDTF statistic
#' @param CI range of confidence interval when using draws input
#' @param npts number of points to use in the integration. Default is 1000
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param Theta_nodes an optional matrix of Theta values to be evaluated in the draws for the
#'   sDTF statistic. However, these values are not averaged across, and instead give the bootstrap
#'   confidence intervals at the respective Theta nodes. Useful when following up a large
#'   uDTF/sDTF statistic to determine where the difference between the test curves are large
#'   (while still accounting for sampling variability). Returns a matrix with observed
#'   variability
#' @param plot a character vector indicating which plot to draw. Possible values are 'none',
#'    'func' for the test score functions, and 'sDTF' for the evaluated sDTF values across the
#'    integration grid. Each plot is drawn with imputed confidence envelopes
#' @param auto.key logical; automatically generate key in lattice plot?
#' @param ... additional arguments to be passed to \code{lattice} and \code{boot}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chalmers, R. P., Counsell, A., and Flora, D. B. (2016). It might not
#'   make a big DIF: Improved Differential Test Functioning statistics that account for
#'   sampling variability. \emph{Educational and Psychological Measurement, 76}, 114-140.
#'   \doi{10.1177/0013164415584576}
#'
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
#' model <- 'F = 1-30
#'           CONSTRAINB = (1-5, a1), (1-5, d)'
#'
#' a <- matrix(1, n)
#' d <- matrix(rnorm(n), n)
#' group <- c(rep('Group_1', N), rep('Group_2', N))
#'
#' ## -------------
#' # groups completely equal
#' dat1 <- simdata(a, d, N, itemtype = '2PL')
#' dat2 <- simdata(a, d, N, itemtype = '2PL')
#' dat <- rbind(dat1, dat2)
#' mod <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                      invariance=c('free_means', 'free_var'))
#' plot(mod)
#'
#' DTF(mod)
#' mirtCluster()
#' DTF(mod, draws = 1000) #95% C.I. for sDTF containing 0. uDTF is very small
#' DTF(mod, draws = 1000, plot='sDTF') #sDTF 95% C.I.'s across Theta always include 0
#'
#' ## -------------
#' ## random slopes and intercepts for 15 items, and latent mean difference
#' ##    (no systematic DTF should exist, but DIF will be present)
#' set.seed(1234)
#' dat1 <- simdata(a, d, N, itemtype = '2PL', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), runif(n-15, -.2, .2)),
#'                 d + c(numeric(15), runif(n-15, -.5, .5)), N, itemtype = '2PL')
#' dat <- rbind(dat1, dat2)
#' mod1 <- multipleGroup(dat, 1, group=group)
#' plot(mod1) #does not account for group differences! Need anchors
#'
#' mod2 <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod2)
#'
#' #significant DIF in multiple items....
#' # DIF(mod2, which.par=c('a1', 'd'), items2test=16:30)
#' DTF(mod2)
#' DTF(mod2, draws=1000) #non-sig DTF due to item cancellation
#'
#' ## -------------
#' ## systematic differing slopes and intercepts (clear DTF)
#' dat1 <- simdata(a, d, N, itemtype = '2PL', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, .5)),
#'                 N, itemtype = '2PL')
#' dat <- rbind(dat1, dat2)
#' mod3 <- multipleGroup(dat, model, group=group, SE=TRUE,
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod3) #visable DTF happening
#'
#' # DIF(mod3, c('a1', 'd'), items2test=16:30)
#' DTF(mod3) #unsigned bias. Signed bias indicates group 2 scores generally higher on average
#' DTF(mod3, draws=1000)
#' DTF(mod3, draws=1000, plot='func')
#' DTF(mod3, draws=1000, plot='sDTF') #multiple DTF areas along Theta
#'
#' # evaluate specific values for sDTF
#' Theta_nodes <- matrix(seq(-6,6,length.out = 100))
#' sDTF <- DTF(mod3, Theta_nodes=Theta_nodes)
#' head(sDTF)
#' sDTF <- DTF(mod3, Theta_nodes=Theta_nodes, draws=100)
#' head(sDTF)
#'
#' }
DTF <- function(mod, draws = NULL, CI = .95, npts = 1000, theta_lim=c(-6,6), Theta_nodes = NULL,
                plot = 'none', auto.key = list(space = 'right', points=FALSE, lines=TRUE), ...){

    bs_range <- function(x, CI){
        ss <- sort(x)
        N <- length(ss)
        ret <- c(upper = ss[ceiling(N * (1 - (1-CI)/2))],
                 lower = ss[floor(N * (1-CI)/2)])
        ret
    }

    panel.bands <- function(x, y, upper, lower, fill, col,
                            subscripts, ..., font, fontface){
        upper <- upper[subscripts]
        lower <- lower[subscripts]
        panel.polygon(c(x, rev(x)), c(upper, rev(lower)), col = fill, border = FALSE,
                      ...)
    }

    fn <- function(x, omod, impute, imputenums, Theta, max_score, Theta_nodes = NULL,
                   plot, integration, theta_lim, type, pre.ev, shortpars, longpars){
        mod <- omod
        if(impute){
            mod <- imputePars2(MGmod=mod, pre.ev=pre.ev, imputenums=imputenums,
                               shortpars=shortpars, longpars=longpars)
        }
        if(!is.null(Theta_nodes)){
            if(type == 'score'){
                T1 <- expected.test(mod, Theta_nodes, group=1L, mins=FALSE)
                T2 <- expected.test(mod, Theta_nodes, group=2L, mins=FALSE)
            }
            D <- T1 - T2
            ret <- c("sDTF." = D)
            return(ret)
        }
        calc_DTFs(mod=mod, Theta=Theta, plot=plot, max_score=max_score, type=type)
    }

    if(missing(mod)) missingMsg('mod')
    stopifnot(is.character(plot))
    stopifnot(mod@Model$nfact == 1L)
    boot <- FALSE
    integration <- 'quad'
    type <- 'score'
    if(!(plot %in% c('none', 'func', 'sDTF')))
        stop('plot type not supported')
    if(!is.null(Theta_nodes)){
        integration <- 'quad'
        if(plot != 'none') message('plots are not drawn when Theta_nodes is included')
        plot <- 'none'
    }
    if(class(mod) != 'MultipleGroupClass')
        stop('mod input was not estimated by multipleGroup()', call.=FALSE)
    if(length(mod@ParObjects$pars) != 2L)
        stop('DTF only supports two group models at a time', call.=FALSE)
    if(!any(sapply(mod@ParObjects$pars, function(x) x@ParObjects$pars[[length(x@ParObjects$pars)]]@est)))
        message('No hyper-parameters were estimated in the DIF model. For effective
                \tDTF testing, freeing the focal group hyper-parameters is recommend.')
    if(!is.null(Theta_nodes)){
        if(!is.matrix(Theta_nodes))
            stop('Theta_nodes must be a matrix', call.=FALSE)
        if(ncol(Theta_nodes) != mod@Model$nfact)
            stop('Theta_nodes input does not have the correct number of factors', call.=FALSE)
        colnames(Theta_nodes) <- if(ncol(Theta_nodes) > 1L)
            paste0('Theta.', 1L:ncol(Theta_nodes)) else 'Theta'
    }
    if(plot != 'none'){
        if(is.null(draws))
            stop('Must specify number of draws to generate plot confidence intervals', call.=FALSE)
    }
    if(length(type) > 1L && (plot != 'none' || !is.null(Theta_nodes)))
        stop('Multiple type arguments cannot be combined with plot or Theta_nodes arguments')

    if(is.null(draws)){
        draws <- 1L
        impute <- FALSE
    } else if(!boot){
        if(length(mod@vcov) == 1L)
            stop('Stop an information matrix must be computed', call.=FALSE)
        if(!mod@OptimInfo$secondordertest)
            stop('ACOV matrix is not positive definite')
        impute <- TRUE
        shortpars <- mod@Internals$shortpars
        covB <- mod@vcov
        names <- colnames(covB)
        imputenums <- sapply(strsplit(names, '\\.'), function(x) as.integer(x[2L]))
        longpars <- c(do.call(c, lapply(mod@ParObjects$pars[[1L]]@ParObjects$pars, function(x) x@par)),
                      do.call(c, lapply(mod@ParObjects$pars[[2L]]@ParObjects$pars, function(x) x@par)))
        pre.ev <- eigen(covB)
    } else impute <- TRUE
    if(is.null(Theta_nodes)){
        if(integration == 'quad'){
            theta <- matrix(seq(theta_lim[1L], theta_lim[2L], length.out=npts))
            Theta <- thetaComb(theta, mod@Model$nfact)
        }
        if(plot == 'sDTF') Theta_nodes <- Theta
    } else Theta <- Theta_nodes
    max_score <- sum(mod@Data$mins + mod@Data$K - 1L)
    list_scores <- myLapply(1L, fn, omod=mod, impute=FALSE, Theta_nodes=Theta_nodes,
                            imputenums=NULL, max_score=max_score, Theta=Theta, plot=plot,
                            integration=integration, theta_lim=theta_lim, type=type)
    if(impute){
        oCM <- list_scores[[1L]]
        list_scores <- myLapply(1L:draws, fn, omod=mod, impute=TRUE, shortpars=shortpars,
                                max_score=max_score, imputenums=imputenums, Theta=Theta,
                                Theta_nodes=Theta_nodes, plot=plot, integration=integration,
                                theta_lim=theta_lim, type=type, pre.ev=pre.ev, longpars=longpars)
        scores <- do.call(rbind, list_scores)
        pars <- list(mod@ParObjects$pars[[1L]]@ParObjects$pars, mod@ParObjects$pars[[2L]]@ParObjects$pars)
        pars <- reloadPars(longpars=longpars, pars=pars, ngroups=2L, J=length(pars[[1L]])-1L)
        if(!is.null(Theta_nodes)){
            CIs <- apply(scores, 2L, bs_range, CI=CI)
            rownames(CIs) <- rownames(CIs) <-
                c(paste0('CI_', round(CI + (1-CI)/2, 3L)*100),
                  paste0('CI_', round((1-CI)/2, 3L)*100))
            ret <- data.frame(Theta=Theta_nodes, "sDTF"=oCM, t(CIs))
            rownames(ret) <- paste(type, 1L:nrow(ret), sep='.')
            if(plot == 'sDTF'){
                lim <- c(min(ret[,4L]), max(ret[,3L]))
                return(xyplot(ret[,2L] ~ ret[,1L], upper=ret[,3L], lower = ret[,4L], alpha = .2,
                              fill = 'darkgrey',
                              panel = function(x, y, lower, upper, fill, alpha, ...){
                                  panel.xyplot(c(min(x), max(x)), c(0,0), col = 'red', type = 'l')
                                  panel.polygon(c(x, rev(x)), c(upper, rev(lower)), col = fill,
                                                border = FALSE, alpha=alpha, ...)
                                  panel.xyplot(x, y, type='l', lty=1, col = 'black', ...)
                                  },
                              ylab = expression(sDTF(theta)), xlab = expression(theta), ylim = lim,
                              main = 'Signed DTF', ...))
            } else {
                return(ret)
            }
        }
        if(plot == 'func'){
            group <- rep(mod@Data$groupNames, each=nrow(Theta))
            CIs <- apply(scores, 2L, bs_range, CI=CI)
            df <- data.frame(Theta=rbind(Theta, Theta), group, TS=oCM, t(CIs))
            lim <- c(min(df$lower), max(df$upper))
            main <- switch(type, score='Expected Total Score')
            ylab <- switch(type, score=expression(T(theta)))
            return(xyplot(TS ~ Theta, data=df, groups=group, auto.key=auto.key,
                   upper=df$upper, lower=df$lower, ylim = lim,
                   panel = function(x, y, ...){
                       panel.superpose(x, y, panel.groups = panel.bands, type='l',  ...)
                       panel.xyplot(x, y, type='l', lty=1,...)
                   },
                   xlab = expression(theta), ylab = ylab,
                   main = main, ...))
        }
        tests <- NULL
        if('score' %in% type){
            ts <- oCM['sDTF.score'] / sd(scores[,'sDTF.score'])
            ps <- pt(abs(ts), df = draws-1, lower.tail = FALSE) * 2
            tests <- c('P(sDTF.score = 0)'=as.vector(ps))
        }
        CIs <- apply(scores, 2L, bs_range, CI=CI)
        if(!is.matrix(CIs) || nrow(CIs) != 2L)
            stop('Too few draws were specified', call.=FALSE)
        rownames(CIs) <-
            c(paste0('CI_', round(CI + (1-CI)/2, 3L)*100),
              paste0('CI_', round((1-CI)/2, 3L)*100))
        ret <- list(observed=oCM, CIs=CIs, tests=tests)
        if(is.null(ret$tests)) ret$tests <- NULL
    } else {
        if(!is.null(Theta_nodes)){
            ret <- data.frame(Theta_nodes, sDTF=list_scores[[1L]])
            rownames(ret) <- paste(type, 1L:nrow(ret), sep='.')
        } else {
            ret <- list_scores[[1L]]
        }
    }
    ret
}

calc_DTFs <- function(mod, Theta, plot, max_score, type){
    Ds <- matrix(NA, nrow(Theta), 2L)
    colnames(Ds) <- c('score', 'null')
    if('score' %in% type){
        T1 <- expected.test(mod, Theta, group=1L, mins=FALSE)
        T2 <- expected.test(mod, Theta, group=2L, mins=FALSE)
        Ds[,1L] <- T1 - T2
    }
    if(plot != 'none') return(c(T1, T2))
    uDTF <- colMeans(abs(Ds))
    uDTF_percent <- uDTF[1L]/max_score * 100
    sDTF <- colMeans(Ds)
    sDTF_percent <- sDTF[1L]/max_score * 100
    ret <- c("sDTF"=sDTF[1L], "sDTF(%)"=sDTF_percent, "sDTF"=sDTF[-1L],
             "uDTF"=uDTF[1L], "uDTF(%)"=uDTF_percent, "uDTF"=uDTF[-1L])
    ret[!is.na(ret)]
}
