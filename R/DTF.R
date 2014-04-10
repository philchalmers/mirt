#' Differential test functioning
#'
#' Function performs various omnibus differential test functioning proceduces on an object
#' estimated with \code{multipleGroup()}. If the latent means/covariances are suspected to differ
#' then the input object should contain a set of 'anchor' items to ensure that only differential
#' test features are being detected rather than group differences. Returns signed (average area
#' above and below) and unsigned (total area) statistics, with descriptives such as the percent
#' average bias between group total scores.
#' 
#'
#' @aliases DTF
#' @param mod a multipleGroup object which estimated only 2 groups
#' @param MI a number indicating how many draws to take to form a suitable multiple imputation
#'   for the expected test scores (100 or more). Requires an estimated parameter
#'   information matrix
#' @param CI range of condfidince interval when using MI
#' @param npts number of points to use in the integration
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is 
#'   used in conjunction with \code{npts}
#' @param digits number of digits to round result to
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
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
#' plot(mod, type = 'score')
#'
#' DTF(mod)
#' DTF(mod, MI = 100) #95% C.I. for DTI containins 0
#' 
#' ## -------------
#' ## random slopes and intercepts for 15 items, and latent mean difference 
#' ##    (no systematic DTF should exist, but DIF will be present)
#' dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), sign(rnorm(n-15))*runif(n-15, .25, .5)),     
#'                 d + c(numeric(15), sign(rnorm(n-15))*runif(n-15, .5, 1)), N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod1 <- multipleGroup(dat, 1, group=group)
#' plot(mod1, type = 'score') #does not account for group differences! Need anchors
#' 
#' mod2 <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type = 'crossprod',
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod2, type = 'score')
#' 
#' mirtCluster()
#' #significant DIF in multiple items....
#' DIF(mod2, which.par=c('a1', 'd'), items2test=16:30) 
#' DTF(mod2) #...but not substantial DTF due to randomness of DIF
#' DTF(mod2, MI=100)
#' 
#' ## -------------
#' ## systematic differing slopes and intercepts (clear DTF)
#' dat1 <- simdata(a, d, N, itemtype = 'dich', mu=.50, sigma=matrix(1.5))
#' dat2 <- simdata(a + c(numeric(15), rnorm(n-15, 1, .25)), d + c(numeric(15), rnorm(n-15, 1, .5)),
#'                 N, itemtype = 'dich')
#' dat <- rbind(dat1, dat2)
#' mod3 <- multipleGroup(dat, model, group=group, SE=TRUE, SE.type='crossprod',
#'                       invariance=c('free_means', 'free_var'))
#' plot(mod3, type = 'score') #visable DTF happening
#' 
#' DIF(mod3, c('a1', 'd'), items2test=16:30) 
#' DTF(mod3) #huge unsigned bias. Signed bias indicates group 2 scores generally lower
#' DTF(mod3, MI=100) 
#' 
#' }
DTF <- function(mod, MI = NULL, CI = .95, npts = 200, digits = 4, theta_lim=c(-6,6)){

    if(class(mod) != 'MultipleGroupClass')
        stop('mod input was not estimated by multipleGroup()')
    if(length(mod@cmods) != 2L)
        stop('DTF only supports two group models at a time')
    J <- length(mod@K)
    if(is.null(MI)){
        MI <- 1L
        impute <- FALSE
    } else {
        if(is(try(chol(mod@information), silent=TRUE), 'try-error')){
            stop('Proper information matrix must be precomputed in model')
        } else {
            impute <- TRUE
            list_scores <- vector('list', MI)
            mod <- assignInformationMG(mod)
            covBs <- lapply(mod@cmods, function(x) solve(x@information))
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
    Theta <- thetaComb(theta, mod@nfact)
    max_score <- sum(apply(mod@data, 2, min) + (mod@K - 1L))
    omod <- mod
    for(mi in 1L:MI){
        if(impute){
            for(g in 1L:2L)
                mod@cmods[[g]]@pars <- imputePars(pars=omod@cmods[[g]]@pars, covB=covBs[[g]],
                                                  imputenums=imputenums[[g]], constrain=omod@constrain)
        }
        T1 <- expected.test(mod, Theta, group=1)
        T2 <- expected.test(mod, Theta, group=2)
        D <- T1 - T2
        uDTF <- mean(D^2)
        uDTF_percent <- sqrt(uDTF)/max_score * 100
        sDTF <- mean(D)
        sDTF_percent <- sDTF/max_score * 100
        max_DTF_percent <- max(abs(D))/max_score * 100
        ret <- list(signed = c(DTF=sDTF, `DTF(%)`=sDTF_percent),
                    unsigned = c(DTF=uDTF, `DTF(%)`=uDTF_percent,
                                 `max.DTF(%)`=max_DTF_percent))
        if(impute) list_scores[[mi]] <- c(ret$signed, ret$unsigned)
    }
    if(impute){
        scores <- do.call(rbind, list_scores)
        CM <- apply(scores, 2, mean)
        SD <- apply(scores, 2, sd)
        t_sDTF <- CM[1L] / SD[1L]
        p_sDTF <- pt(abs(t_sDTF), df=MI-1, lower.tail=FALSE)
        t_DTF <- CM[3L] / SD[3L]
        p_DTF <- pt(t_DTF, df=MI-1, lower.tail=FALSE) * 2 #one-tailed
        tt <- qt(CI + (1-CI)/2, df=MI-1)
        upper <- CM + tt * SD
        lower <- CM - tt * SD
        signed <- rbind(upper[1:2], CM[1:2], lower[1:2])
        unsigned <- rbind(upper[3:5], CM[3:5], lower[3:5])
        tests <- c(p_sDTF, p_DTF)
        names(tests) <- c("P(sDTF = 0)", "P(uDTF = 0)")
        rownames(signed) <- rownames(unsigned) <-
            c(paste0('CI_', round(CI + (1-CI)/2,3)), 'value',
              paste0('CI_', round((1-CI)/2, 3)))
        ret <- list(signed=signed, unsigned=unsigned, tests=tests)
    }
    ret <- lapply(ret, round, digits=digits)
    ret
}