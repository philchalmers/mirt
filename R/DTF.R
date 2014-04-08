#' Differential test functioning
#'
#' Function performs various omnibus differential test functioning proceduces on an object
#' estimated with \code{multipleGroup()}. If the latent means/covariances are suspected to differ
#' then the input object should contain a set of 'anchor' items to ensure that only differential 
#' test features are being detected rather than group differences.
#' 
#' @aliases DTF
#' @param mod a multipleGroup object
#' @param MI a number indicating how many draws to take to form a suitable multiple imputation
#'   for the expected test scores (100 or more). Requires an estimated parameter 
#'   information matrix
#' @param CI range of condfidince interval when using MI
#' @param npts number of points to use in the integration, ranging from -6 to 6
#' @param digits number of digits to round result to
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{multipleGroup}}, \code{\link{DIF}}
#' @keywords DTF
#' @export DTF
#' @examples
#' \dontrun{
#' 
#' 
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
        uDTF <- mean(abs(D))
        uDTF_percent <- uDTF/max_score * 100
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
        tt <- qt(CI + (1-CI)/2, df=MI-1)
        upper <- CM + tt * SD
        lower <- CM - tt * SD
        signed <- rbind(upper[1:2], CM[1:2], lower[1:2])
        unsigned <- rbind(upper[3:5], CM[3:5], lower[3:5])
        rownames(signed) <- rownames(unsigned) <- 
            c(paste0('CI_', round(CI + (1-CI)/2,3)), 'value', 
              paste0('CI_', round((1-CI)/2, 3)))
        ret <- list(signed=signed, unsigned=unsigned)
    }
    ret <- lapply(ret, round, digits=digits)
    ret
}