#' Collapse values from multiple imputation draws
#'
#' This function computes updated parameter and standard error estimates using multiple
#' imputation methodology. Given a set of parameter estimates and their associated standard
#' errors the function returns the weighted average of the overall between and within
#' variability due to the multiple imputations according to Rubin's (1987) methodology.
#'
#' @aliases averageMI
#' @param par a list containing parameter estimates which were computed the imputed datasets
#' @param SEpar a list containing standard errors associated with \code{par}
#' @param as.data.frame logical; return a data.frame instead of a list? Default is TRUE
#' @param digits number of digits to round result. Default is 4
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @return returns a list or data.frame containing the updated averaged parameter estimates,
#'   standard errors, and t-values with the associated degrees of freedom and two tailed p-values
#' @keywords multiple imputation
#' @references
#' Rubin, D.B. (1987) Multiple Imputation for Nonresponse in Surveys. Wiley & Sons, New York.
#' @export averageMI
#' @examples
#'
#' \dontrun{
#'
#' #simulate data
#' set.seed(1234)
#' N <- 1000
#'
#' # covariates
#' X1 <- rnorm(N); X2 <- rnorm(N)
#' covdata <- data.frame(X1, X2)
#' Theta <- matrix(0.5 * X1 + -1 * X2 + rnorm(N, sd = 0.5))
#'
#' #items and response data
#' a <- matrix(1, 20); d <- matrix(rnorm(20))
#' dat <- simdata(a, d, 1000, itemtype = 'dich', Theta=Theta)
#'
#' mod1 <- mirt(dat, 1, 'Rasch', covdata=covdata, formula = ~ X1 + X2)
#' coef(mod1, simplify=TRUE)
#'
#' #draw plausible values for secondary analyses
#' pv <- fscores(mod1, plausible.draws = 10)
#' pvmods <- lapply(pv, function(x, covdata) lm(x ~ covdata$X1 + covdata$X2),
#'                  covdata=covdata)
#'
#' # compute Rubin's multiple imputation average
#' so <- lapply(pvmods, summary)
#' par <- lapply(so, function(x) x$coefficients[, 'Estimate'])
#' SEpar <- lapply(so, function(x) x$coefficients[, 'Std. Error'])
#' averageMI(par, SEpar)
#'
#' }
averageMI <- function(par, SEpar, as.data.frame = TRUE, digits = 4){
    if(missing(par)) missingMsg('par')
    if(missing(SEpar)) missingMsg('SEpar')
    if(!is.list(par)) stop('par must be a list', call.=FALSE)
    if(!is.list(SEpar)) stop('SEpar must be a list', call.=FALSE)
    par <- lapply(par, as.matrix)
    SEpar <- lapply(SEpar, as.matrix)
    MI <- length(par)
    scores <- par[[1L]]
    Ubar <- SEpar[[1L]]^2
    for(i in 2L:MI){
        scores <- par[[i]] + scores
        Ubar <- SEpar[[i]]^2 + Ubar
    }
    scores <- scores / MI
    Ubar <- Ubar / MI
    tmp <- lapply(par, function(x, scores, MI)
        (x - scores)^2, scores=scores, MI=MI)
    B <- tmp[[1L]]
    for(i in 2L:MI) B <- B + tmp[[i]]
    B <- (1 / (MI-1L)) * B
    SEscores <- sqrt(Ubar + (1 + 1/MI) * B)
    df <- (MI - 1) * (1 + MI * Ubar / ((MI + 1) * B))^2
    ret <- list(par=scores, SEpar=SEscores, t = scores/SEscores, df=df)
    ret$p <- (1 - pt(abs(ret$t), ret$df, lower.tail=TRUE))/2
    ret <- lapply(ret, function(x, digits) round(x, digits), digits=digits)
    if(as.data.frame){
        n <- ncol(scores)
        ret <- as.data.frame(ret)
        if(n == 1)
            colnames(ret) <- c('par', 'SEpar', 't', 'df', 'p')
        else colnames(ret) <- paste(c('par', 'SEpar', 't', 'df', 'p'), 1:n, sep='_')
    }
    return(ret)
}
