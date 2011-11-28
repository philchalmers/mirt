#' Full information maximum likelihood estimation of mutlidimensional IRT models
#' 
#' Analysis of dichotomous and polychotomous response data using
#' latent trait models under the Item Response Theory paradigm. Includes the
#' multivariate two- and three-parameter logistic models, confirmatory
#' bifactor analysis, polytomous confirmatory and exploratory item response
#' models, and partially-compensatory item response modeling in conjunction
#' with other IRT models.
#' ...
#' @name mirt-package
#' @docType package
#' @title Full information maximum likelihood estimation of IRT models.
#' @author \email{rphilip.chalmers@@gmail.com}
#' @useDynLib mirt 
#' @importFrom graphics plot
#' @importFrom stats anova coef fitted residuals logLik
#' @exportMethod anova 
#' @exportMethod coef
#' @exportMethod fitted
#' @exportMethod residuals
#' @exportMethod logLik
#' @exportMethod plot
#' @keywords package
NULL

