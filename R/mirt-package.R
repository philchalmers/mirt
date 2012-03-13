#' Full information maximum likelihood estimation of multidimensional IRT models
#' 
#' Analysis of dichotomous and polychotomous response data using
#' latent trait models under the Item Response Theory paradigm. Includes the
#' multivariate two- and three-parameter logistic models, confirmatory
#' bifactor analysis, polytomous confirmatory and exploratory item response
#' models, and partially-compensatory item response modeling in conjunction
#' with other IRT models.
#'
#' 
#' @name mirt-package
#' @docType package
#' @title Full information maximum likelihood estimation of IRT models.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @useDynLib mirt 
#' @import psych
#' @import MASS
#' @import GPArotation
#' @import mvtnorm
#' @import Matrix
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

#' Description of Science data
#' 
#' A 4 item data set borrowed from \code{\link[ltm]{ltm}} package, first example 
#' of the \code{\link[ltm]{grm}} function. See more complete documentation therein.  
#' 
#' 
#' @name Science
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
NULL

#' Description of SAT12 data
#' 
#' Data obtained from the TESTFACT (Woods et al., 2003) manual, with 32 response pattern 
#' scored items for a grade 12 science assessment test (SAT) measuring topics of chemistry, 
#' biology, and physics. The scoring key for these data is 
#' [1, 4, 5, 2, 3, 1, 2, 1, 3, 1, 2, 4, 2, 1, 5, 3, 4, 4, 1, 4, 3, 3, 4, 1, 3, 5, 1, 3, 1, 5, 4, 5], 
#' respectively. 
#' 
#' 
#' @name SAT12
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003).
#' TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis 
#' [Computer software]. Lincolnwood, IL: Scientific Software International.
#'
#' @keywords data
NULL

#' Description of LSAT7 data
#' 
#' Data from Bock & Lieberman (1970); contains 5 dichotomously scored 
#' items obtained from the Law School Admissions Test, section 7. 
#' 
#' 
#' @name LSAT7
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Bock, R. D., & Lieberman, M. (1970). Fitting a response model for \emph{n} 
#' dichotomously scored items. \emph{Psychometrika, 35}(2), 179-197.
#'
#' @keywords data
NULL
