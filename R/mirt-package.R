#' Full information maximum likelihood estimation of multidimensional IRT models
#'
#' Analysis of dichotomous and polytomous response data using unidimensional and
#' multidimensional latent trait models under the Item Response Theory paradigm.
#' Exploratory and confirmatory models can be estimated with quadrature (EM) or
#' stochastic (MHRM) methods. Confirmatory bi-factor and two-tier analyses are available
#' for modeling item testlets. Multiple group analysis and mixed effects designs also
#' are available for detecting differential item functioning and modeling item and
#' person covariates.
#'
#' Users interested in the most recent version of this package can visit
#' \url{https://github.com/philchalmers/mirt} and follow the instructions
#' for installing the package from source. Questions regarding the package can
#' be sent to the mirt-package Google Group, located at
#' \url{https://groups.google.com/forum/#!forum/mirt-package}.
#'
#'
#'
#' @name mirt-package
#' @docType package
#' @title Full information maximum likelihood estimation of IRT models.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @useDynLib mirt
#' @importFrom stats anova fitted residuals
#' @importFrom MASS ginv
#' @import lattice GPArotation mvtnorm Rcpp stats4 methods
#' @exportMethod anova
#' @exportMethod fitted
#' @exportMethod residuals
#' @exportMethod summary
#' @keywords package
NULL

#' Description of Science data
#'
#' A 4-item data set borrowed from \code{ltm} package in R, first example
#' of the \code{grm()} function. See more complete documentation therein.
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
#' respectively. However, careful analysis using the nominal response model suggests that the
#' scoring key for item 32 may be incorrect, and should be changed from 5 to 3.
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
#' @examples
#'
#' \dontrun{
#' #score the data (missing scored as 0)
#' head(SAT12)
#' data <- key2binary(SAT12,
#'     key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' #score the data, missing treated as NA
#' SAT12missing <- SAT12
#' SAT12missing[SAT12missing == '8'] <- NA
#' data <- key2binary(SAT12missing,
#'     key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#'
#' #potentially better scoring for item 32
#' data <- key2binary(SAT12,
#'     key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,3))
#' }
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
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT7)
#' head(dat)
#' (mod <- mirt(dat, 1))
#' coef(mod)
#' }
NULL

#' Description of LSAT6 data
#'
#' Data from Thissen (1982); contains 5 dichotomously scored
#' items obtained from the Law School Admissions Test, section 6.
#'
#'
#' @name LSAT6
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
#' \emph{Psychometrika, 47}, 175-186.
#'
#' @keywords data
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT6)
#' head(dat)
#' model <- mirt.model('F = 1-5
#'                      CONSTRAIN = (1-5, a1)')
#' (mod <- mirt(dat, model))
#' coef(mod)
#'
#'
#' }
NULL

#' Description of deAyala data
#'
#' Mathematics data from de Ayala (2009; pg. 14); 5 item dataset in table format.
#'
#'
#' @name deAyala
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' de Ayala, R. J. (2009). \emph{The theory and practice of item response theory}. Guilford Press.
#'
#' @keywords data
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(deAyala)
#' head(dat)
#' }
NULL

#' Description of Bock 1997 data
#'
#' A 3-item tabulated data set extracted from Table 3 in Chapter Two.
#'
#'
#' @name Bock1997
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords data
#' @references
#' Bock, R. D. (1997). The Nominal Categories Model. In van der Linden, W. J. & Hambleton, R. K.
#' \emph{Handbook of modern item response theory}. New York: Springer.
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(Bock1997)
#' head(dat)
#' mod <- mirt(dat, 1, 'nominal')
#'
#' #reproduce table 3 in Bock (1997)
#' fs <- round(fscores(mod, verbose = FALSE)[,c('F1','SE_F1')],2)
#' fttd <- round(fitted(mod),1)
#' table <- data.frame(fttd, fs)
#' table
#' 
#' #using nominal.highlow matrix to specify lowest and highest categories
#' (nominal.highlow <- matrix(c(4,4,4,4,1,1,1,1), 2, byrow = TRUE))
#' mod <- mirt(dat, 1, 'nominal', nominal.highlow=nominal.highlow)
#' coef(mod)
#'
#'  }
NULL

