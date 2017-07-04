#' Full information maximum likelihood estimation of multidimensional IRT models
#'
#' Analysis of dichotomous and polytomous response data using
#' unidimensional and multidimensional latent trait models under the Item
#' Response Theory paradigm. Exploratory and confirmatory models can be
#' estimated with quadrature (EM) or stochastic (MHRM) methods. Confirmatory
#' bi-factor and two-tier analyses are available for modeling item testlets.
#' Multiple group analysis and mixed effects designs also are available for
#' detecting differential item and test functioning as well as modelling
#' item and person covariates. Finally, latent class models such as the DINA,
#' DINO, multidimensional latent class, and several other discrete variable
#' models are supported.
#'
#' Users interested in the most recent version of this package can visit
#' \url{https://github.com/philchalmers/mirt} and follow the instructions
#' for installing the package from source. Questions regarding the package can
#' be sent to the mirt-package Google Group, located at
#' \url{https://groups.google.com/forum/#!forum/mirt-package}. User contributed files,
#' workshop files, and evaluated help files are also available on the package wiki
#' (\url{https://github.com/philchalmers/mirt/wiki}).
#'
#'
#'
#' @name mirt-package
#' @docType package
#' @title Full information maximum likelihood estimation of IRT models.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @useDynLib mirt
#' @import stats lattice GPArotation Rcpp stats4 methods sfsmisc mgcv splines
#' @importFrom utils write.table flush.console packageVersion
#' @importFrom graphics symbols
#' @exportMethod anova residuals summary logLik vcov
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
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
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @keywords data
#' @examples
#'
#' \dontrun{
#' mod <- mirt(Science, 1)
#' plot(mod, type = 'trace')
#' }
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
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
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
#' head(data)
#'
#' #score the data, missing (value of 8) treated as NA
#' SAT12missing <- SAT12
#' SAT12missing[SAT12missing == 8] <- NA
#' data <- key2binary(SAT12missing,
#'     key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' head(data)
#'
#' #potentially better scoring for item 32 (based on nominal model finding)
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
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
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
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Thissen, D. (1982). Marginal maximum likelihood estimation for the one-parameter logistic model.
#' \emph{Psychometrika, 47}, 175-186.
#'
#' @keywords data
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(LSAT6)
#' head(dat)
#' model <- 'F = 1-5
#'          CONSTRAIN = (1-5, a1)'
#' (mod <- mirt(dat, model))
#' M2(mod)
#' itemfit(mod)
#' coef(mod, simplify=TRUE)
#'
#' #equivalentely, but with a different parameterization
#' mod2 <- mirt(dat, 1, itemtype = 'Rasch')
#' anova(mod, mod2) #equal
#' M2(mod2)
#' itemfit(mod2)
#' coef(mod2, simplify=TRUE)
#' sqrt(coef(mod2)$GroupPars[2]) #latent SD equal to the slope in mod
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
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
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
#' 
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(Bock1997)
#' head(dat)
#' mod <- mirt(dat, 1, 'nominal')
#'
#' #reproduce table 3 in Bock (1997)
#' fs <- round(fscores(mod, verbose = FALSE, full.scores = FALSE)[,c('F1','SE_F1')],2)
#' fttd <- residuals(mod, type = 'exp')
#' table <- data.frame(fttd[,-ncol(fttd)], fs)
#' table
#'
#' mod <- mirt(dat, 1, 'nominal')
#' coef(mod)
#'
#'  }
NULL

