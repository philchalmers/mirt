#' Full information maximum likelihood estimation of multidimensional IRT models
#'
#' Analysis of dichotomous and polytomous response data using
#' unidimensional and multidimensional latent trait models under the Item
#' Response Theory (IRT) paradigm. Exploratory and confirmatory models can be
#' estimated with quadrature (EM) or stochastic (MHRM) methods. Confirmatory
#' bi-factor and two-tier analyses are available for modeling item testlets.
#' Multiple group analysis and mixed effects designs also are available for
#' detecting differential item and test functioning as well as modeling
#' item and person covariates. Finally, latent class models such as the DINA,
#' DINO, multidimensional latent class, mixture and zero-inflated IRT models,
#' and several other discrete variable models are supported.
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
#' @title Full information maximum likelihood estimation of IRT models.
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @useDynLib mirt
#' @import stats lattice GPArotation Rcpp stats4 methods mgcv splines vegan dcurver pbapply
#' @importFrom utils write.table flush.console packageVersion capture.output head
#' @importFrom gridExtra grid.arrange
#' @importFrom graphics symbols
#' @importFrom Deriv Deriv
#' @importFrom SimDesign manageWarnings
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
#' of the \code{grm()} function. See more complete documentation therein, as
#' well as Karlheinz and Melich (1992).
#'
#'
#' @name Science
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#'
#' Karlheinz, R. and Melich, A. (1992). Euro-Barometer 38.1:
#' \emph{Consumer Protection and Perceptions of Science and Technology}.
#' INRA (Europe), Brussels. [computer file]
#'
#' @keywords data
#' @examples
#'
#' \dontrun{
#' itemstats(Science)
#'
#' mod <- mirt(Science, 1)
#' plot(mod, type = 'trace')
#' }
NULL

#' Social Life Feelings Data
#'
#' A 5-item data set analyzed by Bartholomew (1998). Data contains
#' dichotomous responses (endorsement vs non-endorsement) from 1490 German
#' respondents to five statements on perceptions of social life.
#'
#' @name SLF
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Bartholomew, D., J. (1998). Scaling unobservable constructs in social science. Journal of the Royal
#' Statistical Society - Series C, 47, 1-13.
#' @keywords data
#' @examples
#'
#' \dontrun{
#' # tabular format
#' data(SLF)
#' SLF
#'
#' # full dataset
#' full <- expand.table(SLF)
#' itemstats(full)
#'
#' mod <- mirt(full)
#' plot(mod, type = 'trace')
#'
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
#'
#' Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003).
#' TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis
#' [Computer software]. Lincolnwood, IL: Scientific Software International.
#'
#' @keywords data
#' @examples
#'
#' \dontrun{
#'
#' itemstats(SAT12, use_ts = FALSE)
#'
#' # score the data (missing scored as 0)
#' head(SAT12)
#' dat <- key2binary(SAT12,
#'     key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' head(dat)
#' itemstats(dat)
#'
#' # score the data, missing (value of 8) treated as NA
#' SAT12missing <- SAT12
#' SAT12missing[SAT12missing == 8] <- NA
#' dat <- key2binary(SAT12missing,
#'     key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' head(dat)
#'
#' # potentially better scoring for item 32 (based on nominal model finding)
#' dat <- key2binary(SAT12,
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
#' itemstats(dat)
#'
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
#' itemstats(dat)
#'
#' model <- 'F = 1-5
#'          CONSTRAIN = (1-5, a1)'
#' (mod <- mirt(dat, model))
#' M2(mod)
#' itemfit(mod)
#' coef(mod, simplify=TRUE)
#'
#' # equivalentely, but with a different parameterization
#' mod2 <- mirt(dat, 1, itemtype = 'Rasch')
#' anova(mod, mod2) #equal
#' M2(mod2)
#' itemfit(mod2)
#' coef(mod2, simplify=TRUE)
#' sqrt(coef(mod2)$GroupPars[2]) #latent SD equal to the slope in mod
#'
#' }
NULL

#' Description of ASVAB data
#'
#' Data from
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
#' itemstats(dat)
#'
#' (mod <- mirt(dat, 1))
#' coef(mod)
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
#'
#' de Ayala, R. J. (2009). \emph{The theory and practice of item response theory}. Guilford Press.
#'
#' @keywords data
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(deAyala)
#' head(dat)
#' itemstats(dat)
#'
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
#' @examples
#'
#' \dontrun{
#' dat <- expand.table(Bock1997)
#' head(dat)
#' itemstats(dat, use_ts=FALSE)
#'
#' mod <- mirt(dat, 1, 'nominal')
#'
#' # reproduce table 3 in Bock (1997)
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

#' Description of ASVAB data
#'
#' Table of counts extracted from Mislvey (1985). Data the 16 possible
#' response patterns observed for four items from the arithmetic reasoning
#' test of the Armed Services Vocational Aptitude Battery (ASVAB), Form 8A,
#' from samples of white males and females and black males and females.
#'
#' @name ASVAB
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#'
#' Mislevy, R. J. (1985). Estimation of latent group effects.
#' \emph{Journal of the American Statistical Association, 80}, 993-997.
#'
#' @keywords data
#' @examples
#'
#' data(ASVAB)
#' datWM <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, White_Male)))
#' datWF <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, White_Female)))
#' datBM <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, Black_Male)))
#' datBF <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, Black_Female)))
#'
#' dat <- rbind(datWM, datWF, datBM, datBF)
#' sex <- rep(c("Male", "Female", "Male", "Female"),
#'   times=c(nrow(datWM), nrow(datWF), nrow(datBM), nrow(datBF))) |> factor()
#' color <- rep(c("White", "Black"),
#'   times=c(nrow(datWM) + nrow(datWF), nrow(datBM) + nrow(datBF))) |> factor()
#' group <- sex:color
#'
#' itemstats(dat, group=group)
#'
NULL

#' Description of Attitude data
#'
#' Table of counts extracted from Andrich (1988). Data the
#' response patterns observed for an eight item survey.
#'
#' The items in this survey were:
#' \enumerate{
#'   \item Capital punishment is one of the most hideous practices of our time.
#'   \item The state cannot teach the sacredness of human life by destroying it.
#'   \item Capital punishment is not an effective deterrent to crime.
#'   \item I don't believe in capital punishment but I am not sure it
#'     isn't necessary.
#'   \item I think capital punishment is necessary but I wish it were not.
#'   \item Until we find a more civilized way to prevent crime we must have capital
#'    punishment.
#'   \item Capital punishment is justified because it does act as a
#'    deterrent to crime.
#'   \item Capital punishment gives the criminal what he deserves.
#' }
#'
#'
#' @name Attitude
#' @docType data
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#'
#' Andrich, D. (1988). The Application of an Unfolding Model of the PIRT
#' Type to the Measurement of Attitude. \emph{Applied Psychological Measurement, 12}, 33-51.
#'
#' @keywords data
#' @examples
#'
#' head(Attitude)
#' df <- expand.table(Attitude)
#' itemstats(df)
#'
#' # estimate SSLM with estimated " latitude of acceptance" (rho)
#' mod.rho <- mirt(df, 1, itemtype = 'sslm')
#' coef(mod.rho)
#' coef(mod.rho, simplify=TRUE)  # slope-intercept-log_rho
#' coef(mod.rho, simplify=TRUE, IRTpars=TRUE)  # discrimination-difficulty-rho
#' plot(mod.rho)
#' plot(mod.rho, type = 'trace')
#'
#' # without estimating rho, and fixing to rho^2 = 1  (hence, log_rho = -exp(1) =
#'   -2.718282 in order to obtain (exp(exp(log_rho))) = 1)
#' syntax <- "Theta = 1-8
#'            FIXED = (1-8, log_rho1)
#'            START = (1-8, log_rho1, -2.71828)"
#' mod <- mirt(df, syntax, itemtype = 'sslm')  # model found in Andrich (1988)
#' coef(mod)
#' coef(mod, simplify=TRUE)  # slope-intercept-log_rho
#' coef(mod, simplify=TRUE, IRTpars=TRUE)  # discrimination-difficulty-rho
#' plot(mod)
#' plot(mod, type = 'trace') # notice that all curves have a fixed height of .5
#'
#' # goodness of fit (less constrained model fits better)
#' anova(mod, mod.rho) # original model fits much worse
#' M2(mod)
#' M2(mod.rho)
#' itemfit(mod, p.adjust='fdr')
#' itemfit(mod.rho, p.adjust='fdr')
#'
#'
NULL
