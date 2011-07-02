\name{confmirt}
\alias{confmirt}
\alias{coef,confmirt-method}
\alias{summary,confmirt-method}
\alias{residuals,confmirt-method}
\alias{anova,confmirt-method}

\title{Confirmatory Full-Information Item Factor Analysis for Mixed Data Formats}
\description{
  \code{confmirt} fits a conditional (i.e., confirmatory) maximum likelihood factor analysis model to dichotomous and polychotomous data under the item response theory paradigm using Cai's (2010) Metropolis-Hastings Robbins-Monro algorithm.
}

\usage{
confmirt(data, sem.mod, guess = 0, gmeans = 0, ncycles = 2000,
	burnin = 100, SEM.cycles = 50, kdraws = 1, tol = .001, printcycles = TRUE,
    calcLL = TRUE, draws = 2000, debug = FALSE, ...)

\S4method{coef}{confmirt}(object, SE = TRUE, print.gmeans = FALSE, digits = 3, ...)

\S4method{summary}{confmirt}(object, digits = 3, ...)

\S4method{residuals}{confmirt}(object, digits = 3, ...)

\S4method{anova}{confmirt}(object, object2, ...)

}

\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that consists of numerically ordered data.}
  \item{sem.mod}{an object resembling the \code{sem} declaration of a confirmatory factor analysis, but with omitted residual variances.}
  \item{guess}{starting values for the pseudo-guessing parameter. Can be entered as a single value to assign a global guessing parameter or may be entered as a numeric vector for each item. When the value is greater than zero the parameter will be estimated.}
  \item{gmeans}{starting values for the group means. Can be entered as a single value to assign a global parameter or may be entered as a numeric vector for each factor. When the value is not equal to zero the parameter will be estimated.}
  \item{ncycles}{the maximum number of iterations to be performed}
  \item{burnin}{number of burn-in cycles to perform before beginning the SEM stage}
  \item{SEM.cycles}{number of stochastic EM cycles to perform before beginning the MH-RM algorithm}
  \item{kdraws}{number of Metropolis-Hastings imputations of the factor scores at each iteration. Default is 1.}
  \item{tol}{tolerance that must be reached before the model is returned}
  \item{printcycles}{logical; display iteration history during estimation?}
  \item{calcLL}{logical; calculate the log-likelihood via Monte Carlo integration?}
  \item{draws}{the number of Monte Carlo draws to estimate the log-likelihood}
  \item{debug}{logical; turn on debugging features?}
  \item{object}{an object of class \code{confmirt}.}
  \item{object2}{an object of class \code{confmirt}.}
  \item{SE}{logical; print standard errors?}
  \item{print.gmeans}{logical; print factor means?}
  \item{digits}{the number of significant digits to be rounded}
  \item{...}{additional arguments to be passed}
}

\details{

  \code{confmirt} follows the confirmatory item factor analysis strategy by a stochastic version of maximum likelihood estimation described by Cai (2010). The general equation used for multidimensional item response theory in this function is in the logistic form with a scaling correction of 1.702. This correction is applied to allow comparison to mainstream programs such as TESTFACT (2003) and POLYFACT. Missing data are treated as 'missing at random' so that each response vector is included in the estimation. Residuals are computed using the LD statistic (Chen & Thissen, 1997) in the lower diagonal of the matrix returned by \code{residuals}, and Cramer's V above the diagonal. For computing the log-likelihood more accurately see \code{\link{logLik}}.

  Specification of the confirmatory item factor analysis model follows many of the rules in the SEM framework for confirmatory factor analysis. The variances of the latent factors should be fixed to a constant for identification (commonly fixed to 1), and this is the recommended practice since the code is optimized for this method. Parameters may be fixed to specific values by placing a constant in the 4th column of the \code{specify.model} function call with the parameter named \code{NA} in the 3rd column. Item slopes may also be constrained to be equal by labeling them with the same parameter name in the 3rd column. Guessing parameters may be specified for dichotomous items, but are treated as constants, and if a guessing parameter is declared for a polychotomous item it is ignored.

  \code{coef} displays the item parameters with their associated standard errors, while use of \code{summary} transforms the slopes into a factor loadings metric. Also, nested models may be compared by using the \code{anova} function, where a Chi-squared difference test and AIC difference values are displayed.

}

\value{
  \code{confmirt} returns an object of class \code{confmirt}, with the following elements:
  \item{pars}{estimated parameters of the model, the rightmost columns being the intercepts}
  \item{guess}{a vector of the estimated pseudo-guessing parameters}
  \item{gpars}{means and covariance matrix for the group parameters}
  \item{SEpars}{a matrix of standard errors for the parameters estimated}
  \item{SEg}{a vector of standard errors for the pseudo-guessing parameters}
  \item{SEgpars}{a list of standard errors for the group parameters}
  \item{estpars}{a list of logical sets used to identify which parameters were estimated}
  \item{Theta}{last iteration of the Metropolis-Hastings imputation of the factor scores}
  \item{fulldata}{complete data complete with the original item names and values}
  \item{K}{a vector indicating how many unique values each item contains}
  \item{cycles}{number of MHRM cycles performed; will be less than \code{ncycles} if the tolerance is reached}
  \item{itemloc}{a vector indicating where the item location begin}
  \item{converge}{logical; did the model converge?}
  \item{Call}{the function call}
}

\references{

Cai, L. (2010). Metropolis-Hastings Robbins-Monro Algorithm for confirmatory item factor analysis. \emph{Journal of Educational and Behavioral Statistics, 35}, 307-335.

Wood, R., Wilson, D. T., Gibbons, R. D., Schilling, S. G., Muraki, E., & Bock, R. D. (2003). TESTFACT 4 for Windows: Test Scoring, Item Statistics, and Full-information Item Factor Analysis [Computer software]. Lincolnwood, IL: Scientific Software International.

}

\author{Phil Chalmers \email{rphilip.chalmers@gmail.com}}

\seealso{\code{\link{expand.table}},\code{\link{key2binary}},\code{\link{simdata}}}

\examples{
\dontrun{
#simulate data
a <- matrix(c(
1.5,NA,
0.5,NA,
1.0,NA,
1.0,0.5,
 NA,1.5,
 NA,0.5,
 NA,1.0,
 NA,1.0),ncol=2,byrow=TRUE)

d <- matrix(c(
-1.0,NA,NA,
-1.5,NA,NA,
 1.5,NA,NA,
 0.0,NA,NA,
3.0,2.0,-0.5,
2.5,1.0,-1,
2.0,0.0,NA,
1.0,NA,NA),ncol=3,byrow=TRUE)

sigma <- diag(2)
sigma[1,2] <- sigma[2,1] <- .4
dataset <- simdata(a,d,2000,sigma)

#analyses
#CIFA for 2 factor crossed structure

model.1 <- specify.model()
    F1 -> Item_1,   lam11, NA
    F1 -> Item_2,   lam12, NA
    F1 -> Item_3,   lam13, NA
    F1 -> Item_4,   lam14, NA
    F2 -> Item_4,   lam24, NA
    F2 -> Item_5,   lam25, NA
    F2 -> Item_6,   lam26, NA
    F2 -> Item_7,   lam27, NA
    F2 -> Item_8,   lam28, NA
    F1 <-> F2,      phi12, NA
    F1 <-> F1,      NA,     1
    F2 <-> F2,      NA,     1

mod1 <- confmirt(dataset,model.1)
coef(mod1)
summary(mod1)
residuals(mod1)

#fix first slope at 1.5, and set slopes 7 & 8 to be equal
model.2 <- specify.model()
    F1 -> Item_1,      NA, 1.5
    F1 -> Item_2,   lam12, NA
    F1 -> Item_3,   lam13, NA
    F1 -> Item_4,   lam14, NA
    F2 -> Item_4,   lam24, NA
    F2 -> Item_5,   lam25, NA
    F2 -> Item_6,   lam26, NA
    F2 -> Item_7,   lam27, NA
    F2 -> Item_8,   lam27, NA
    F1 <-> F2,      phi12, NA
    F1 <-> F1,      NA,     1
    F2 <-> F2,      NA,     1

mod2 <- confmirt(dataset,model.2)
anova(mod2,mod1)

#####
#bifactor 	
model.3 <- specify.model()
    F1 -> Item_1,   lam11, NA
    F1 -> Item_2,   lam12, NA
    F1 -> Item_3,   lam13, NA
    F1 -> Item_4,   lam14, NA
    F1 -> Item_5,   lam15, NA
    F1 -> Item_6,   lam16, NA
    F1 -> Item_7,   lam17, NA
    F1 -> Item_8,   lam18, NA
    F2 -> Item_1,   lam21, NA
    F2 -> Item_2,   lam22, NA
    F2 -> Item_3,   lam23, NA
    F2 -> Item_4,   lam24, NA
    F3 -> Item_5,   lam35, NA
    F3 -> Item_6,   lam36, NA                	
    F3 -> Item_7,   lam37, NA
    F3 -> Item_8,   lam38, NA
    F1 <-> F1,      NA,     1
    F2 <-> F2,      NA,     1
    F3 <-> F3,      NA,     1
	
mod3 <- confmirt(dataset,model.3)
coef(mod3)
summary(mod3)
residuals(mod3)
anova(mod1,mod3)

     }
}

\keyword{models}
