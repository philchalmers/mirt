\name{multipleGroup}
\alias{anova,MultipleGroupClass-method}
\alias{coef,MultipleGroupClass-method}
\alias{fitted,MultipleGroupClass-method}
\alias{multipleGroup}
\alias{plot,MultipleGroupClass-method}
\alias{residuals,MultipleGroupClass-method}
\alias{summary,MultipleGroupClass-method}
\title{Multiple Group Estimation}
\usage{
  multipleGroup(data, model, group, itemtype = NULL, guess
    = 0, upper = 1, SE = FALSE, SE.type = 'SEM', invariance
    = '', pars = NULL, method = 'EM', constrain = NULL,
    parprior = NULL, calcNull = TRUE, draws = 5000, quadpts
    = NULL, grsm.block = NULL, rsm.block = NULL, key =
    NULL, technical = list(), accelerate = TRUE, verbose =
    TRUE, ...)

  \S4method{coef}{MultipleGroupClass}(object, CI = .95,
    digits = 3, verbose = TRUE, ...)

  \S4method{summary}{MultipleGroupClass}(object, digits =
    3, verbose = TRUE, ...)

  \S4method{anova}{MultipleGroupClass}(object, object2)

  \S4method{residuals}{MultipleGroupClass}(object, ...)

  \S4method{fitted}{MultipleGroupClass}(object, ...)

  \S4method{plot}{MultipleGroupClass}(x, y, type = 'info',
    npts = 50, theta_angle = 45, rot = list(xaxis = -70,
    yaxis = 30, zaxis = 10), ...)
}
\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that
  consists of numerically ordered data, with missing data
  coded as \code{NA}}

  \item{model}{a signel model object returned from
  \code{mirt.model()} declaring how the factor model is to
  be estimated. See \code{\link{mirt.model}} for more
  details}

  \item{group}{a character vector indicating group
  membership}

  \item{invariance}{a character vector containing the
  following possible options: \describe{
  \item{\code{'free_means'}}{for freely estimating all
  latent means (reference group constrained to 0)}
  \item{\code{'free_var'}}{for freely estimating all latent
  variances (reference group constrained to 1's)}
  \item{\code{'free_cov'}}{for freely estimating all latent
  covariances (reference group constrained to an Identity
  matrix)} \item{\code{'free_varcov'}}{calls both
  \code{'free_var'} and \code{'free_cov'}}
  \item{\code{'slopes'}}{to constrain all the slopes to be
  equal across all groups} \item{\code{'intercepts'}}{to
  constrain all the intercepts to be equal across all
  groups, note for nominal models this also includes the
  category specific slope parameters}}

  Additionally, specifying specific item name bundles (from
  \code{colnames(data)}) will constrain all freely
  estimated parameters in each item to be equal across
  groups. This is useful for selecting 'anchor' items for
  vertical and horizontal scaling, and for detecting
  differential item functioning (DIF) across groups}

  \item{guess}{initial (or fixed) values for the
  pseudo-guessing parameter. Can be entered as a single
  value to assign a global guessing parameter or may be
  entered as a numeric vector for each item}

  \item{upper}{initial (or fixed) upper bound parameters
  for 4-PL model. Can be entered as a single value to
  assign a global upper bound parameter or may be entered
  as a numeric vector corresponding to each item}

  \item{accelerate}{see \code{\link{mirt}} for more
  details}

  \item{SE}{logical; estimate the information matrix for
  standard errors?}

  \item{SE.type}{see \code{\link{mirt}} for more details}

  \item{verbose}{logical; display iteration history during
  estimation?}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood}

  \item{quadpts}{the number of quadratures to be used per
  dimensions when \code{method = 'EM'}}

  \item{calcNull}{logical; calculate the Null model for fit
  statics (e.g., TLI)?}

  \item{method}{a character indicating whether to use the
  EM (\code{'EM'}) or the MH-RM (\code{'MHRM'}) algorithm}

  \item{type}{type of plot to view; can be \code{'info'} to
  show the test information function, \code{'infocontour'}
  for the test information contours, \code{'SE'} for the
  test standard error function, \code{'RE'} for the
  relative efficiency plot, and \code{'score'} for the
  expected total score plot}

  \item{theta_angle}{numeric values ranging from 0 to 90
  used in \code{plot}}

  \item{npts}{number of quadrature points to be used for
  plotting features. Larger values make plots look
  smoother}

  \item{rot}{allows rotation of the 3D graphics}

  \item{x}{an object of class \code{mirt} to be plotted or
  printed}

  \item{y}{an unused variable to be ignored}

  \item{key}{see \code{\link{mirt}} for details}

  \item{itemtype}{see \code{\link{mirt}} for details}

  \item{CI}{see \code{\link{mirt}} for details}

  \item{constrain}{see \code{\link{mirt}} for details}

  \item{grsm.block}{see \code{\link{mirt}} for details}

  \item{rsm.block}{see \code{\link{mirt}} for details}

  \item{parprior}{see \code{\link{mirt}} for details}

  \item{pars}{see \code{\link{mirt}} for details}

  \item{object}{an object of class \code{confmirtClass}}

  \item{object2}{an object of class \code{confmirtClass}}

  \item{digits}{the number of significant digits to be
  rounded}

  \item{...}{additional arguments to be passed}

  \item{technical}{list specifying subtle parameters that
  can be adjusted. These values are \describe{
  \item{NCYCLES}{max number of cycles; default 2000 for
  MHRM and 500 for EM} \item{MAXQUAD}{maximum number of
  quadratures; default 10000} \item{BURNIN}{number of burn
  in cycles (stage 1); default 150} \item{SEMCYCLES}{number
  of SEM cycles (stage 2); default 50} \item{TOL}{minimum
  threshold tolerance for convergence. If MH-RM, must occur
  on three consecutive occations; default is .001 for MH-RM
  and .0001 for EM} \item{set.seed}{seed number used during
  estimation. Default is 12345} \item{gain}{a vector of two
  values specifying the numerator and exponent values for
  the RM gain function \eqn{(val1 / cycle)^val2}. Default
  is \code{c(0.15,0.65)}} \item{return_newconstrain}{if
  \code{TRUE} returns a list consisting of the constraints
  to be used just before estimation begins} }}
}
\description{
  \code{multipleGroup} performs a full-information
  maximum-likelihood multiple group analysis for
  dichotomous and polytomous data under the item response
  theory paradigm using either Cai's (2010)
  Metropolis-Hastings Robbins-Monro (MHRM) algorithm or
  with an EM algorithm approach.
}
\details{
  By default the estimation in \code{multipleGroup} assumes
  that the models are maximally independent, and therefore
  could initially be performed by sub-setting the data and
  running identical models with \code{mirt} and aggregating
  the results (e.g., log-likelihood). However, constrains
  may be imposed across groups by invoking various
  \code{invariance} keywords and \code{constrain = ...}
  arguments, by inputting user specified design matrix from
  \code{mod2values} or from passing \code{pars = 'values'},
  or by supplying a \code{constrain} list for user defined
  equality constraints between parameters.
}
\examples{
\dontrun{
#single factor
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
itemtype <- rep('dich', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))
models <- mirt.model('F1 = 1-15')

mod_configural <- multipleGroup(dat, models, group = group) #completely separate analyses

mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes')) #equal slopes
#equal intercepts, free variance and means
mod_scalar2 <- multipleGroup(dat, models, group = group,
                             invariance=c('slopes', 'intercepts', 'free_var','free_means'))
mod_scalar1 <- multipleGroup(dat, models, group = group,  #fixed means
                             invariance=c('slopes', 'intercepts', 'free_var'))
mod_fullconstrain <- multipleGroup(dat, models, group = group,
                             invariance=c('slopes', 'intercepts'))

summary(mod_scalar2)
coef(mod_scalar2)
residuals(mod_scalar2)
plot(mod_configural)
plot(mod_configural, type = 'score')
itemplot(mod_configural, 2)
itemplot(mod_configural, 2, type = 'RE')

anova(mod_metric, mod_configural) #equal slopes only
anova(mod_scalar2, mod_metric) #equal intercepts, free variance and mean
anova(mod_scalar1, mod_scalar2) #fix mean
anova(mod_fullconstrain, mod_scalar1) #fix variance


#test whether first 6 slopes should be equal across groups
values <- multipleGroup(dat, models, group = group, pars = 'values')
values
constrain <- list(c(1, 63), c(5,67), c(9,71), c(13,75), c(17,79), c(21,83))
equalslopes <- multipleGroup(dat, models, group = group, constrain = constrain)
anova(equalslopes, mod_configural)

#same as above, but using mirt.model syntax
newmodel <- mirt.model('
    F = 1-15
    CONSTRAINB = (1-6, a1)')
equalslopes <- multipleGroup(dat, newmodel, group = group)
coef(equalslopes)

#############
#DIF test for each item (using all other items as anchors)
itemnames <- colnames(dat)
refmodel <- multipleGroup(dat, models, group = group,
                             invariance=c('free_means', 'free_varcov', itemnames))

#loop over items (in practice, run in parallel to increase speed)
estmodels <- vector('list', ncol(dat))
for(i in 1:ncol(dat))
    estmodels[[i]] <- multipleGroup(dat, models, group = group, verbose = FALSE, calcNull=FALSE,
                             invariance=c('free_means', 'free_varcov', itemnames[-i]))

(anovas <- lapply(estmodels, anova, object2=refmodel))

#family-wise error control
p <- do.call(rbind, lapply(anovas, function(x) x[2,9]))
p.adjust(p, method = 'BH')

#same as above, except only test if slopes vary (1 df)
#constrain all intercepts
estmodels <- vector('list', ncol(dat))
for(i in 1:ncol(dat))
    estmodels[[i]] <- multipleGroup(dat, models, group = group, verbose = FALSE, calcNull=FALSE,
                             invariance=c('free_means', 'free_varcov', 'intercepts',
                             itemnames[-i]))


(anovas <- lapply(estmodels, anova, object2=refmodel))

#############
#multiple factors

a <- matrix(c(abs(rnorm(5,1,.3)), rep(0,15),abs(rnorm(5,1,.3)),
rep(0,15),abs(rnorm(5,1,.3))), 15, 3)
d <- matrix(rnorm(15,0,.7),ncol=1)
mu <- c(-.4, -.7, .1)
sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
itemtype <- rep('dich', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = mu, sigma = sigma)
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

#group models
model <- mirt.model('
   F1 = 1-5
   F2 = 6-10
   F3 = 11-15')

#EM approach (not as accurate with 3 factors, but generally good for quick model comparisons)
mod_configural <- multipleGroup(dat, model, group = group) #completely separate analyses
mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes')) #equal slopes
mod_scalar <- multipleGroup(dat, model, group = group, #equal means, slopes, intercepts
                             invariance=c('slopes', 'intercepts', 'free_varcov'))
mod_fullconstrain <- multipleGroup(dat, model, group = group, #equal means, slopes, intercepts
                             invariance=c('slopes', 'intercepts'))

anova(mod_metric, mod_configural)
anova(mod_scalar, mod_metric)
anova(mod_fullconstrain, mod_scalar)

#same as above, but with MHRM (more accurate with 3 factors, but slower)
#define mirt cluster to compute log-likelihoods faster
mirtCluster()

#completely separate analyses
mod_configural <- multipleGroup(dat, models, group = group, method = 'MHRM')
#equal slopes
mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'MHRM')
#equal means, slopes, intercepts
mod_scalar <- multipleGroup(dat, models, group = group, method = 'MHRM',
                             invariance=c('slopes', 'intercepts', 'free_varcov'))
#equal means, slopes, intercepts
mod_fullconstrain <- multipleGroup(dat, models, group = group, method = 'MHRM',
                             invariance=c('slopes', 'intercepts'))

anova(mod_metric, mod_configural)
anova(mod_scalar, mod_metric)
anova(mod_fullconstrain, mod_scalar)

############
#polytomous item example
set.seed(12345)
a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d <- matrix(rnorm(15,0,.7),ncol=1)
d <- cbind(d, d-1, d-2)
itemtype <- rep('graded', nrow(a))
N <- 1000
dataset1 <- simdata(a, d, N, itemtype)
dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))
models <- mirt.model('F1 = 1-15')
models2 <- mirt.model('
   F1 = 1-10
   F2 = 10-15')

mod_configural <- multipleGroup(dat, models, group = group)
plot(mod_configural)
plot(mod_configural, type = 'SE')
itemplot(mod_configural, 1)
itemplot(mod_configural, 1, type = 'info')
fs <- fscores(mod_configural)
head(fs[["D1"]])
fscores(mod_configural, method = 'EAPsum')

mod_configural2 <- multipleGroup(dat, models2, group = group)
plot(mod_configural2, type = 'infocontour')
plot(mod_configural2, type = 'SE')
plot(mod_configural2, type = 'RE')
itemplot(mod_configural2, 10)

}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{expand.table}}, \code{\link{key2binary}},
  \code{\link{mirt.model}}, \code{\link{mirt}},
  \code{\link{bfactor}}, \code{\link{multipleGroup}},
  \code{\link{mixedmirt}}, \code{\link{wald}},
  \code{\link{itemplot}}, \code{\link{fscores}},
  \code{\link{fitIndices}}, \code{\link{extract.item}},
  \code{\link{iteminfo}}, \code{\link{testinfo}},
  \code{\link{probtrace}}, \code{\link{boot.mirt}},
  \code{\link{imputeMissing}}, \code{\link{itemfit}},
  \code{\link{mod2values}}, \code{\link{simdata}},
  \code{\link{createItem}}, \code{\link{mirtCluster}}
}
\keyword{models}

