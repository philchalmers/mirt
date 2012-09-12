\name{multipleGroup}
\alias{anova,MultipleGroupClass-method}
\alias{coef,MultipleGroupClass-method}
\alias{multipleGroup}
\title{Multiple Group Estimation}
\usage{
  multipleGroup(data, model, group, itemtype = NULL, guess
    = 0, upper = 1, free.start = NULL, invariance = '',
    method = 'MHRM', constrain = NULL, startvalues = NULL,
    parprior = NULL, freepars = NULL, draws = 2000, quadpts
    = NULL, technical = list(), debug = FALSE, verbose =
    TRUE)

  \S4method{coef}{MultipleGroupClass}(object, digits = 3,
    verbose = TRUE, ...)

  \S4method{anova}{MultipleGroupClass}(object, object2)
}
\arguments{
  \item{data}{a \code{matrix} or \code{data.frame} that
  consists of numerically ordered data, with missing data
  coded as \code{NA}}

  \item{model}{an object or named list of objects returned
  from \code{confmirt.model()} declaring how the factor
  model is to be estimated. The names of the list input
  must correspond to the unique values in the \code{group}
  variable. See \code{\link{confmirt.model}} for more
  details}

  \item{group}{a character vector indicating group
  membership}

  \item{invariance}{a character vector containing the
  following possible options: \describe{
  \item{\code{'free_means'}}{for freely estimating all
  latent means (reference group constrained to 0)}
  \item{\code{'free_varcov'}}{for freely estimating the
  variance-covariance matrix accross groups (reference
  group has variances equal to 1, but freely estimated
  covariance terms if specified in the model)}
  \item{\code{'covariances'}}{to constrain all the
  covariance parameters to be equal, note that this only
  makes sense if the factor variances are the same (i.e.,
  unity)} \item{\code{'slopes'}}{to constrain all the
  slopes to be equal across all groups}
  \item{\code{'intercepts'}}{to constrain all the
  intercepts to be equal across all groups, note for
  nominal models this also includes the category specific
  slope parameters} }}

  \item{guess}{initial (or fixed) values for the
  pseudo-guessing parameter. Can be entered as a single
  value to assign a global guessing parameter or may be
  entered as a numeric vector for each item}

  \item{upper}{initial (or fixed) upper bound parameters
  for 4-PL model. Can be entered as a single value to
  assign a global upper bound parameter or may be entered
  as a numeric vector corresponding to each item}

  \item{free.start}{a list containing the start value and
  logical indicating whether a given parameter is to be
  freely estimated. Each element of the list consists of
  three components, the parameter number, the starting (or
  fixed) value, and a logical to indicate whether the
  parameter is free. For example, \code{free.start =
  list(c(20,0.0,FALSE), c(10,1.5,TRUE))} would fix
  parameter 20 to 0.0 while parameter 10 would be freely
  estimated from a starting value of 1.5. Note that this
  will override the values specified by a user defined
  \code{startvalues} or \code{freepars} input for the
  specified parameters, and this may conflict with the
  \code{invariance} input (e.g., freeing slopes manually
  while specifying \code{invariance = 'slopes'} is
  ambiguous and should be avoided).}

  \item{verbose}{logical; display iteration history during
  estimation?}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood}

  \item{quadpts}{the number of quadratures to be used per
  dimensions when \code{method = 'EM'}}

  \item{method}{a character indicating whether to use the
  EM (\code{'EM'}) or the MH-RM (\code{'MHRM'}) algorithm}

  \item{itemtype}{type of items to be modeled, declared as
  a vector for each item or a single value which will be
  repeated globally. The NULL default assumes that the
  items are ordinal or 2PL, however they may be changed to
  the following: 'Rasch', '1PL', '2PL', '3PL', '3PLu',
  '4PL', 'graded', 'gpcm', 'nominal', 'mcm', and
  'partcomp', for the Rasch/partial credit, 1 and 2
  parameter logistic, 3 parameter logistic (lower asymptote
  and upper), 4 parameter logistic, graded response model,
  generalized partial credit model, nominal model, multiple
  choice model, and partially compensatory model,
  respectively. The default assumes that items follow a
  '2PL' or 'graded' format If \code{NULL} the default
  assumes that the data follow a '2PL' or 'graded' format}

  \item{constrain}{a list of user declared equality
  constraints. To see how to define the parameters
  correctly use \code{constrain = 'index'} initially to see
  how the parameters are labeled. To constrain parameters
  to be equal create a list with separate concatenated
  vectors signifying which parameters to constrain. For
  example, to set parameters 1 and 5 equal, and also set
  parameters 2, 6, and 10 equal use \code{constrain =
  list(c(1,5), c(2,6,10))}}

  \item{parprior}{a list of user declared prior item
  probabilities. To see how to define the parameters
  correctly use \code{parprior = 'index'} initially to see
  how the parameters are labeled. Can define either normal
  (normally for slopes and intercepts) or beta (for
  guessing and upper bounds) prior probabilities. Note that
  for upper bounds the value used in the prior is 1 - u so
  that the lower and upper bounds can function the same. To
  specify a prior the form is c('priortype', ...), where
  normal priors are \code{parprior = list(c(parnumber,
  'norm', mean, sd))} and betas are \code{parprior =
  list(c(parnumber, 'beta', alpha, beta))}.}

  \item{freepars}{a list of user declared logical values
  indicating which parameters to estimate. To see how to
  define the parameters correctly use \code{freepars =
  'index'} initially to see how the parameters are labeled.
  These values may be modified and input back into the
  function by using \code{freepars=newfreepars}. Note that
  user input values must match what the default structure
  would have been}

  \item{startvalues}{a list of user declared start values
  for parameters. To see how to define the parameters
  correctly use \code{startvalues = 'index'} initially to
  see what the defaults would noramlly be. These values may
  be modified and input back into the function by using
  \code{startavlues=newstartvalues}. Note that user input
  values must match what the default structure would have
  been}

  \item{debug}{logical; turn on debugging features?}

  \item{object}{an object of class \code{confmirtClass}}

  \item{object2}{an object of class \code{confmirtClass}}

  \item{digits}{the number of significant digits to be
  rounded}

  \item{technical}{list specifying subtle parameters that
  can be adjusted. These values are \describe{
  \item{NCYCLES}{max number of cycles; default 2000 for
  MHRM and 300 for EM} \item{MAXQUAD}{maximum number of
  quadratures; default 10000} \item{MSTEPMAXIT}{number of
  M-step iterations; default 15} \item{BURNIN}{number of
  burn in cycles (stage 1); default 150}
  \item{SEMCYCLES}{number of SEM cycles (stage 2); default
  50} \item{KDRAWS}{number of parallel MH sets to be drawn;
  default 1} \item{TOL}{minimum threshold tolerance for
  convergence of MH-RM, must occur on three consecutive
  occations; default .001} \item{set.seed}{seed number used
  during estimation. Default is 12345} \item{gain}{a vector
  of three values specifying the numerator, exponent, and
  subtracted values for the RM gain value. Default is
  \code{c(0.05,0.5,0.004)}} \item{return_newconstrain}{if
  \code{TRUE} returns a list consisting of the constraints
  to be used just before estimation begins} }}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{multipleGroup} performes a full-information
  maximum-likelihood multiple group analysis for
  dichotomous and polytomous data under the item response
  theory paradigm using either Cai's (2010)
  Metropolis-Hastings Robbins-Monro algorithm or with an EM
  approach.
}
\details{
  By default the estimation in \code{multipleGroup} assumes
  that the models are maximimally independent, and
  therefore could initially be performed by subsetting the
  data and running identical models with \code{confmirt} or
  \code{mirt} and aggregating the results (e.g.,
  log-likelihood). However, constrains may be imposed
  accross groups by invoking various \code{invariance}
  keywords or by inputing user defined \code{freepars},
  \code{constrain}, and \code{startvalues} lists.
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
models <- confmirt.model()
   F1 = 1-15


mod_configural <- multipleGroup(dat, models, group = group, method = 'EM') #completely seperate analyses
mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes'), method = 'EM') #equal slopes
mod_scalar2 <- multipleGroup(dat, models, group = group, method = 'EM',  #equal intercepts, free variance and means
                             invariance=c('slopes', 'intercepts', 'free_varcov','free_means'))
mod_scalar1 <- multipleGroup(dat, models, group = group, method = 'EM', #fixed means
                             invariance=c('slopes', 'intercepts', 'free_varcov'))
mod_fullconstrain <- mirt(data, models) #fix variance (equivelent to full constrain)

anova(mod_metric, mod_configural) #equal slopes only
anova(mod_scalar2, mod_metric) #equal intercepts, free variance and mean
anova(mod_scalar1, mod_scalar2) #fix mean
anova(mod_fullconstrain, mod_scalar1) #fix variance


#Wald test can be useful here too
#compare whether intercepts should be equal
index <- multipleGroup(dat, models, group = group, constrain = 'index')
index
nitems <- ncol(dat)
L <- matrix(0, nitems, 124)
for(i in 1:nitems){
     L[i, index[[1]][[i]][2]] <- 1
     L[i, index[[2]][[i]][2]] <- -1
 }
wald(L, mod_configural)


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
model1 <- confmirt.model()
   F1 = 1-5
   F2 = 6-10
   F3 = 11-15


model2 <- confmirt.model()
   F1 = 1-5
   F2 = 6-10
   F3 = 11-15
   COV = F1*F2, F1*F3, F2*F3


models <- list(D1=model1, D2=model2) #note the names match the groups

mod_configural <- multipleGroup(dat, models, group = group) #completely seperate analyses
mod_metric <- multipleGroup(dat, models, group = group, invariance=c('slopes')) #equal slopes
mod_scalar <- multipleGroup(dat, models, group = group, #equal means, slopes, intercepts
                             invariance=c('slopes', 'intercepts', 'free_varcov'))
mod_fullconstrain <- confmirt(data, models)

anova(mod_metric, mod_configural)
anova(mod_scalar, mod_metric)
anova(mod_fullconstrain, mod_scalar)
}
}
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{expand.table}}, \code{\link{key2binary}},
  \code{\link{simdata}}, \code{\link{confmirt.model}}
}
\keyword{models}

