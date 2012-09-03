\name{multipleGroup}
\alias{anova,multipleGroup-method}
\alias{coef,multipleGroup-method}
\alias{multipleGroup}
\title{Multiple Group Estimation}
\usage{
  multipleGroup(data, model, group, itemtype = NULL, guess
    = 0, upper = 1, invariance = '', constrain = NULL,
    startvalues = NULL, parprior = NULL, freepars = NULL,
    draws = 2000, technical = list(), debug = FALSE,
    verbose = TRUE)

  \S4method{coef}{multipleGroup}(object, digits = 3,
    verbose = TRUE, ...)

  \S4method{anova}{multipleGroup}(object, object2, ...)
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

  \item{verbose}{logical; display iteration history during
  estimation?}

  \item{draws}{the number of Monte Carlo draws to estimate
  the log-likelihood}

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
  \item{NCYCLES}{max number of MH-RM cycles; default 2000}
  \item{BURNIN}{number of burn in cycles (stage 1); default
  150} \item{SEMCYCLES}{number of SEM cycles (stage 2);
  default 50} \item{KDRAWS}{number of parallel MH sets to
  be drawn; default 1} \item{TOL}{minimum threshold
  tolerance for convergence of MH-RM, must occur on three
  consecutive occations; default .001} \item{set.seed}{seed
  number used during estimation. Default is 12345}
  \item{guess.prior.n}{a scalar or vector for the weighting
  of the beta priors for guessing parameters (default is
  50, typical ranges are from 2 to 500). If a scalar is
  specified this is used globally, otherwise a numeric
  vector of size \code{ncol(data)} can be used to
  correspond to particular items (NA values use the
  default)} \item{gain}{a vector of three values specifying
  the numerator, exponent, and subtracted values for the RM
  gain value. Default is \code{c(0.05,0.5,0.004)}}
  \item{return_newconstrain}{if \code{TRUE} returns a list
  consisting of the constraints to be used just before
  estimation begins} }}

  \item{...}{additional arguments to be passed}
}
\description{
  \code{multipleGroup} performes a full-information
  maximum-likelihood multiple group analysis for
  dichotomous and polytomous data under the item response
  theory paradigm using Cai's (2010) Metropolis-Hastings
  Robbins-Monro algorithm.
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
\author{
  Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
  \code{\link{expand.table}}, \code{\link{key2binary}},
  \code{\link{simdata}}, \code{\link{confmirt.model}}
}
\keyword{models}

