#' Multiple Group Estimation
#'
#' \code{multipleGroup} performs a full-information
#' maximum-likelihood multiple group analysis for any combination of dichotomous and polytomous
#' data under the item response theory paradigm using either Cai's (2010)
#' Metropolis-Hastings Robbins-Monro (MHRM) algorithm or with an EM algorithm approach. This
#' function may be used for detecting differential item functioning (DIF), thought the
#' \code{\link{DIF}} function may provide a more convenient approach. If the grouping
#' variable is not specified then the \code{dentype} input can be modified to fit
#' mixture models to estimate any latent group components.
#'
#' By default the estimation in \code{multipleGroup} assumes that the models are maximally
#' independent, and therefore could initially be performed by sub-setting the data and running
#' identical models with \code{mirt} and aggregating the results (e.g., log-likelihood).
#' However, constrains may be automatically imposed across groups by invoking various
#' \code{invariance} keywords. Users may also supply a list of parameter equality constraints
#' to by \code{constrain} argument, of define equality constraints using the
#' \code{\link{mirt.model}} syntax (recommended).
#'
#' @return function returns an object of class \code{MultipleGroupClass}
#'   (\link{MultipleGroupClass-class}).
#'
#' @aliases multipleGroup
#' @param data a \code{matrix} or \code{data.frame} that consists of
#'   numerically ordered data, organized in the form of integers,
#'    with missing data coded as \code{NA}
#' @param model string to be passed to, or a model object returned from, \code{\link{mirt.model}}
#'   declaring how the global model is to be estimated (useful to apply constraints here)
#' @param group a \code{character} or \code{factor} vector indicating group membership. If a \code{character}
#'   vector is supplied this will be automatically transformed into a \code{\link{factor}} variable.
#'   As well, the first level of the (factorized) grouping variable will be treated as the "reference" group
#' @param itemtype can be same type of input as is documented in \code{\link{mirt}}, however may also be a
#'   \code{ngroups} by \code{nitems} matrix specifying the type of IRT models for each group, respectively.
#'   Rows of this input correspond to the levels of the \code{group} input. For mixture models the rows correspond
#'   to the respective mixture grouping variables to be constructed, and the IRT models should be within these
#'   mixtures
#' @param invariance a character vector containing the following possible options:
#'   \describe{
#'     \item{\code{'free_mean'} or \code{'free_means'}}{freely estimate all latent means in all focal groups
#'       (reference group constrained to a vector of 0's)}
#'     \item{\code{'free_var'}, \code{'free_vars'}, \code{'free_variance'}, or \code{'free_variances'}}{
#'       freely estimate all latent variances in focal groups
#'       (reference group variances all constrained to 1)}
#'     \item{\code{'slopes'}}{to constrain all the slopes to be equal across all groups}
#'     \item{\code{'intercepts'}}{to constrain all the intercepts to be equal across all
#'       groups, note for nominal models this also includes the category specific slope parameters}
#'    }
#'
#'   Additionally, specifying specific item name bundles (from \code{colnames(data)}) will
#'   constrain all freely estimated parameters in each item to be equal across groups. This is
#'   useful for selecting 'anchor' items for vertical and horizontal scaling, and for detecting
#'   differential item functioning (DIF) across groups
#' @param method a character object that is either \code{'EM'}, \code{'QMCEM'}, or \code{'MHRM'}
#'   (default is \code{'EM'}). See \code{\link{mirt}} for details
#' @param itemdesign see \code{\link{mirt}} for details
#' @param item.formula see \code{\link{mirt}} for details
#' @param dentype type of density form to use for the latent trait parameters. Current options include
#'   all of the methods described in \code{\link{mirt}}, as well as
#'
#'   \itemize{
#'     \item \code{'mixture-#'} estimates mixtures of Gaussian distributions,
#'       where the \code{#} placeholder represents the number of potential grouping variables
#'       (e.g., \code{'mixture-3'} will estimate 3 underlying classes). Each class is
#'       assigned the group name \code{MIXTURE_#}, where \code{#} is the class number.
#'
#'       Note that internally the mixture coefficients are stored as log values where
#'       the first mixture group coefficient is fixed at 0. Additionally, it is recommended
#'       to use the \code{nruns} argument as mixture IRT models are known to contain
#'       local maximums
#'    }
#'
#' @param nruns a numeric value indicating how many times the model should be fit to the data
#'   when using random starting values, which is particularly useful
#'   when evaluating mixture IRT Models. If greater than 1, \code{GenRandomPars} is set to \code{TRUE}
#'   by default. Using this returns a list of fitted model objects, where the model
#'   with the highest log-likelihood should generally be selected as the model
#'   best associated with the MLE (this is done automatically if \code{return_max = TRUE}).
#'   Note that if a \code{\link{mirtCluster}} was
#'   defined earlier then the runs will be run in parallel
#'
#' @param return_max logical; when \code{nruns > 1}, return the model that has the most optimal
#'   maximum likelihood criteria? If FALSE, returns a list of all the estimated objects
#' @param GenRandomPars see \code{\link{mirt}} for details
#' @param verbose see \code{\link{mirt}} for details
#' @param ... additional arguments to be passed to the estimation engine. See \code{\link{mirt}}
#'   for details and examples
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Magnus, B. E. and  Garnier-Villarreal (2022). A multidimensional zero-inflated
#' graded response model for ordinal symptom data. \emph{Psychological Methods}, 27, 261-279.
#'
#' Wall, M., M., Park, J., Y., and Moustaki I. (2015). IRT modeling in the presence of
#' zero-inflation with application to psychiatric disorder severity.
#' \emph{Applied Psychological Measurement} 39: 583-597.
#'
#' @seealso \code{\link{mirt}}, \code{\link{DIF}}, \code{\link{extract.group}}, \code{\link{DRF}}
#' @keywords models
#' @export multipleGroup
#' @examples
#' \donttest{
#'
#' # single factor
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' # marginal information
#' itemstats(dat)
#'
#' # conditional information
#' itemstats(dat, group=group)
#'
#' mod_configural <- multipleGroup(dat, 1, group = group) #completely separate analyses
#' # limited information fit statistics
#' M2(mod_configural)
#'
#' mod_metric <- multipleGroup(dat, 1, group = group, invariance=c('slopes')) #equal slopes
#' # equal intercepts, free variance and means
#' mod_scalar2 <- multipleGroup(dat, 1, group = group,
#'                              invariance=c('slopes', 'intercepts', 'free_var','free_means'))
#' mod_scalar1 <- multipleGroup(dat, 1, group = group,  #fixed means
#'                              invariance=c('slopes', 'intercepts', 'free_var'))
#' mod_fullconstrain <- multipleGroup(dat, 1, group = group,
#'                              invariance=c('slopes', 'intercepts'))
#' extract.mirt(mod_fullconstrain, 'time') #time of estimation components
#'
#' # optionally use Newton-Raphson for (generally) faster convergence in the
#' #  M-step's, though occasionally less stable
#' mod_fullconstrain <- multipleGroup(dat, 1, group = group, optimizer = 'NR',
#'                              invariance=c('slopes', 'intercepts'))
#' extract.mirt(mod_fullconstrain, 'time') #time of estimation components
#'
#' summary(mod_scalar2)
#' coef(mod_scalar2, simplify=TRUE)
#' residuals(mod_scalar2)
#' plot(mod_configural)
#' plot(mod_configural, type = 'info')
#' plot(mod_configural, type = 'trace')
#' plot(mod_configural, type = 'trace', which.items = 1:4)
#' itemplot(mod_configural, 2)
#' itemplot(mod_configural, 2, type = 'RE')
#'
#' anova(mod_metric, mod_configural) #equal slopes only
#' anova(mod_scalar2, mod_metric) #equal intercepts, free variance and mean
#' anova(mod_scalar1, mod_scalar2) #fix mean
#' anova(mod_fullconstrain, mod_scalar1) #fix variance
#'
#' # compared all at once (in order of most constrained to least)
#' anova(mod_fullconstrain, mod_scalar2, mod_configural)
#'
#'
#' # test whether first 6 slopes should be equal across groups
#' values <- multipleGroup(dat, 1, group = group, pars = 'values')
#' values
#' constrain <- list(c(1, 63), c(5,67), c(9,71), c(13,75), c(17,79), c(21,83))
#' equalslopes <- multipleGroup(dat, 1, group = group, constrain = constrain)
#' anova(equalslopes, mod_configural)
#'
#' # same as above, but using mirt.model syntax
#' newmodel <- '
#'     F = 1-15
#'     CONSTRAINB = (1-6, a1)'
#' equalslopes <- multipleGroup(dat, newmodel, group = group)
#' coef(equalslopes, simplify=TRUE)
#'
#' ############
#' # vertical scaling (i.e., equating when groups answer items others do not)
#' dat2 <- dat
#' dat2[group == 'D1', 1:2] <- dat2[group != 'D1', 14:15] <- NA
#' head(dat2)
#' tail(dat2)
#'
#' # items with missing responses need to be constrained across groups for identification
#' nms <- colnames(dat2)
#' mod <- multipleGroup(dat2, 1, group, invariance = nms[c(1:2, 14:15)])
#'
#' # this will throw an error without proper constraints (SEs cannot be computed either)
#' # mod <- multipleGroup(dat2, 1, group)
#'
#' # model still does not have anchors, therefore need to add a few (here use items 3-5)
#' mod_anchor <- multipleGroup(dat2, 1, group,
#'                             invariance = c(nms[c(1:5, 14:15)], 'free_means', 'free_var'))
#' coef(mod_anchor, simplify=TRUE)
#'
#' # check if identified by computing information matrix
#' mod_anchor <- multipleGroup(dat2, 1, group, pars = mod2values(mod_anchor), TOL=NaN, SE=TRUE,
#'                             invariance = c(nms[c(1:5, 14:15)], 'free_means', 'free_var'))
#' mod_anchor
#' coef(mod_anchor)
#' coef(mod_anchor, printSE=TRUE)
#'
#'
#' #############
#' # DIF test for each item (using all other items as anchors)
#' itemnames <- colnames(dat)
#' refmodel <- multipleGroup(dat, 1, group = group, SE=TRUE,
#'                           invariance=c('free_means', 'free_var', itemnames))
#'
#' # loop over items (in practice, run in parallel to increase speed). May be better to use ?DIF
#' estmodels <- vector('list', ncol(dat))
#' for(i in 1:ncol(dat))
#'     estmodels[[i]] <- multipleGroup(dat, 1, group = group, verbose = FALSE,
#'                              invariance=c('free_means', 'free_var', itemnames[-i]))
#' anova(refmodel, estmodels[[1]])
#' (anovas <- lapply(estmodels, function(x, refmodel) anova(refmodel, x),
#'    refmodel=refmodel))
#'
#' # family-wise error control
#' p <- do.call(rbind, lapply(anovas, function(x) x[2, 'p']))
#' p.adjust(p, method = 'BH')
#'
#' # same as above, except only test if slopes vary (1 df)
#' # constrain all intercepts
#' estmodels <- vector('list', ncol(dat))
#' for(i in 1:ncol(dat))
#'     estmodels[[i]] <- multipleGroup(dat, 1, group = group, verbose = FALSE,
#'                              invariance=c('free_means', 'free_var', 'intercepts',
#'                              itemnames[-i]))
#'
#' (anovas <- lapply(estmodels, function(x, refmodel) anova(refmodel, x),
#'    refmodel=refmodel))
#'
#' # quickly test with Wald test using DIF()
#' mod_configural2 <- multipleGroup(dat, 1, group = group, SE=TRUE)
#' DIF(mod_configural2, which.par = c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr')
#'
#'
#'
#' #############
#' # Three group model where the latent variable parameters are constrained to
#' # be equal in the focal groups
#'
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dataset3 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2, dataset3)
#' group <- rep(c('D1', 'D2', 'D3'), each=N)
#'
#' # marginal information
#' itemstats(dat)
#'
#' # conditional information
#' itemstats(dat, group=group)
#'
#' model <- 'F1 = 1-15
#'           FREE[D2, D3] = (GROUP, MEAN_1), (GROUP, COV_11)
#'           CONSTRAINB[D2,D3] = (GROUP, MEAN_1), (GROUP, COV_11)'
#'
#' mod <- multipleGroup(dat, model, group = group, invariance = colnames(dat))
#' coef(mod, simplify=TRUE)
#'
#' #############
#' # Testing main effects in multiple independent grouping variables
#' set.seed(1234)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 500
#'
#' # generated data have interaction effect for latent means, as well as a
#' # main effect across D but no main effect across G
#' d11 <- simdata(a, d, N, itemtype, mu = 0)
#' d12 <- simdata(a, d, N, itemtype, mu = 0)
#' d13 <- simdata(a, d, N, itemtype, mu = 0)
#' d21 <- simdata(a, d, N, itemtype, mu = 1/2)
#' d22 <- simdata(a, d, N, itemtype, mu = 1/2)
#' d23 <- simdata(a, d, N, itemtype, mu = -1)
#' dat <- do.call(rbind, list(d11, d12, d13, d21, d22, d23))
#' group <- rep(c('G1.D1', 'G1.D2', 'G1.D3', 'G2.D1', 'G2.D2', 'G2.D3'), each=N)
#' table(group)
#'
#' \dontrun{
#' # in practice, group would be organized in a data.frame as follows
#' df <- data.frame(group)
#' dfw <- tidyr::separate_wider_delim(df, group, delim='.', names = c('G', 'D'))
#' head(dfw)
#'
#' # for use with multipleGroup() combine into a single long group
#' group <- with(dfw, factor(G):factor(D))
#'
#' # conditional information
#' itemstats(dat, group=group)
#'
#' mod <- multipleGroup(dat, group = group, SE=TRUE,
#'                      invariance = c(colnames(dat), 'free_mean', 'free_var'))
#' coef(mod, simplify=TRUE)
#' sapply(coef(mod, simplify=TRUE), \(x) unname(x$means)) # mean estimates
#' wald(mod) # view parameter names for later testing
#'
#' # test for main effect over G group (manually compute marginal mean)
#' wald(mod, "0 + MEAN_1.123 + MEAN_1.185 = MEAN_1.247 + MEAN_1.309 + MEAN_1.371")
#'
#' # test for main effect over D group  (manually compute marginal means)
#' wald(mod, c("0 + MEAN_1.247 = MEAN_1.123 + MEAN_1.309",
#'             "0 + MEAN_1.247 = MEAN_1.185 + MEAN_1.371"))
#'
#' # post-hoc tests (better practice would include p.adjust() )
#' wald(mod, "0 + MEAN_1.247 = MEAN_1.123 + MEAN_1.309") # D1 vs D2
#' wald(mod, "0 + MEAN_1.247 = MEAN_1.185 + MEAN_1.371") # D1 vs D3
#' wald(mod, "MEAN_1.123 + MEAN_1.309 = MEAN_1.185 + MEAN_1.371") # D2 vs D3
#' }
#'
#' #############
#' # multiple factors
#'
#' a <- matrix(c(abs(rnorm(5,1,.3)), rep(0,15),abs(rnorm(5,1,.3)),
#'      rep(0,15),abs(rnorm(5,1,.3))), 15, 3)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' mu <- c(-.4, -.7, .1)
#' sigma <- matrix(c(1.21,.297,1.232,.297,.81,.252,1.232,.252,1.96),3,3)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = mu, sigma = sigma)
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#'
#' # group models
#' model <- '
#'    F1 = 1-5
#'    F2 = 6-10
#'    F3 = 11-15'
#'
#' # define mirt cluster to use parallel architecture
#' if(interactive()) mirtCluster()
#'
#' # EM approach (not as accurate with 3 factors, but generally good for quick model comparisons)
#' mod_configural <- multipleGroup(dat, model, group = group) #completely separate analyses
#' mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes')) #equal slopes
#' mod_fullconstrain <- multipleGroup(dat, model, group = group, #equal means, slopes, intercepts
#'                              invariance=c('slopes', 'intercepts'))
#'
#' anova(mod_metric, mod_configural)
#' anova(mod_fullconstrain, mod_metric)
#'
#' # same as above, but with MHRM (generally  more accurate with 3+ factors, but slower)
#' mod_configural <- multipleGroup(dat, model, group = group, method = 'MHRM')
#' mod_metric <- multipleGroup(dat, model, group = group, invariance=c('slopes'), method = 'MHRM')
#' mod_fullconstrain <- multipleGroup(dat, model, group = group, method = 'MHRM',
#'                              invariance=c('slopes', 'intercepts'))
#'
#' anova(mod_metric, mod_configural)
#' anova(mod_fullconstrain, mod_metric)
#'
#' ############
#' # polytomous item example
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' d <- cbind(d, d-1, d-2)
#' itemtype <- rep('graded', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- c(rep('D1', N), rep('D2', N))
#' model <- 'F1 = 1-15'
#'
#' mod_configural <- multipleGroup(dat, model, group = group)
#' plot(mod_configural)
#' plot(mod_configural, type = 'SE')
#' plot(mod_configural, type = 'gen.difficulty')
#' itemplot(mod_configural, 1)
#' itemplot(mod_configural, 1, type = 'info')
#' plot(mod_configural, type = 'trace') # messy, score function typically better
#' plot(mod_configural, type = 'itemscore')
#'
#' fs <- fscores(mod_configural, full.scores = FALSE)
#' head(fs[["D1"]])
#' fscores(mod_configural, method = 'EAPsum', full.scores = FALSE)
#'
#' # constrain slopes within each group to be equal (but not across groups)
#' model2 <- 'F1 = 1-15
#'            CONSTRAIN = (1-15, a1)'
#' mod_configural2 <- multipleGroup(dat, model2, group = group)
#' plot(mod_configural2, type = 'SE')
#' plot(mod_configural2, type = 'RE')
#' itemplot(mod_configural2, 10)
#'
#' ############
#' ## empirical histogram example (normal and bimodal groups)
#' set.seed(1234)
#' a <- matrix(rlnorm(50, .2, .2))
#' d <- matrix(rnorm(50))
#' ThetaNormal <- matrix(rnorm(2000))
#' ThetaBimodal <- scale(matrix(c(rnorm(1000, -2), rnorm(1000,2)))) #bimodal
#' Theta <- rbind(ThetaNormal, ThetaBimodal)
#' dat <- simdata(a, d, 4000, itemtype = '2PL', Theta=Theta)
#' group <- rep(c('G1', 'G2'), each=2000)
#'
#' EH <- multipleGroup(dat, 1, group=group, dentype="empiricalhist", invariance = colnames(dat))
#' coef(EH, simplify=TRUE)
#' plot(EH, type = 'empiricalhist', npts = 60)
#'
#' # DIF test for item 1
#' EH1 <- multipleGroup(dat, 1, group=group, dentype="empiricalhist", invariance = colnames(dat)[-1])
#' anova(EH, EH1)
#'
#' #--------------------------------
#' # Mixture model (no prior group variable specified)
#'
#' set.seed(12345)
#' nitems <- 20
#' a1 <- matrix(.75, ncol=1, nrow=nitems)
#' a2 <- matrix(1.25, ncol=1, nrow=nitems)
#' d1 <- matrix(rnorm(nitems,0,1),ncol=1)
#' d2 <- matrix(rnorm(nitems,0,1),ncol=1)
#' itemtype <- rep('2PL', nrow(a1))
#' N1 <- 500
#' N2 <- N1*2 # second class twice as large
#'
#' dataset1 <- simdata(a1, d1, N1, itemtype)
#' dataset2 <- simdata(a2, d2, N2, itemtype)
#' dat <- rbind(dataset1, dataset2)
#' # group <- c(rep('D1', N1), rep('D2', N2))
#'
#' # Mixture Rasch model (Rost, 1990)
#' models <- 'F1 = 1-20
#'            CONSTRAIN = (1-20, a1)'
#' mod_mix <- multipleGroup(dat, models, dentype = 'mixture-2', GenRandomPars = TRUE)
#' coef(mod_mix, simplify=TRUE)
#' summary(mod_mix)
#' plot(mod_mix)
#' plot(mod_mix, type = 'trace')
#' plot(mod_mix, type = 'gen.difficulty')
#' itemplot(mod_mix, 1, type = 'info')
#'
#' head(fscores(mod_mix)) # theta estimates
#' head(fscores(mod_mix, method = 'classify')) # classification probability
#' itemfit(mod_mix)
#'
#' # Above works fine, but its generally a good idea to evaluate models
#' # with multiple random starting values in case local maximums are an issue.
#' # To do this use the argument "nruns", which returns the best model
#' if(interactive()) mirtCluster()
#' mod_mix <- multipleGroup(dat, models, dentype = 'mixture-2', nruns=5)
#' mod_mix
#'
#' # For obtaining isolated estimates within each mixture, use extract.group()
#' #   to construct single-group extractions of the mixtures
#' mix1 <- extract.group(mod_mix, group = "MIXTURE_1")
#' mix2 <- extract.group(mod_mix, group = "MIXTURE_2")
#'
#' # EAP estimates per mixture group, ignoring the original mixture structure.
#' #   Used to demonstrate the behaviour of how the individuals would have
#' #   been scored if they (deterministically) belonged to one class
#' data.frame(EAP_mix1=unname(fscores(mix1)),
#'            EAP_mix2=unname(fscores(mix2)),
#'            EAP=unname(fscores(mod_mix))) |> head()
#'
#' ############
#' # Mixture 2PL model
#' mod_mix2 <- multipleGroup(dat, 1, dentype = 'mixture-2', nruns=5)
#' anova(mod_mix, mod_mix2)
#' coef(mod_mix2, simplify=TRUE)
#' itemfit(mod_mix2)
#'
#' # Compare to single group
#' mod <- mirt(dat)
#' anova(mod, mod_mix2)
#'
#' ########################################
#' # Zero-inflated 2PL IRT model (Wall, Park, and Moustaki, 2015)
#'
#' n <- 1000
#' nitems <- 20
#'
#' a <- rep(2, nitems)
#' d <- rep(c(-2,-1,0,1,2), each=nitems/5)
#' zi_p <- 0.2 # Proportion of people in zero class
#'
#' theta <- rnorm(n, 0, 1)
#' zeros <- matrix(0, n*zi_p, nitems)
#' nonzeros <- simdata(a, d, n*(1-zi_p), itemtype = '2PL',
#'                    Theta = as.matrix(theta[1:(n*(1-zi_p))]))
#' data <- rbind(nonzeros, zeros)
#'
#' # define class with extreme theta but fixed item parameters
#' zi2PL <- "F = 1-20
#'           START [MIXTURE_1] = (GROUP, MEAN_1, -100), (GROUP, COV_11, .00001),
#'                               (1-20, a1, 1.0), (1-20, d, 0)
#'           FIXED [MIXTURE_1] = (GROUP, MEAN_1), (GROUP, COV_11),
#'                               (1-20, a1), (1-20, d)"
#'
#' # define custom Theta integration grid that contains extreme theta + normal grid
#' technical <- list(customTheta = matrix(c(-100, seq(-6,6,length.out=61))))
#'
#' # fit ZIM-IRT
#' zi2PL.fit <- multipleGroup(data, zi2PL, dentype = 'mixture-2', technical=technical)
#' coef(zi2PL.fit, simplify=TRUE)
#'
#' # classification estimates
#' pi_hat <- fscores(zi2PL.fit, method = 'classify')
#' head(pi_hat)
#' tail(pi_hat)
#'
#' # EAP estimates (not useful for zip class)
#' fs <- fscores(zi2PL.fit)
#' head(fs)
#' tail(fs)
#'
#' ########################################
#' # Zero-inflated graded response model (Magnus and Garnier-Villarreal, 2022)
#'
#' n <- 1000
#' nitems <- 20
#'
#' a <- matrix(rlnorm(20,.2,.3))
#'
#' # for the graded model, ensure that there is enough space between the intercepts,
#' # otherwise closer categories will not be selected often (minimum distance of 0.3 here)
#' diffs <- t(apply(matrix(runif(20*4, .3, 1), 20), 1, cumsum))
#' diffs <- -(diffs - rowMeans(diffs))
#' d <- diffs + rnorm(20)
#'
#' zi_p <- 0.2 # Proportion of people in zero/lowest category class
#'
#' theta <- rnorm(n, 0, 1)
#' zeros <- matrix(0, n*zi_p, nitems)
#' nonzeros <- simdata(a, d, n*(1-zi_p), itemtype = 'graded',
#'                     Theta = as.matrix(theta[1:(n*(1-zi_p))]))
#' data <- rbind(nonzeros, zeros)
#'
#' # intercepts will be labelled as d1 through d4
#' apply(data, 2, table)
#'
#' # ignoring zero inflation (bad idea)
#' modGRM <- mirt(data)
#' coef(modGRM, simplify=TRUE)
#'
#' # Define class with extreme theta but fixed item parameters
#' #   For GRM in zero-inflated class the intercept values are arbitrary
#' #   as the model forces the responses all into the first category (hence,
#' #   spacing arbitrarily set to 1)
#' ziGRM <- "F = 1-20
#'           START [MIXTURE_1] = (GROUP, MEAN_1, -100), (GROUP, COV_11, .00001),
#'                               (1-20, a1, 1.0),
#'                               (1-20, d1, 2), (1-20, d2, 1), (1-20, d3, 0), (1-20, d4, -1)
#'           FIXED [MIXTURE_1] = (GROUP, MEAN_1), (GROUP, COV_11),
#'                               (1-20, a1),
#'                               (1-20, d1), (1-20, d2), (1-20, d3), (1-20, d4)"
#'
#' # define custom Theta integration grid that contains extreme theta + normal grid
#' technical <- list(customTheta = matrix(c(-100, seq(-6,6,length.out=61))))
#'
#' # fit zero-inflated GRM
#' ziGRM.fit <- multipleGroup(data, ziGRM, dentype = 'mixture-2', technical=technical)
#' coef(ziGRM.fit, simplify=TRUE)
#'
#' # classification estimates
#' pi_hat <- fscores(ziGRM.fit, method = 'classify')
#' head(pi_hat)
#' tail(pi_hat)
#'
#' # EAP estimates (not useful for zip class)
#' fs <- fscores(ziGRM.fit)
#' head(fs)
#' tail(fs)
#'
#' }
multipleGroup <- function(data, model = 1, group, itemtype = NULL,
                          invariance = '', method = 'EM',
                          dentype = 'Gaussian', itemdesign=NULL, item.formula = NULL,
                          nruns = 1, return_max = TRUE, GenRandomPars = FALSE,
                          verbose = interactive(), ...)
{
    Call <- match.call()
    dots <- list(...)
    if(nruns > 1) GenRandomPars <- TRUE
    mixed.design <- make.mixed.design(item.formula=item.formula,
                                      itemdesign=itemdesign, data=data)
    if(is.character(model)) model <- mirt.model(model)
    if(!is.null(dots$formula))
        stop('latent regression models not supported for multiple group yet', call.=FALSE) #TODO
    invariance[invariance %in% c("free_mean")] <- 'free_means'
    invariance[invariance %in% c("free_vars", 'free_variance', 'free_variances')] <- 'free_var'
    constrain <- dots$constrain
    invariance.check <- sum(invariance %in% c('free_means', 'free_var'))
    if(any(invariance == 'intercepts')) invariance.check <- invariance.check - 1L
    if(any(invariance == 'slopes')) invariance.check <- invariance.check - 1L
    if(any(invariance %in% colnames(data)))
        invariance.check <- invariance.check - 2L
    if(!is.null(dots$dentype))
        if(dots$dentype == "empiricalhist" && any(invariance.check))
            stop('freeing group parameters not meaningful when estimating empirical histograms',
                 call.=FALSE)
    if(invariance.check > 0L && length(constrain) == 0){
        warn <- TRUE
        if(is(model, 'mirt.model')){
            if(any(model$x[,1L] == 'CONSTRAINB'))
                warn <- FALSE
        }
        if(warn)
            stop('Model is not identified without further constraints (may require additional
                 anchoring items).', call.=FALSE)
    }
    if(grepl('mixture', dentype)) group <- rep('full', nrow(data))
    mods <- myLapply(1:nruns, function(x, ...) return(ESTIMATION(...)),
                     progress=verbose && nruns > 1L,
                     data=data, model=model, group=group, invariance=invariance, method=method,
                     itemtype=itemtype, dentype=dentype, mixed.design=mixed.design,
                     GenRandomPars=GenRandomPars, ...)
    is_model <- is(mods[[1]], 'MultipleGroupClass') ||
        is(mods[[1]], 'MixtureClass')
    if(is_model){
        for(i in 1:length(mods)) mods[[i]]@Call <- Call
    }
    if(!return_max){
        return(mods)
    } else {
        if(is_model){
            LL <- sapply(mods, function(x) x@Fit$logLik)
            if(verbose && nruns > 1L){
                cat('Model log-likelihoods:\n')
                print(round(LL, 4))
            }
            mods <- mods[[which(max(LL) == LL)[1L]]]
        }
    }
    if(!is_model) mods <- mods[[1L]]
    mods
}
