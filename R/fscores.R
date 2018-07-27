#' Compute factor score estimates (a.k.a, ability estimates, latent trait estimates, etc)
#'
#' Computes MAP, EAP, ML (Embretson & Reise, 2000), EAP for sum-scores (Thissen et al., 1995),
#' or WLE (Warm, 1989) factor scores with a multivariate normal
#' prior distribution using equally spaced quadrature. EAP scores for models with more than
#' three factors are generally not recommended since the integration grid becomes very large,
#' resulting in slower estimation and less precision if the \code{quadpts} are too low.
#' Therefore, MAP scores should be used instead of EAP scores for higher dimensional models.
#' Multiple imputation variants are possible for each estimator if a parameter
#' information matrix was computed, which are useful if the sample size/number of items were small.
#' As well, if the model contained latent regression predictors this information will
#' be used in computing MAP and EAP estimates (for these models, \code{full.scores=TRUE}
#' will always be used). Finally, plausible value imputation is also available, and will also account
#' for latent regression predictor effects.
#'
#' The function will return either a table with the computed scores and standard errors,
#' the original data matrix with scores appended to the rightmost column, or the scores only. By
#' default the latent means and covariances are determined from the estimated object,
#' though these can be overwritten. Iterative estimation methods can be estimated
#' in parallel to decrease estimation times if a \code{\link{mirtCluster}} object is available.
#'
#' If the input object is a discrete latent class object estimated from \code{\link{mdirt}}
#' then the returned results will be with respect to the posterior classification for each
#' individual. The method inputs for \code{'DiscreteClass'} objects may only be \code{'EAP'},
#' for posterior classification of each response pattern, or \code{'EAPsum'} for posterior
#' classification based on the raw sum-score. For more information on these algorithms refer to
#' the \code{mirtCAT} package and the associated JSS paper (Chalmers, 2016).
#'
#'
#' @aliases fscores
#' @param object a computed model object of class \code{SingleGroupClass},
#'   \code{MultipleGroupClass}, or \code{DiscreteClass}
#' @param full.scores if \code{FALSE} then a summary table with
#'   factor scores for each unique pattern is displayed. Otherwise, a matrix of factor scores
#'   for each response pattern in the data is returned (default)
#' @param rotate prior rotation to be used when estimating the factor scores. See
#'   \code{\link{summary-method}} for details. If the object is not an exploratory model
#'   then this argument is ignored
#' @param Target target rotation; see \code{\link{summary-method}} for details
#' @param plausible.draws number of plausible values to draw for future researchers
#'   to perform secondary analyses of the latent trait scores. Typically used in conjunction
#'   with latent regression predictors (see \code{\link{mirt}} for details), but can
#'   also be generated when no predictor variables were modeled. If \code{plausible.draws}
#'   is greater than 0 a list of plausible values will be returned
#' @param plausible.type type of plausible values to obtain. Can be either \code{'normal'} (default)
#'   to use a normal approximation based on the ACOV matrix, or \code{'MH'} to obtain Metropolis-Hastings
#'   samples from the posterior (silently passes object to \code{\link{mirt}}, therefore arguments like
#'   \code{technical} can be supplied to increase the number of burn-in draws and discarded samples)
#' @param method type of factor score estimation method. Can be:
#' \itemize{
#'     \item \code{"EAP"} for the expected a-posteriori (default)
#'     \item \code{"MAP"} for the maximum a-posteriori (i.e, Bayes modal)
#'     \item \code{"ML"} for maximum likelihood
#'     \item \code{"WLE"} for weighted likelihood estimation
#'     \item \code{"EAPsum"} for the expected a-posteriori for each sum score
#'     \item \code{"plausible"} for a single plausible value imputation for each case.
#'       This is equivalent to setting \code{plausible.draws = 1}
#'     \item \code{"classify"} for the posteriori classification probabilities (only
#'       applicable when the input model was of class \code{MixtureClass})
#'  }
#' @param quadpts number of quadratures to use per dimension. If not specified, a suitable
#'   one will be created which decreases as the number of dimensions increases
#'   (and therefore for estimates such as EAP, will be less accurate). This is determined from
#'   the switch statement
#'   \code{quadpts <- switch(as.character(nfact), '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, 3)}
#' @param theta_lim lower and upper range to evaluate latent trait integral for each dimension. If
#'   omitted, a range will be generated automatically based on the number of dimensions
#' @param max_theta the maximum/minimum value any given factor score estimate will achieve using
#'   any modal estimator method (e.g., MAP, WLE, ML)
#' @param mean a vector for custom latent variable means. If NULL, the default for 'group' values
#'   from the computed mirt object will be used
#' @param cov a custom matrix of the latent variable covariance matrix. If NULL, the default for
#'   'group' values from the computed mirt object will be used
#' @param MI a number indicating how many multiple imputation draws to perform. Default is 0,
#'   indicating that no MI draws will be performed
#' @param use_dentype_estimate logical; if the density of the latent trait was estimated in the model
#'   (e.g., via Davidian curves or empirical histograms), should this information be used to compute
#'   the latent trait estimates? Only applicable for EAP-based estimates (EAP, EAPsum, and plausible)
#' @param response.pattern an optional argument used to calculate the factor scores and standard
#'   errors for a given response vector or matrix/data.frame
#' @param append_response.pattern logical; should the inputs from \code{response.pattern} also be
#'   appended to the factor score output?
#' @param na.rm logical; remove rows with any missing values? This is generally not required due to
#'   the nature of computing factors scores, however for the "EAPsum" method this may be necessary
#'   to ensure that the sum-scores correspond to the same composite score
#' @param returnER logical; return empirical reliability (also known as marginal reliability)
#'   estimates as a numeric values?
#' @param return.acov logical; return a list containing covariance matrices instead of factors
#'   scores? \code{impute = TRUE} not supported with this option
#' @param start a matrix of starting values to use for iterative estimation methods. Default
#'   will start at a vector of 0's for each response pattern, or will start at the EAP
#'   estimates (unidimensional models only). Must be in the form that matches
#'   \code{full.scores = FALSE} (mostly used in the \code{mirtCAT} package)
#' @param full.scores.SE logical; when \code{full.scores == TRUE}, also return the
#'   standard errors associated with each respondent? Default is \code{FALSE}
#' @param verbose logical; print verbose output messages?
#' @param QMC logical; use quasi-Monte Carlo integration? If \code{quadpts} is omitted the
#'   default number of nodes is 5000
#' @param custom_den a function used to define the integration density (if required). The NULL default
#'   assumes that the multivariate normal distribution with the 'GroupPars' hyper-parameters are
#'   used. At the minimum must be of the form:
#'
#'   \code{function(Theta, ...)}
#'
#'   where Theta is a matrix of latent trait values (will be a grid of values
#'   if \code{method == 'EAPsum'} or \code{method == 'EAP'}, otherwise Theta will have only 1 row).
#'   Additional arguments may included and are caught through the \code{fscores(...)} input. The
#'   function \emph{must} return a numeric vector of density weights (one for each row in Theta)
#' @param custom_theta a matrix of custom integration nodes to use instead of the default, where
#'   each column corresponds to the respective dimension in the model
#' @param min_expected when computing goodness of fit tests when \code{method = 'EAPsum'}, this value is used
#'   to collapse across the conditioned total scores until the expected values are greater than this value. Note
#'   that this only affect the goodness of fit tests and not the returned EAP for sum scores table
#' @param converge_info logical; include a column in the return objects containing a logical for each
#'   response pattern indicating whether a maximum value was found (not relevant non-iterative methods,
#'   such as EAP and EAPsum). Value is a reflection of the \code{code} element from \code{\link{nlm}}
#'   (e.g., 1 indicates convergence)
#' @param ... additional arguments to be passed to \code{nlm}
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords factor.scores
#' @seealso \code{\link{averageMI}}
#' @export fscores
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Chalmers, R. P. (2016). Generating Adaptive and Non-Adaptive Test Interfaces for
#' Multidimensional Item Response Theory Applications. \emph{Journal of Statistical Software, 71}(5),
#' 1-39. \doi{10.18637/jss.v071.i05}
#'
#' Embretson, S. E. & Reise, S. P. (2000). Item Response Theory for Psychologists. Erlbaum.
#'
#' Thissen, D., Pommerich, M., Billeaud, K., & Williams, V. S. L. (1995).
#' Item Response Theory for Scores on Tests Including Polytomous Items with Ordered Responses.
#' \emph{Applied Psychological Measurement, 19}, 39-49.
#'
#' Warm, T. A. (1989). Weighted likelihood estimation of ability in item response theory.
#' \emph{Psychometrika, 54}, 427-450.
#'
#' @examples
#'
#' mod <- mirt(Science, 1)
#' tabscores <- fscores(mod, full.scores = FALSE)
#' head(tabscores)
#'
#' \dontrun{
#' fullscores <- fscores(mod)
#' fullscores_with_SE <- fscores(mod, full.scores.SE=TRUE)
#' head(fullscores)
#' head(fullscores_with_SE)
#'
#' #change method argument to use MAP estimates
#' fullscores <- fscores(mod, method='MAP')
#' head(fullscores)
#'
#' #calculate MAP for a given response vector
#' fscores(mod, method='MAP', response.pattern = c(1,2,3,4))
#' #or matrix
#' fscores(mod, method='MAP', response.pattern = rbind(c(1,2,3,4), c(2,2,1,3)))
#'
#' # return only the scores and their SEs
#' fscores(mod, method='MAP', response.pattern = c(1,2,3,4),
#'   append_response.pattern=FALSE)
#'
#' #use custom latent variable properties (diffuse prior for MAP is very close to ML)
#' fscores(mod, method='MAP', cov = matrix(1000), full.scores = FALSE)
#' fscores(mod, method='ML', full.scores = FALSE)
#'
#' # EAPsum table of values based on total scores
#' fscores(mod, method = 'EAPsum', full.scores = FALSE)
#'
#' #WLE estimation, run in parallel using available cores
#' mirtCluster()
#' head(fscores(mod, method='WLE', full.scores = FALSE))
#'
#' #multiple imputation using 30 draws for EAP scores. Requires information matrix
#' mod <- mirt(Science, 1, SE=TRUE)
#' fs <- fscores(mod, MI = 30)
#' head(fs)
#'
#' # plausible values for future work
#' pv <- fscores(mod, plausible.draws = 5)
#' lapply(pv, function(x) c(mean=mean(x), var=var(x), min=min(x), max=max(x)))
#'
#' ## define a custom_den function. EAP with a uniform prior between -3 and 3
#' fun <- function(Theta, ...) as.numeric(dunif(Theta, min = -3, max = 3))
#' head(fscores(mod, custom_den = fun))
#'
#' # custom MAP prior: standard truncated normal between 5 and -2
#' library(msm)
#' # need the :: scope for parallel to see the function (not require if no mirtCluster() defined)
#' fun <- function(Theta, ...) msm::dtnorm(Theta, mean = 0, sd = 1, lower = -2, upper = 5)
#' head(fscores(mod, custom_den = fun, method = 'MAP', full.scores = FALSE))
#'
#'}
fscores <- function(object, method = "EAP", full.scores = TRUE, rotate = 'oblimin', Target = NULL,
                    response.pattern = NULL, append_response.pattern=TRUE, na.rm=FALSE,
                    plausible.draws = 0, plausible.type = 'normal', quadpts = NULL,
                    returnER = FALSE, return.acov = FALSE, mean = NULL, cov = NULL, verbose = TRUE,
                    full.scores.SE = FALSE, theta_lim = c(-6,6), MI = 0, use_dentype_estimate=FALSE,
                    QMC = FALSE, custom_den = NULL, custom_theta = NULL, min_expected = 1,
                    converge_info = FALSE, max_theta = 20, start = NULL, ...)
{
    if(!is(object, 'DiscreteClass')){
        if(QMC && is.null(quadpts)) quadpts <- 5000
        if(is.null(quadpts))
            quadpts <- switch(as.character(object@Model$nfact),
                              '1'=61, '2'=31, '3'=15, '4'=9, '5'=7, 3)
    } else quadpts <- 1
    if(method == 'plausible'){
        plausible.draws <- 1
        method <- 'EAP'
    }
    if(any(extract.mirt(object, 'itemtype') == 'spline') && !(method %in% c('EAP', 'EAPsum')))
        stop('Only EAP and EAPsum method supported when spline items are modeled', call.=FALSE)
    if(returnER) full.scores <- FALSE
    if(na.rm)
        object <- removeMissing(object)
    if(theta_lim[1L] < -max_theta) theta_lim[1L] <- -max_theta
    if(theta_lim[2L] > max_theta) theta_lim[2L] <- max_theta
    if(is(object, 'DiscreteClass') && plausible.draws > 0L){
        fs <- fscores(object)
        ret <- lapply(seq_len(plausible.draws), function(ind, fs){
            mat <- matrix(0L, nrow(fs), ncol(fs))
            for(i in seq_len(ncol(fs))){
                if(all(fs[,i] > 1-1e-10)){
                    mat[,i] <- 1L
                } else if(all(fs[,i] < 1e-10)){
                    mat[,i] <- 0L
                } else {
                    mat[,i] <- sapply(fs[,i], function(prob)
                        sample(c(0L,1L), 1L, prob = c(1-prob, prob)))
                }
            }
            mat
        }, fs=fs)
        return(ret)
    }
    ret <- fscores.internal(object=object, rotate=rotate, full.scores=full.scores, method=method,
                            quadpts=quadpts, response.pattern=response.pattern, QMC=QMC,
                            verbose=verbose, returnER=returnER, gmean=mean, gcov=cov,
                            theta_lim=theta_lim, MI=MI, converge_info=converge_info,
                            full.scores.SE=full.scores.SE, return.acov=return.acov,
                            plausible.draws = plausible.draws, custom_den=custom_den,
                            custom_theta=custom_theta, Target=Target, min_expected=min_expected,
                            plausible.type=plausible.type, max_theta=max_theta, start=start,
                            use_dentype_estimate=use_dentype_estimate,
                            append_response.pattern=append_response.pattern, ...)
    ret
}
