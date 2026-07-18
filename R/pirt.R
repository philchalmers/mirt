#' Projective IRT (PIRT) models
#'
#' Computes the projective IRT model parameters either using the
#' logistic kernel approximation approaches (Stucky et al. 2013;
#' Ip 2010/Doebler and Doebler, 2022) or the MML-PIRT approach by
#' Chalmers et al. (in review). These return information pertaining to a
#' lower-dimensional (often unidimensional) IRT model
#' that has marginalized one or more latent traits.
#' When the method is the MML-PIRT a working \code{mirt} object will be
#' returned, otherwise if the target is the logistic kernel
#' approximations then a list
#' of the PIRT estimates (and potentially their delta-method SEs) are returned.
#'
#' Logistic kernel approximation approaches are limited only to the
#' M2PL (and related models, like the M4PL) and MGRM family of
#' ordinal models, while MML-PIRT can be applied to any MIRT model supported
#' by the package.
#'
#' Standard error are computed when \code{FAapprox = TRUE} via the delta method
#' only when the parent model estimated the ACOV matrix, while the
#' MML-PIRT approach does not require this prerequisite
#' as the ACOV can be computed directly. By default,
#' MML-PIRT computes the associated ACOV matrix using the Oakes'
#' identity given then marginalized E-table.
#'
#' @param mod fitted model from \code{mirt}
#'
#' @param model type of PIRT model to fit to marginalized table
#'
#' @param itemtype type of model to fit to the projected marginalized factor.
#'   Defaults to the same type as the parent object \code{mod}.
#'   See \code{\link{mirt}} for details
#'
#' @param SE logical; compute the ACOV matrix to obtain standard errors
#'   (SE)?
#'
#' @param project which dimension to project to in the stored E-table. Default
#'   projects all dimensions to the first dimension (\code{project = 1}),
#'   as this is the assumed primary dimension in, for instance, \code{bfactor()}
#'
#' @param invariance type of group invariance to specify (see \code{\link{multipleGroup}}).
#'   Included here as the multiple-group PIRT approach needs to also calibrate the
#'   scale and location of the reference group, which is done silently internally whenever
#'   \code{'free_means'} or \code{'free_vars'} are specified
#'
#' @param estiamtor type of PIRT estimator to use. Default is \code{'MML'} to fit and
#'   return the MML-PIRT variable described by Chalmers et al. (in review).
#'   \code{'FA'} will use Stucky et al.'s (2012) factor analysis logistic
#'   approximation to project to the first dimension, and \code{'LKA'}
#'   will use the use logistic kernel approximation formula
#'   presented in Ip (2010)/Doebler and Doebler(2020), however these only project
#'   to the first dimension? The latter two estimators are also only
#'   applicable when input model is from the M*PL and MGRM family (see
#'   for polytomous MGRM derivation), and only return lists of the resulting
#'   parameters and SE estimates
#'
#'   Note that for methods other than \code{'MML'} the SEs will be returned
#'   only if the parent model  as a suitably estimated ACOV matrix to perform
#'   the delta method (e.g., \code{bfactor(..., SE=TRUE})
#'
#' @param IRTpars logical; report classical IRT parameterization and SEs?
#'   Only used when \code{FAapprox = TRUE}
#'
#' @param shortform character vector of item names, or numeric vector of item locations,
#'   used create a short-form by specifying which items to extract
#'   from the PIRT model. For instance, selecting one item from each specific
#'   factor has the property of identical marginal and conditional bifactors
#'   (see Stucky, Thissen, and Elden, 2013). Only supported when using the MML-PIRT
#'   approach as this returns the associated model object
#'
#' @param ... extra information passed to the estimation engine
#'
#' @references
#'
#' Chalmers, R. P., Falk, C. F., Reise, S. P. (in review).
#' A General Approach for Estimating Projective IRT Models.
#' Applied Psychological Measurement.
#'
#' Doebler, A., and Doebler, P. (2022). Rotate and Project: Measurement of the
#' Intended Concept with Unidimensional Item Response Theory.
#' from Multidimensional Ordinal Items.
#' Multivariate Behavioral Research, 57 (1), 40-56.
#'
#' Ip, H. (2010). Empirically indistinguishable multidimensional IRT and
#' locally dependent unidimensional item response models.
#' British Journal of Mathematical and Statistical Psychology, 63(2),395-416.
#'
#' Stucky, B. D., Thissen, D., & Orlando E., M. (2012). Using Logistic
#' Approximations of Marginal Trace Lines to Develop Short Assessments.
#' Applied Psychological Measurement, 37(1), 41-57.
#'
#' @export
#' @examples
#'
#' \donttest{
#'
#' # Table 3 in Ip (2010)
#' as <- cbind(rep(c(2.63, 1.63, 1.38), each=5),
#'             c(rep(3.07, 5), numeric(10)),
#'             c(numeric(5), rep(1.36, 5), numeric(5)),
#'             c(numeric(10), rep(.69,5)))
#' as
#'
#' set.seed(8675309)
#' d <- round(rnorm(15, 0, sd=1.5), 2)
#' nitems <- length(d)
#' nfact <- ncol(as)
#' dat <- simdata(as, d, 10000, itemtype='2PL')
#' itemstats(dat)
#'
#'
#' # bifactor model, storing E-table
#' mod <- bfactor(dat,
#'                "S1 = 1-5
#'                 S2 = 6-10
#'                 S3 = 11-15")
#' summary(mod)
#' coef(mod, simplify=TRUE)$items
#'
#' # PIRT model via MML-PIRT
#' pmod <- pirt(mod)
#' coef(pmod)
#' coef(pmod, printSE=TRUE)
#' coef(pmod, simplify=TRUE)$items
#'
#' # Logistic approximations (parameters only)
#' pirt(mod, estimator = 'FA')
#' pirt(mod, estimator = 'LKA')  # effectively the same
#'
#' # plots
#' plot(pmod)
#' plot(pmod, type = 'info')
#' plot(pmod, type = 'itemscore')
#' plot(pmod, type = 'infotrace')
#' itemplot(pmod, 1, type = 'info', CE=TRUE)
#'
#' # standardized loadings
#' summary(pmod)
#'
#' # factor scores (EAP). Not optimal as the estimates/SEs ignore
#' #   the marginal components
#' fscores(pmod, method = 'EAPsum', full.scores=FALSE) # EAPs for sum-scores
#'
#' # EAPs using Lord-Wingersky 2.0 approach (recommended if using sum-scores)
#' fscores(mod, method = 'EAPsum_2.0', full.scores=FALSE)
#'
#' # EAP scores from full response data (not generally recommended)
#' fscores(pmod, full.scores=FALSE) |> head()
#' fs <- fscores(pmod)
#'
#' # compare to bifactor scores (marginal vs conditional)
#' bfs <- fscores(mod, method = 'MAP')
#' cor(bfs[,1], fs)
#'
#' # item fit (silly, but possible)
#' itemfit(pmod)
#'
#' # three item shortform by taking the highest PIRT discrimination/slope terms
#' # per specific factor to ensure that the marginal and joint likelihood
#' # functions are comparable (see Stucky et al. 2013 for details)
#' slopes <- coef(pmod, simplify=TRUE)$items[,'a1']
#' slopes
#' namemax <- \(slp) names(slp)[which.max(slp)]
#' shortform <- c(namemax(slopes[1:5]),
#'                namemax(slopes[6:10]), namemax(slopes[11:15]))
#' pirt_short <- pirt(mod, shortform=shortform)
#' coef(pirt_short, simplify=TRUE)$items
#'
#' # EAP-sum for short form
#' fscores(pirt_short, method='EAPsum', full.scores=FALSE)
#'
#' ###########################
#' # same example, but with polytomous data
#'
#' # Table 3 in Ip (2010)
#' as <- cbind(rep(c(2.63, 1.63, 1.38), each=5),
#'             c(rep(3.07, 5), numeric(10)),
#'             c(numeric(5), rep(1.36, 5), numeric(5)),
#'             c(numeric(10), rep(.69,5)))
#' as
#'
#' set.seed(8675309)
#' diffs <- t(apply(matrix(runif(20*4, .5, 1), 20), 1, cumsum))
#' diffs <- -(diffs - rowMeans(diffs))
#' d <- diffs + rnorm(20)
#' nitems <- nrow(d)
#' nfact <- ncol(as)
#' dat <- simdata(as, d, 10000, itemtype='graded')
#' itemstats(dat)
#'
#' # bifactor model, storing E-table
#' mod <- bfactor(dat,
#'                "S1 = 1-5
#'                 S2 = 6-10
#'                 S3 = 11-15")
#' summary(mod)
#' coef(mod, simplify=TRUE)$items
#'
#' # PIRT model via MML-PIRT
#' pmod <- pirt(mod)
#' coef(pmod)
#' coef(pmod, printSE=TRUE)
#' coef(pmod, simplify=TRUE)$items
#'
#' # Logistic approximations (parameters only)
#' pirt(mod, estimator = 'FA')
#' pirt(mod, estimator = 'LKA')  # see Doebler and Doebler (2020)
#'
#' # plots
#' plot(pmod)
#' plot(pmod, type = 'info')
#' plot(pmod, type = 'trace')
#' plot(pmod, type = 'infotrace')
#' itemplot(pmod, 1, type = 'info', CE=TRUE)
#'
#' # standardized loadings
#' summary(pmod)
#'
#' # factor scores (EAP)
#' fscores(pmod, method = 'EAPsum', full.scores=FALSE) # EAPs for sum-scores
#'
#' # EAP scores from full response data
#' fs <- fscores(pmod)
#'
#' # compare to bifactor scores (marginal vs conditional)
#' bfs <- fscores(mod, method = 'MAP')
#' cor(bfs[,1], fs)
#'
#'
#' ###########################
#' # Two tier example projected to simple structure
#'
#' # simulate data
#' set.seed(1234)
#' a <- matrix(c(
#'     0,1,0.5,NA,NA,
#'     0,1,0.5,NA,NA,
#'     0,1,0.5,NA,NA,
#'     0,1,0.5,NA,NA,
#'     0,1,0.5,NA,NA,
#'     0,1,NA,0.5,NA,
#'     0,1,NA,0.5,NA,
#'     0,1,NA,0.5,NA,
#'     1,0,NA,0.5,NA,
#'     1,0,NA,0.5,NA,
#'     1,0,NA,0.5,NA,
#'     1,0,NA,NA,0.5,
#'     1,0,NA,NA,0.5,
#'     1,0,NA,NA,0.5,
#'     1,0,NA,NA,0.5,
#'     1,0,NA,NA,0.5),ncol=5,byrow=TRUE)
#'
#' d <- matrix(rnorm(16))
#' items <- rep('2PL', 16)
#'
#' sigma <- diag(5)
#' sigma[1,2] <- sigma[2,1] <- .4
#' dataset <- simdata(a,d,2000,itemtype=items,sigma=sigma)
#' itemstats(dataset)
#'
#' specific <- "
#' S1 = 1-5
#' S2 = 6-11
#' S3 = 12-16"
#'
#' model <- '
#'     G1 = 1-8
#'     G2 = 9-16
#'     COV = G1*G2'
#'
#' # quadpts dropped for faster estimation, but not as precise
#' simmod <- bfactor(dataset, specific, model, quadpts = 15, TOL = 1e-3)
#' coef(simmod, simplify=TRUE)
#'
#' # simple structure projection
#' pirt2 <- pirt(simmod, model=model, project=1:2)
#' coef(pirt2, simplify=TRUE)
#'
#' }
#'
pirt <- function(mod, model=1, project=1,
                 itemtype = extract.mirt(mod, 'itemtype'),
                 SE=TRUE, estimator = 'MML', invariance='',
                 IRTpars = FALSE, shortform = NULL, ...){

    stopifnot(estimator %in% c('MML', "FA", 'LKA'))

    if(estimator == 'LKA')
        return(Ip2010(mod=mod, IRTpars=IRTpars))
    if(estimator == 'FA'){
        vcov <- vcov(mod)
        ret <- if(length(vcov) == 1 && is.na(vcov))
            FA.projection(mod, IRTpars=IRTpars)
        else FA.projection_SE(mod, IRTpars=IRTpars)
    } else {
        if(model == 1 && length(project) > 1)
            stop('Must specify MIRT model when using marginalizing multiple columns')
        lst <- marginal_etable(mod, column=project)
        dat <- extract.mirt(mod, 'data')
        ngroups <- extract.mirt(mod, 'ngroups')
        customTheta <- lst$theta
        fixedEtable <- lst$etable
        if(!is.null(shortform)){
            SE <- FALSE
            pick <- if(is.character(shortform))
                which(colnames(dat) %in% shortform) else shortform
            if(length(itemtype) > length(shortform))
                itemtype <- itemtype[pick]
            dat <- dat[ , shortform, drop=FALSE]
            K <- extract.mirt(mod, 'K')
            cK <- c(0, cumsum(K))
            for(g in 1:ngroups){
                newfixed <- vector('list', length(shortform))
                i <- 1
                for(p in pick){
                    from <- cK[p] + 1
                    to <- cK[p + 1]
                    newfixed[[i]] <- fixedEtable[[g]][,from:to]
                    i <- i + 1
                }
                fixedEtable[[g]] <- do.call(cbind, newfixed)
            }
        }
        if(extract.mirt(mod, 'exploratory') && ncol(customTheta) == 1){
            itemloc <- extract.mirt(mod, 'itemloc')
            itemloc <- itemloc[-length(itemloc)]
            lowertab <- fixedEtable[[1]][,itemloc]
            cors <- apply(lowertab, 2, \(x) cor(customTheta, x))
            if(mean(cors > 0) > .5)
                customTheta <- customTheta[nrow(customTheta):1, ,drop=FALSE]
        }
        if(ngroups == 1){
            ret <- mirt(dat, model=model, itemtype=itemtype, optimizer='nlminb', SE=SE,
                      technical=list(customTheta=customTheta,
                                     fixedEtable=fixedEtable), verbose=FALSE, ...)
        } else {
            old_invariance <- invariance
            if(any(c('free_means', 'free_mean') %in% invariance))
                invariance[invariance %in% c('free_means', 'free_mean')] <- 'free_all_means'
            if(any(c('free_vars', 'free_var') %in% invariance))
                invariance[invariance %in% c('free_vars', 'free_var')] <- 'free_all_vars'
            if(any(c('free_all_means',  'free_all_vars') %in% invariance)){
                tmpmod <- multipleGroup(dat, model=model, group=extract.mirt(mod, 'group'),
                                        itemtype=itemtype, optimizer='nlminb', SE=FALSE, invariance=invariance,
                                        technical=list(customTheta=customTheta,
                                                       fixedEtable=fixedEtable), verbose=FALSE, ...)
                sv <- multipleGroup(dat, model=model, group=extract.mirt(mod, 'group'), pars='values',
                                    itemtype=itemtype, optimizer='nlminb', SE=SE, invariance=old_invariance,
                                    technical=list(customTheta=customTheta,
                                                   fixedEtable=fixedEtable), verbose=FALSE, ...)
                sv$value <- mod2values(tmpmod)$value
                ret <- multipleGroup(dat, model=model, group=extract.mirt(mod, 'group'), pars=sv,
                                     itemtype=itemtype, optimizer='nlminb', SE=SE, invariance=old_invariance,
                                     technical=list(customTheta=customTheta,
                                                    fixedEtable=fixedEtable), verbose=FALSE, ...)
            } else {
                ret <- multipleGroup(dat, model=model, group=extract.mirt(mod, 'group'),
                                     itemtype=itemtype, optimizer='nlminb', SE=SE, invariance=invariance,
                                     technical=list(customTheta=customTheta,
                                                    fixedEtable=fixedEtable), verbose=FALSE, ...)
            }
        }
    }
    ret
}

marginal_etable <- function(mod, column = 1){
    nitems <- extract.mirt(mod, 'nitems')
    ngroups <- extract.mirt(mod, 'ngroups')
    grp_etable <- vector('list', ngroups)
    for(g in 1:ngroups){
        Etable <- mod@Internals$Etable[[g]]$r1
        Theta <- mod@Internals$Theta
        theta <- vector('list', length(column))
        for(i in 1:length(theta))
            theta[[i]] <- unique(Theta[,column[i]])
        theta <- as.matrix(expand.grid(theta))
        etable <- matrix(0, nrow(theta), ncol(Etable))
        for(i in 1:nrow(theta)){
            pick <- colSums(theta[i,] == t(Theta[ ,column])) == length(column)
            etable[i, ] <- colSums(Etable[pick, ])
        }
        grp_etable[[g]] <- etable
    }
    list(theta=theta, etable=grp_etable)
}

FA.projection <- function(mod, IRTpars=FALSE){
    cfs <- coef(mod, simplify=TRUE)$items
    a1 <- cfs[,'a1']
    so <- summary(mod, verbose=FALSE)
    lam <- so$rotF
    lam1 <- lam[,1]
    sigma1 <- sqrt(1 - lam1^2)
    D <- 1.702
    a1star <- lam1/sigma1*D

    # Equation 9.16 (though this equation incorrectly has an
    # extra negative as they use the
    # form a1*theta - c rather than a1*theta + c)
    cs <- cfs[,grepl('d', colnames(cfs))]
    bstar <- -cs / a1
    if(!IRTpars){
        ds <- -bstar * a1star
        ret <- data.frame(a=unname(a1star), ds)
        if(ncol(ret) == 2) colnames(ret) <- c('a', 'd')
    } else {
        ret <- data.frame(alpha=unname(a1star), beta=unname(bstar))
    }
    as.matrix(ret)
}

FA.projection_SE <- function(mod, IRTpars=FALSE){

    FA.projection.item <- function(ipars, IRTpars=FALSE){
        lam <- mirt:::Lambdas(ipars, Names = NULL)
        as <- ipars[[1]]@par[grepl('^a', ipars[[1]]@parnames)]
        cs <- ipars[[1]]@par[grepl('^d', ipars[[1]]@parnames)]
        lam1 <- lam[,1]
        sigma1 <- sqrt(1 - lam1^2)
        D <- 1.702
        a1star <- lam1/sigma1*D

        # Equation 9.16 (though this equation incorrectly has an extra negative as they
        #   use the form a1*theta - c rather than a1*theta + c)
        bstar <- -cs / as[1]
        if(!IRTpars){
            ds <- -bstar * a1star
            ret <- c(a1star, ds)
        } else {
            ret <- c(a1star, bstar)
        }
        ret
    }

    fn_tpars <- function(par, ipars, index=1){
        ipars[[1]]@par <- par
        tpars <- FA.projection.item(ipars)
        tpars[index]
    }

    pars <- mod@ParObjects$pars
    acov <- vcov(mod)
    tpars <- FA.projection(mod)
    nitems <- length(pars) - 1
    SEs <- matrix(NA, nitems, ncol(tpars))
    for(i in 1:nitems){
        for(j in 1:ncol(tpars)){
            ipars <- pars[c(i,nitems+1)]
            par <- ipars[[1]]@par
            iacov <- mirt:::subset_vcov(ipars[[1]], acov)
            SEs[i,j] <- DeltaMethod(fn_tpars, par, iacov,
                                    ipars=ipars, index=j)$se
        }
        if(is.na(SEs[i,j])){
            SEs[i,j] <- DeltaMethod(fn_tpars, par, iacov,
                                    ipars=ipars, index=j)$se
        }
    }
    list(items=tpars, SE=SEs)
}

Ip2010 <- function(mod, IRTpars=FALSE){

    # Compare to Ip (2010) transformation
    Ip2010.item <- function(alpha_gen, alpha_grp, lam,
                            rho = 0, sigma_gen=1, sigma_grp=1){
        k <- 16 * sqrt(3) / (15 * pi)
        logit <- (k^2 * alpha_grp^2 * (1-rho^2) * sigma_grp^2 +1 )^{-1/2}
        alphastar <- logit *
            (alpha_gen + alpha_grp * rho * sigma_grp / sigma_gen)
        lamstar <- logit * lam
        betastar <- unname(-lamstar/alphastar)
        c(alpha=alphastar, beta=betastar)
    }

    nfact <- extract.mirt(mod, 'nfact')
    K <- extract.mirt(mod, 'K')
    # stopifnot(all(K == 2))
    cfs <- coef(mod, simplify=TRUE)$items
    nitems <- nrow(cfs)
    alpha_gen <- cfs[,'a1']
    alpha_grp <- cfs[,paste0('a', 2:nfact), drop=FALSE] |> rowSums()
    lam <- cfs[,grepl('d', colnames(cfs)), drop=FALSE]
    ret <- sapply(1:nitems, \(i)
                  Ip2010.item(alpha_gen=alpha_gen[i],
                              alpha_grp=alpha_grp[i], lam=lam[i,])) |> t()
    if(ncol(ret) == 2) c('alpha', 'beta')
    else colnames(ret) <- c('alpha', paste0('beta', 2:ncol(ret)-1))
    if(!IRTpars){
        ret[,-1] <- -ret[,1] * ret[,-1]
        if(ncol(ret) == 2) colnames(ret) <- c('a', 'd')
        else colnames(ret) <- c('a', paste0('d', 2:ncol(ret)-1))
    }
    ret <- as.matrix(ret)
    rownames(ret) <- extract.mirt(mod, 'itemnames')
    ret
}
