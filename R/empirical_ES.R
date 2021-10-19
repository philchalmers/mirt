#' Empirical effect sizes based on latent trait estimates
#'
#' Computes effect size measures of differential item functioning and differential
#' test/bundle functioning based on expected scores from Meade (2010).
#' Item parameters from both reference and focal group are used in conjunction with
#' focal group empirical theta estimates (and an assumed normally distributed theta)
#' to compute expected scores.
#'
#' @section DIF:
#'
#' The default \code{DIF = TRUE} produces several effect sizes indices at the item level.
#' Signed indices allow DIF favoring the focal group at one point on the theta
#' distribution to cancel DIF favoring the reference group at another point on the theta
#' distribution. Unsigned indices take the absolute value before summing or averaging,
#' thus not allowing cancellation of DIF across theta.
#'
#' \describe{
#'   \item{SIDS}{Signed Item Difference in the Sample. The average difference in expected scores
#' across the focal sample using both focal and reference group item parameters.}
#'   \item{UIDS}{Unsigned Item Difference in the Sample. Same as SIDS except absolute value of
#' expected scores is taken prior to averaging across the sample.}
#'   \item{D-Max}{The maximum difference in expected scores in the sample.}
#'   \item{ESSD}{Expected Score Standardized Difference. Cohen's D for difference in expected scores.}
#'   \item{SIDN}{Signed Item Difference in a Normal distribution. Identical to SIDS but
#' averaged across a normal distribution rather than the sample.}
#'   \item{UIDN}{Unsigned Item Difference in a Normal distribution. Identical to UIDS but
#' averaged across a normal distribution rather than the sample.}
#' }
#'
#' @section DBF/DTF:
#'
#' \code{DIF = FALSE} produces a series of test/bundle-level indices that are based on item-level
#' indices.
#'
#' \describe{
#'   \item{STDS}{Signed Test Differences in the Sample. The sum of the SIDS across items.}
#'   \item{UTDS}{Unsigned Test Differences in the Sample. The sum of the UIDS across items.}
#'   \item{Stark's DTFR}{Stark's version of STDS using a normal distribution rather than
#' sample estimated thetas.}
#'   \item{UDTFR}{Unsigned Expected Test Scores Differences in the Sample. The difference
#' in observed summed scale scores expected, on average, across a hypothetical focal
#' group with a normally distributed theta, had DF been uniform in nature for all items}
#'   \item{UETSDS}{Unsigned Expected Test Score Differences in the Sample.
#' The hypothetical difference expected scale scores that would have been present if
#' scale-level DF had been uniform across respondents (i.e., always favoring the
#' focal group).}
#'   \item{UETSDN}{Identical to UETSDS but computed using a normal distribution.}
#'   \item{Test D-Max}{Maximum expected test score differences in the sample.}
#'   \item{ETSSD}{Expected Test Score Standardized Difference. Cohen's D for expected
#' test scores.}
#' }
#'
#' @param mod a multipleGroup object which estimated only 2 groups. The first group in this object
#'   is assumed to be the reference group by default (i.e., \code{ref.group = 1}), which conforms to the
#'    \code{invariance} arguments in \code{\link{multipleGroup}}
# @param ref.group default group to use for the reference group (default is 1, conforming to the
#   default setup in \code{\link{multipleGroup}}), and is passed to \code{\link{extract.group}};
#   hence, can be either a number to extract or the name of the grouping variable
#' @param focal_items a numeric vector indicating which items to include the tests. The
#'   default uses all of the items. Selecting fewer items will result in tests of
#'   'differential bundle functioning' when \code{DIF = FALSE}
#' @param npts number of points to use in the integration. Default is 61
#' @param theta_lim lower and upper limits of the latent trait (theta) to be evaluated, and is
#'   used in conjunction with \code{npts}
#' @param Theta.focal an optional matrix of Theta values from the focal group to be evaluated. If not supplied
#'   the default values to \code{\link{fscores}} will be used in conjunction with the \code{...}
#'   arguments passed
#' @param DIF logical; return a data.frame of item-level imputation properties? If \code{FALSE},
#'   only DBF and DTF statistics will be reported
#' @param plot logical; plot expected scores of items/test where expected scores are computed
#'  using focal group thetas and both focal and reference group item parameters
#' @param type type of objects to draw in \code{lattice}; default plots both points and lines
#' @param par.strip.text plotting argument passed to \code{\link{lattice}}
#' @param par.settings plotting argument passed to \code{\link{lattice}}
#' @param ... additional arguments to be passed to \code{\link{fscores}} and \code{\link{xyplot}}
#'
#' @author Adam Meade, with contributions by Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @references
#' Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
#' Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
#' \doi{10.18637/jss.v048.i06}
#'
#' Meade, A. W. (2010). A taxonomy of effect size measures for the differential functioning
#' of items and scales. \emph{Journal of Applied Psychology, 95}, 728-743. \doi{10.1037/a0018966}
#' @export empirical_ES
#' @examples
#' \dontrun{
#'
#' #no DIF
#' set.seed(12345)
#' a <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d <- matrix(rnorm(15,0,.7),ncol=1)
#' itemtype <- rep('2PL', nrow(a))
#' N <- 1000
#' dataset1 <- simdata(a, d, N, itemtype)
#' dataset2 <- simdata(a, d, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#'
#' # ensure 'Ref' is the first group (and therefore reference group during estimation)
#' group <- factor(c(rep('Ref', N), rep('Focal', N)), levels = c('Ref', 'Focal'))
#'
#' mod <- multipleGroup(dat, 1, group = group,
#'    invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#'
#' empirical_ES(mod)
#' empirical_ES(mod, DIF=FALSE)
#' empirical_ES(mod, DIF=FALSE, focal_items = 10:15)
#'
#' empirical_ES(mod, plot=TRUE)
#' empirical_ES(mod, plot=TRUE, DIF=FALSE)
#'
#' ###---------------------------------------------
#' # DIF
#' set.seed(12345)
#' a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
#' d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
#' a2[10:15,] <- a2[10:15,] + rnorm(6, 0, .3)
#' d2[10:15,] <- d2[10:15,] + rnorm(6, 0, .3)
#' itemtype <- rep('dich', nrow(a1))
#' N <- 1000
#' dataset1 <- simdata(a1, d1, N, itemtype)
#' dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
#' dat <- rbind(dataset1, dataset2)
#' group <- factor(c(rep('Ref', N), rep('Focal', N)), levels = c('Ref', 'Focal'))
#'
#' mod <- multipleGroup(dat, 1, group = group,
#'    invariance = c(colnames(dat)[1:5], 'free_means', 'free_var'))
#' coef(mod, simplify=TRUE)
#'
#' empirical_ES(mod)
#' empirical_ES(mod, DIF = FALSE)
#' empirical_ES(mod, plot=TRUE)
#' empirical_ES(mod, plot=TRUE, DIF=FALSE)
#'
#' }
empirical_ES <- function(mod, Theta.focal = NULL, focal_items = 1L:extract.mirt(mod, 'nitems'),
                 DIF = TRUE, npts = 61, theta_lim=c(-6,6), plot=FALSE, type = 'b',
                 par.strip.text = list(cex = 0.7),
                 par.settings = list(strip.background = list(col = '#9ECAE1'),
                                     strip.border = list(col = "black")), ...){
    stopifnot(extract.mirt(mod, 'nfact') == 1L)
    stopifnot(extract.mirt(mod, 'ngroups') == 2L)
    ref.group <- 1
    ref <- extract.group(mod, ref.group)
    focal <- extract.group(mod, ifelse(ref.group == 1, 2, 1))
    focal_select <- extract.mirt(mod, 'group') != extract.mirt(mod, 'groupNames')[ref.group]
    if(is.null(Theta.focal)){
        Theta <- fscores(mod, full.scores = TRUE, full.scores.SE = FALSE,
                         leave_missing=TRUE, ...)
        Theta.focal <- Theta[focal_select, , drop = FALSE]
    } else Theta.focal <- as.matrix(Theta.focal)
    if(sum(focal_select) != nrow(Theta.focal))
        stop('Theta elements do not match the number of individuals in the focal group')

    ############# helper function -  Cohen D ###########
    f.cohen.d <- function (vector.1,vector.2){
      f.mean.v1 <- mean(vector.1)
      f.mean.v2 <- mean(vector.2)
      f.sd.v1 <- sd(vector.1)
      f.sd.v2 <- sd(vector.2)
      f.n.1 <- length(vector.1)
      f.n.2 <- length(vector.2)
      f.numerator <- (f.sd.v1^2)*(f.n.1 - 1) + (f.sd.v2^2)*(f.n.2 - 1)
      f.denom <- f.n.1 + f.n.2 - 2
      f.sd.pooled <- sqrt(f.numerator / f.denom)
      f.return <- (f.mean.v1 - f.mean.v2) / f.sd.pooled
      return(f.return)
    }
    ############# helper function -  Max D ###########
    get.max.D <- function(f.d){
      f.d.abs <- abs(f.d)
      f.max.D.location <- which.max(f.d.abs) #which.max returns location
      f.max.D <- f.d[f.max.D.location]
      f.max.D.theta <-  Theta.focal[f.max.D.location] # this is a global
      f.df.d.max <- c(f.max.D.theta,f.max.D)
      return(f.df.d.max)
    }

    # item level DF
    focal.theta.mean.obs <- mean(Theta.focal)
    #### quadrature nodes ####
    theta.normal   <- seq(theta_lim[1],theta_lim[2],length=npts)
    theta.den   <- dnorm(theta.normal,mean=focal.theta.mean.obs, sd=1)
    theta.density <- theta.den / sum(theta.den)
    nitems <- length(focal_items)
    list.item_ES_foc.obs <- list()
    list.item_ES_ref.obs <- list()
    list.item_ES_foc.nrm <- list()
    list.item_ES_ref.nrm <- list()
    ###### compute the expected scores (ES)
    for(i in 1:nitems){
      foc.extract<-extract.item(focal, i)
      ref.extract<-extract.item(ref, i)
      foc.ES.obs <- expected.item(foc.extract,Theta.focal)
      ref.ES.obs <- expected.item(ref.extract,Theta.focal)
      foc.ES.nrm <- expected.item(foc.extract,theta.normal)
      ref.ES.nrm <- expected.item(ref.extract,theta.normal)
      list.item_ES_foc.obs[[i]] <- foc.ES.obs
      list.item_ES_ref.obs[[i]] <- ref.ES.obs
      list.item_ES_foc.nrm[[i]] <- foc.ES.nrm
      list.item_ES_ref.nrm[[i]] <- ref.ES.nrm
    }
    rm(foc.extract,ref.extract,foc.ES.obs,ref.ES.obs,foc.ES.nrm,ref.ES.nrm)
    # put these in dataframe
    df.ref.obs <- do.call("cbind",list.item_ES_ref.obs)
    df.foc.obs <- do.call("cbind",list.item_ES_foc.obs)
    df.ref.nrm <- do.call("cbind",list.item_ES_ref.nrm)
    df.foc.nrm <- do.call("cbind",list.item_ES_foc.nrm)
    #### means for each item in dataframe
    mean.ES.foc <- colMeans(df.foc.obs)
    mean.ES.ref <- colMeans(df.ref.obs)
    ### dataframe of observed difference scores
    df.dif.obs <- df.foc.obs - df.ref.obs
    df.abs.dif.obs <- abs(df.dif.obs)
    SIDS <- colMeans(df.dif.obs)
    UIDS <- colMeans(df.abs.dif.obs)
    #Item level
    ################## MAX D section ################
    item.max.D <- apply(df.dif.obs,2,get.max.D)
    mat.item.max.d <- t(item.max.D)
    colnames(mat.item.max.d) <- c('theta of max D',"max D")
    mat.item.max.d
    list.item_CohenD.obs <- list()
    for(i in 1:nitems){
      list.item_CohenD.obs[i] <- f.cohen.d(df.foc.obs[,i],df.ref.obs[,i])
    }
    ESSD <- unlist(list.item_CohenD.obs) #vector
    ############# normal distribution #########
    df.dif.nrm <- df.foc.nrm - df.ref.nrm        # DF in ES at each level of theta
    weighted.dif.nrm <- apply(df.dif.nrm,2, function(x) x*theta.density)
    SIDN <- colSums(weighted.dif.nrm)           #now sum these
    df.abs.dif.nrm <- abs(df.dif.nrm)            # abs(DF in ES at each level of theta)
    weighted.dif.abs.nrm <- apply(df.abs.dif.nrm,2, function(x) x*theta.density)
    UIDN <- colSums(weighted.dif.abs.nrm)
    df.item.output <- data.frame(SIDS,UIDS,SIDN,UIDN,ESSD,mat.item.max.d,mean.ES.foc,mean.ES.ref)
    row.names(df.item.output)<-paste0("item.",1:nrow(df.item.output))
    class(df.item.output) <- c('mirt_df', 'data.frame')
    if(!plot && DIF) return(df.item.output[focal_items, ])

    ##################DTF####################
    STDS <- sum(SIDS)
    UTDS <- sum(UIDS)
    ##### UETSDS
    ETS.foc.obs <- rowSums(df.foc.obs)  ### test ES (avg across items)
    ETS.ref.obs <- rowSums(df.ref.obs)  ### test ES (avg across items)
    ETS.dif.obs <- ETS.foc.obs - ETS.ref.obs  ### df in test scores
    UETSDS <- mean(abs(ETS.dif.obs))    ### no cancel across theta, yes across items
    #####
    #cohen D and D max
    ETSSD <- f.cohen.d(ETS.foc.obs, ETS.ref.obs) # cohen's D
    test.Dmax <- get.max.D(ETS.dif.obs)
    ############# test level
    ETS.abs.dif.nrm <- rowSums(df.abs.dif.nrm) # [abs(DF in ES at each level of theta)] summed across items
    UDTFR <- sum(ETS.abs.dif.nrm*theta.density)
    UDTFR.b <- sum(UIDN)
    ETS.dif.nrm <- rowSums(df.dif.nrm)  ### ETS dif at each theta, summed across items
    Starks.DTFR <- sum(ETS.dif.nrm*theta.density)
    Starks.DTFR.b <- sum(SIDN)
    ETS.foc.nrm <- rowSums(df.foc.nrm)  ### ETS at each theta (df cancels across items)
    ETS.ref.nrm <- rowSums(df.ref.nrm)  ### ETS at each theta (df cancels across items)
    ETS.dif.nrm <- ETS.foc.nrm - ETS.ref.nrm   ### DF in ETS at each theta
    ETS.abs.dif.nrm <- abs(ETS.dif.nrm)
    UETSDN <- sum(ETS.abs.dif.nrm*theta.density)
    out.test.stats <- c(STDS,UTDS,UETSDS,ETSSD,Starks.DTFR,UDTFR,UETSDN,test.Dmax)
    out.test.names <- c("STDS","UTDS","UETSDS","ETSSD","Starks.DTFR","UDTFR","UETSDN","theta.of.max.test.D","Test.Dmax")
    df.test.output <- data.frame(out.test.names,out.test.stats)
    names(df.test.output) <- c("Effect Size","Value")
    class(df.item.output) <- c('mirt_df', 'data.frame')
    if(!plot && !DIF) return(df.test.output)

    # plots
    order <- order(Theta.focal[,1])
    grps <- extract.mirt(mod, 'groupNames')
    mykey <- list(space = 'top',
                  columns = 2,
                  text = list(sort(grps)),
                  points = list(pch = 1, col=c("black","red"))
    )
    # fix stupid lattice aesthetic garbage
    if(all(mykey$text[[1]] == grps)){
      mykey$text[[1]] <- c(grps[2], grps[1])
      mykey$points$col <- c("red", "black")
    }
    if(DIF){
        nms <- rep(extract.mirt(mod, 'itemnames')[focal_items], each = nrow(Theta.focal)*2)
        nms <- factor(nms, levels = extract.mirt(mod, 'itemnames')[focal_items])
        plt <- data.frame(S=c(df.ref.obs[order,1],df.foc.obs[order,1]),
                             Theta=c(Theta.focal[order,1], Theta.focal[order,1]),
                             group=c(rep(extract.mirt(mod, 'groupNames'), each = nrow(Theta.focal))))
        if(ncol(df.foc.obs)>1){
          for(i in 2:ncol(df.foc.obs)){
            plt.df <- data.frame(S=c(df.ref.obs[order,i],df.foc.obs[order,i]),
                                 Theta=c(Theta.focal[order,1],Theta.focal[order,1]),
                                 group=rep(extract.mirt(mod, 'groupNames'), each = nrow(Theta.focal)))
            plt <- rbind(plt, plt.df)
          }
        }
        plt$Item <- nms
        return(xyplot(S ~ Theta|Item, plt, type = type, cex = .5,
                      xlab="Focal Group Theta",
                      ylab="Expected Score",
                      groups=plt$group ,
                      col=c("black","red"),
                      key=mykey,
                      main='Expected Scores',
                      par.settings=par.settings,
                      par.strip.text=par.strip.text, ...))
      } else {
          plot.df1 <- data.frame(Theta=Theta.focal[order,1],ETS=ETS.foc.obs[order],group=as.character(grps[2]))
          plot.df2 <- data.frame(Theta=Theta.focal[order,1],ETS=ETS.ref.obs[order],group=as.character(grps[1]))
          plot.df <- rbind(plot.df1,plot.df2)
          main <- if(extract.mirt(mod, 'nitems') == length(focal_items))
              "Expected Test Scores" else "Expected Bundle Scores"
          return(xyplot(ETS~Theta, plot.df, type = type, cex = .6,
                              xlab="Focal Group Theta",
                              ylab="Expected Test Score",
                              groups=plot.df$group ,
                              col=c("black","red"),
                              key = mykey,
                              main=main))
    }# end if plot
    invisible()
}
