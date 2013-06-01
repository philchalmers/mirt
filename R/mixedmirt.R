#' Mixed effects modeling for MIRT models
#'
#' \code{mixedmirt} fits MIRT models using FIML estimation to dichotomous and polytomous
#' IRT models conditional on fixed and random effect of person and item level covariates. 
#' The method uses the MH-RM
#' algorithm exclusively. The D scaling parameter is automatically fixed to 1 so that all
#' coefficients can be interpreted on the exponential metric.
#'
#' @aliases mixedmirt coef,MixedClass-method summary,MixedClass-method anova,MixedClass-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param covdata a \code{data.frame} that consists of the \code{nrow(data)} by \code{K}
#' 'person level' fixed and random predictors
#' @param model an object returned from \code{confmirt.model()} declaring how
#' the factor model is to be estimated. See \code{\link{confmirt.model}} for
#' more details
#' @param fixed a standard R formula for specifying the fixed effect predictors from \code{covdata} 
#' and \code{itemdesign}. To estimate the intercepts for each item the keyword \code{items} is 
#' reserved and automatically added to the \code{itemdesign} input
#' @param random a formula similar to the \code{nlme} random variable specifications for declaring
#' the random slope and intercept predictors. Not currently available, but will be available some 
#' time in the future 
#' @param itemtype same as itemtype in \code{\link{mirt}}, expect limited only to the following 
#' item types: \code{c('Rasch', '2PL', '3PL', '3PLu', '4PL', 'gpcm', 'rsm')}
#' @param itemdesign a \code{data.frame} object used to create a design matrix for the items, where 
#' each \code{nrow(itemdesign) == nitems} and the number of columns is equal to the number of fixed 
#' effect predictors (i.e., item intercepts)
#' @param constrain a list indicating parameter equality constrains. See \code{\link{mirt}} for 
#' more detail
#' @param pars used for parameter starting values. See \code{\link{mirt}} for more detail
#' @param ... additinonal arguments to be passed to the MH-RM estimation engine. See 
#' \code{\link{confmirt}} for more detail
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export mixedmirt
#' @examples
#'
#' \dontrun{
#' 
#' #make some data
#' set.seed(1234)
#' N <- 750
#' a <- matrix(rlnorm(10,.2,.5),10,1)
#' d <- matrix(rnorm(10), 10)
#' Theta <- matrix(sort(rnorm(N)))
#' pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
#' group <- factor(rep(c('G1','G2','G3'), each = N/3))
#' data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
#' covdata <- data.frame(group, pseudoIQ)
#' #use cl for parallel computing
#' library(parallel)
#' cl <- makeCluster(detectCores())
#'
#' #specify IRT model
#' model <- confmirt.model('Theta = 1-10')
#'
#' #model with no person predictors
#' mod0 <- mirt(data, model, itemtype = 'Rasch', D = 1)
#'
#' #group as a fixed effect predictor (aka, uniform dif) 
#' mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items, cl=cl)
#' anova(mod0, mod1)
#' summary(mod1)
#' coef(mod1)
#'
#' #same model as above in lme4
#' wide <- data.frame(id=1:nrow(data),data,covdata)
#' long <- reshape2::melt(wide, id.vars = c('id', 'group', 'pseudoIQ'))
#' library(lme4)
#' lmod0 <- lmer(value ~ 0 + variable + (1|id), long, family = binomial)
#' lmod1 <- lmer(value ~ 0 + group + variable + (1|id), long, family = binomial)
#' anova(lmod0, lmod1)
#'
#' #model using 2PL items instead of Rasch
#' mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items, itemtype = '2PL', cl=cl)
#' anova(mod1, mod1b) #much better with 2PL models using all criteria (as expected, given simdata pars)
#'
#' #continuous predictor and interaction model with group in Rasch model
#' mod2 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group + pseudoIQ, cl=cl)
#' mod3 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group * pseudoIQ, cl=cl)
#' summary(mod2)
#' anova(mod1b, mod2)
#' anova(mod2, mod3)
#'
#' ###################################################
#' ##LLTM, and 2PL version of LLTM
#' data(SAT12)
#' data <- key2binary(SAT12,
#'                    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' model <- confmirt.model('Theta = 1-32')
#'
# #Suppose that the first 16 items were suspected to be easier than the last 16 items, and we wish
# #to test this item structure hypothesis (more intercept designs are possible by including more columns).
#' itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))
#'
#' #notice that the 'fixed = ~ ... + items' argument is ommited
#' LLTM <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemdesign = itemdesign, cl=cl)
#' coef(LLTM)
#' wald(LLTM)
#' L <- matrix(c(-1, 1), 1)
#' wald(LLTM, L) #first half different from second
#'
#' #compare to items with estimated slopes (2PL)
#' twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL', 
#'                    itemdesign = itemdesign, cl=cl)
#' coef(twoPL)
#' wald(twoPL)
#' L <- matrix(0, 1, 34)
#' L[1, 1] <- 1
#' L[1, 2] <- -1
#' wald(twoPL, L) #n.s.
#' anova(twoPL, LLTM)
#' #twoPL model better than LLTM, and don't draw the incorrect conclusion that the first
#' #    half of the test is any easier/harder than the last
#' 
#' ### Polytomous example
#' #make an arbitrary group difference
#' covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
#' 
#' mod <- mixedmirt(Science, covdat, model=confmirt.model('F1 = 1-4'), fixed = ~ 0 + group + items)
#' coef(mod)
#' 
#' #gpcm to estimate slopes 
#' mod2 <- mixedmirt(Science, covdat, model=confmirt.model('F1 = 1-4'), fixed = ~ 0 + group + items),
#'                  itemtype = 'gpcm')
#' coef(mod2)
#' anova(mod, mod2)
#' 
#' }
mixedmirt <- function(data, covdata = NULL, model, fixed = ~ 1, random = NULL, itemtype = 'Rasch',
                      itemdesign = NULL, constrain = NULL, pars = NULL, ...)
{
    Call <- match.call()       
    if(length(itemtype) == 1) itemtype <- rep(itemtype, ncol(data))
    if(!all(itemtype %in% c('Rasch', '2PL', '3PL', '3PLu', '4PL', 'gpcm', 'rsm')))
        stop('itemtype contains unsupported classes of items')
    if(!is.null(random)) stop('random effects not yet supported')    
    if(is.null(covdata)) covdata <- data.frame(UsElEsSvAR = factor(rep(1L, nrow(data))))
    if(is.null(itemdesign)){
        itemdesign <- data.frame(items = factor(1L:ncol(data)))
    } else itemdesign$items <- factor(1L:ncol(data))
    if(!is.null(random)) stop('random effects not yet supported')
    if(!is.data.frame(covdata) || ! is.data.frame(itemdesign))
        stop('Predictor variable inputs must be data.frame objects')
    longdata <- reshape(data.frame(ID=1L:nrow(data), data, covdata), idvar='ID', 
                        varying=list(1L:ncol(data) + 1L), direction='long')
    colnames(longdata) <- c('ID', colnames(covdata), 'items', 'response')
    for(i in 1L:ncol(itemdesign))
        longdata[, colnames(itemdesign)[i]] <- rep(itemdesign[ ,i], each=nrow(data))
    mf <- model.frame(fixed, longdata)
    mm <- model.matrix(fixed, mf)   
    itemindex <- colnames(mm) %in% paste0('items', 1L:ncol(data))
    mmitems <- mm[ , itemindex]
    mm <- mm[ ,!itemindex, drop = FALSE]
    mixed.design <- list(fixed=mm, random=NaN)
    if(is.null(constrain)) constrain <- list()    
    sv <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype,
                     D=1, mixed.design=mixed.design, method='MIXED', constrain=NULL, pars='values')
    mmnames <- colnames(mm)
    N <- nrow(data)
    for(i in 1L:ncol(mm)){
        mmparnum <- sv$parnum[sv$name == mmnames[i]]            
        constrain[[length(constrain) + 1L]] <- mmparnum
    }        
    if(ncol(mmitems) > 0L){                
        itemindex <- colnames(data)[which(paste0('items', 1L:ncol(data)) %in% colnames(mmitems))]        
        for(i in itemindex){
            tmp <- sv[i == sv$item, ] 
            tmp$est[tmp$name %in% paste0('d', 1L:50L)] <- TRUE
            tmp$est[tmp$name == 'd'] <- TRUE
            sv[i == sv$item, ] <- tmp
        }
        attr(sv, 'values') <- pars
        pars <- sv
    }
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype,
                      D=1, mixed.design=mixed.design, method='MIXED', constrain=constrain, pars=pars, ...)
    if(is(mod, 'MixedClass'))
        mod@Call <- Call
    return(mod)
}
