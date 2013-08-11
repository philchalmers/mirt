#' Mixed effects modeling for MIRT models
#'
#' \code{mixedmirt} fits MIRT models using FIML estimation to dichotomous and polytomous
#' IRT models conditional on fixed and random effect of person and item level covariates. 
#' This can also be understood as 'expalanatory IRT' if only fixed effects are modeled, or 
#' multilevel/mixed IRT if random and fixed effects are included. The method uses the MH-RM
#' algorithm exclusively. Additionally, computation of the log-likelihood can be sped up by
#' using parallel estimation via \code{\link{mirtCluster}}.
#' 
#' For dichotomous response models, \code{mixedmirt} follows the general form
#' 
#'  \deqn{P(x = 1|\theta, \psi) = g + \frac{(u - g)}{1 + exp(-1 * [\mathbf{\theta a} + 
#'  \mathbf{X \beta} + \mathbf{Z \gamma}])}} 
#'  
#'  where X is a design matrix with associated \eqn{\beta} fixed effect coefficients, and Z is a 
#'  design matrix with associated \eqn{\gamma} random effects. For simplicity and easier 
#'  interpretation, the unique item intercept values typically found in \eqn{\mathbf{X \beta}} 
#'  are extracted and reassigned within mirt's 'intercept' parameters (e.g., \code{'d'}). 
#'  To observe how the design matrices are structured prior to reassignment and estimation pass 
#'  the argument \code{return.design = TRUE}.
#'  
#'  Polytomous IRT models follow a similar format except the item intercepts are automatically estimated
#'  internally, rendering the \code{items} argument in the fixed formula redundent and therefore 
#'  must be ommited from the specification. If there are a mixture of dichotomous and polytomous 
#'  items the intercepts for the dichotomous models are also estimated for consistency.
#'  
#'  To simulate maximum a posteriori estimates for the random effects use the \code{\link{randef}}
#'  function.
#'
#' @aliases mixedmirt coef,MixedClass-method summary,MixedClass-method anova,MixedClass-method
#' @param data a \code{matrix} or \code{data.frame} that consists of
#' numerically ordered data, with missing data coded as \code{NA}
#' @param covdata a \code{data.frame} that consists of the \code{nrow(data)} by \code{K}
#' 'person level' fixed and random predictors
#' @param model an object returned from \code{mirt.model()} declaring how
#' the factor model is to be estimated. See \code{\link{mirt.model}} for
#' more details
#' @param fixed a right sided R formula for specifying the fixed effect (aka 'explanatory') 
#' predictors from \code{covdata} and \code{itemdesign}. To estimate the intercepts for 
#' each item the keyword \code{items} is reserved and automatically added to the \code{itemdesign} input.
#' If any polytomous items are being model the \code{items} are argument is not valid since all
#' intercept parameters are freely estimated and identified with the parameterizations found in 
#' \code{\link{mirt}}, and the first column in the fixed design matrix (commonly the intercept or a reference
#' group) is ommited
#' @param random a right sided formula or list of formulas containing crossed random effects 
#' of the form \code{v1 + ... v_n | G}, where \code{G} is the grouping variable and \code{v_n} are 
#' random numeric predictors within each group. If no intercept value is specified then by default the 
#' correlations between the \code{v}'s and \code{G} are estimated, but can be supressed by including 
#' the \code{~ -1 + ...} constant 
#' @param itemtype same as itemtype in \code{\link{mirt}}, expect does not support the following 
#' item types: \code{c('PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')}
#' @param itemdesign a \code{data.frame} object used to create a design matrix for the items, where 
#' each \code{nrow(itemdesign) == nitems} and the number of columns is equal to the number of fixed 
#' effect predictors (i.e., item intercepts). By default an \code{items} variable is reserved for 
#' modeling the item intercept parameters
#' @param constrain a list indicating parameter equality constrains. See \code{\link{mirt}} for 
#' more detail
#' @param pars used for parameter starting values. See \code{\link{mirt}} for more detail
#' @param return.design logical; return the design matrices before they have (potentially) 
#' been reassigned? 
#' @param draws the number of Monte Carlo draws to estimate the log-likelihood for the MH-RM algorithm. Default
#' is 5000
#' @param ... additional arguments to be passed to the MH-RM estimation engine. See 
#' \code{\link{mirt}} for more detail
#'
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @seealso \code{\link{randef}}, \code{\link{calcLogLik}}, \code{\link{mirtCluster}}
#' @export mixedmirt
#' @examples
#'
#' \dontrun{
#' 
#' #make some data
#' set.seed(1234)
#' N <- 750
#' a <- matrix(rlnorm(10,.3,1),10,1)
#' d <- matrix(rnorm(10), 10)
#' Theta <- matrix(sort(rnorm(N)))
#' pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
#' group <- factor(rep(c('G1','G2','G3'), each = N/3))
#' data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
#' covdata <- data.frame(group, pseudoIQ)
#' #use parallel computing
#' mirtCluster()
#'
#' #specify IRT model
#' model <- mirt.model('Theta = 1-10')
#'
#' #model with no person predictors
#' mod0 <- mirt(data, model, itemtype = 'Rasch')
#'
#' #group as a fixed effect predictor (aka, uniform dif) 
#' mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items)
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
#' mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + items, itemtype = '2PL')
#' anova(mod1, mod1b) #much better with 2PL models using all criteria (as expected, given simdata pars)
#'
#' #continuous predictor and interaction model with group 
#' mod2 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group + pseudoIQ)
#' mod3 <- mixedmirt(data, covdata, model, fixed = ~ 0 + group * pseudoIQ)
#' summary(mod2)
#' anova(mod1b, mod2)
#' anova(mod2, mod3)
#' 
#' #view fixed design matrix with and without unique item level intercepts
#' withint <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, return.design = TRUE)
#' withoutint <- mixedmirt(data, covdata, model, fixed = ~ 0 + group, return.design = TRUE)
#' 
#' #notice that in result above, the intercepts 'items1 to items 10' were reassigned to 'd'
#' head(withint$X) 
#' tail(withint$X) 
#' head(withoutint$X) #no intercepts design here to be reassigned into item intercepts
#' tail(withoutint$X)  
#' 
#' ###################################################
#' ### random effects
#' #make the number of groups much larger
#' covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
#' 
#' #random groups
#' rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group)
#' summary(rmod1)
#' coef(rmod1)
#' 
#' #random groups and random items 
#' rmod2 <- mixedmirt(data, covdata, 1, random = list(~ 1|group, ~ 1|items))
#' summary(rmod2)
#' eff <- randef(rmod2) #estimate random effects
#' 
#' #random slopes with fixed intercepts (suppressed correlation)
#' rmod3 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ -1 + pseudoIQ|group)
#' summary(rmod3)
#' (eff <- randef(rmod3)) 
#' 
#' ###################################################
#' ##LLTM, and 2PL version of LLTM
#' data(SAT12)
#' data <- key2binary(SAT12,
#'                    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' model <- mirt.model('Theta = 1-32')
#'
# #Suppose that the first 16 items were suspected to be easier than the last 16 items, and we wish
# #to test this item structure hypothesis (more intercept designs are possible by including more columns).
#' itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))
#'
#' #notice that the 'fixed = ~ ... + items' argument is ommited
#' LLTM <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemdesign = itemdesign)
#' summary(LLTM)
#' coef(LLTM)
#' wald(LLTM)
#' L <- c(-1, 1, 0)
#' wald(LLTM, L) #first half different from second
#'
#' #compare to items with estimated slopes (2PL)
#' twoPL <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, itemtype = '2PL', 
#'                    itemdesign = itemdesign)
#' anova(twoPL, LLTM) #much better fit
#' summary(twoPL)
#' coef(twoPL)
#' 
#' wald(twoPL)
#' L <- matrix(0, 1, 34)
#' L[1, 1] <- 1
#' L[1, 2] <- -1
#' wald(twoPL, L) #n.s. 
#' #twoPL model better than LLTM, and don't draw the incorrect conclusion that the first
#' #    half of the test is any easier/harder than the last
#' 
#' ##LLTM with item error term 
#' LLTMwithError <- mixedmirt(data, model = model, fixed = ~ 0 + itemorder, random = ~ 1|items,
#'     itemdesign = itemdesign)
#' summary(LLTMwithError) 
#' #large item level variance after itemorder is regressed; not a great predictor of item difficulty
#' coef(LLTMwithError) 
#' 
#' ###################################################
#' ### Polytomous example
#' 
#' #make an arbitrary group difference
#' covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
#' 
#' #partial credit model
#' mod <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group)
#' coef(mod)
#' 
#' #gpcm to estimate slopes 
#' mod2 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
#'                  itemtype = 'gpcm')
#' summary(mod2)
#' anova(mod, mod2)
#'
#' #graded model
#' mod3 <- mixedmirt(Science, covdat, model=1, fixed = ~ 0 + group,
#'                  itemtype = 'graded')
#' coef(mod3)
#' 
#' }
mixedmirt <- function(data, covdata = NULL, model, fixed = ~ 1, random = NULL, itemtype = 'Rasch',
                      itemdesign = NULL, constrain = NULL, pars = NULL, return.design = FALSE, 
                      draws = 5000, ...)
{
    Call <- match.call()       
    svinput <- pars
    if(length(itemtype) == 1L) itemtype <- rep(itemtype, ncol(data))    
    if(any(itemtype %in% c('PC2PL', 'PC3PL', '2PLNRM', '3PLNRM', '3PLuNRM', '4PLNRM')))
        stop('itemtype contains unsupported classes of items')    
    if(is(random, 'formula')) {
        random <- list(random)
    } else if(is.null(random)) random <- list()
    RETVALUES <- ifelse(is.character(pars), TRUE, FALSE) 
    if(!is.list(random)) stop('Incorrect input for random argument')
    if(is.null(covdata)) covdata <- data.frame(UsElEsSvAR = factor(rep(1L, nrow(data))))
    if(is.null(itemdesign)){
        itemdesign <- data.frame(items = factor(1L:ncol(data)))
    } else itemdesign$items <- factor(1L:ncol(data))    
    if(!is.data.frame(covdata) || ! is.data.frame(itemdesign))
        stop('Predictor variable inputs must be data.frame objects')    
    dropcases <- which(rowSums(is.na(covdata)) != 0)
    if(length(dropcases) > 0L){
        data <- data[-dropcases, ]
        covdata <- covdata[-dropcases, ]
    }
    longdata <- reshape(data.frame(ID=1L:nrow(data), data, covdata), idvar='ID', 
                        varying=list(1L:ncol(data) + 1L), direction='long')
    colnames(longdata) <- c('ID', colnames(covdata), 'items', 'response')
    for(i in 1L:ncol(itemdesign))
        longdata[, colnames(itemdesign)[i]] <- rep(itemdesign[ ,i], each=nrow(data))
    mf <- model.frame(fixed, longdata)
    mm <- model.matrix(fixed, mf)   
    K <- sapply(as.data.frame(data), function(x) length(na.omit(unique(x))))
    if(any(K > 2)){
        if(any(colnames(mm) %in% paste0('items', 1:ncol(data))))
            stop('fixed formulas do no support the \'items\' internal variable for 
                 polytomous items. Please remove')
        mm <- mm[ , -1L, drop = FALSE]
    }
    if(return.design) return(list(X=mm, Z=NaN))
    itemindex <- colnames(mm) %in% paste0('items', 1L:ncol(data))
    mmitems <- mm[, itemindex]
    mm <- mm[ ,!itemindex, drop = FALSE]
    if(length(random) > 0L){
        mr <- make.randomdesign(random=random, longdata=longdata, covnames=colnames(covdata), 
                                itemdesign=itemdesign, N=nrow(covdata))
    } else mr <- list()       
    mixed.design <- list(fixed=mm, random=mr)    
    if(is.null(constrain)) constrain <- list()      
    sv <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype,
                     D=1, mixed.design=mixed.design, method='MIXED', constrain=NULL, pars='values')
    mmnames <- colnames(mm)
    N <- nrow(data)
    if(ncol(mm) > 0L){
        for(i in 1L:ncol(mm)){
            mmparnum <- sv$parnum[sv$name == mmnames[i]]            
            constrain[[length(constrain) + 1L]] <- mmparnum
        }            
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
    } else pars <- sv
    if(RETVALUES){
        attr(pars, 'values') <- NULL
        return(pars)
    }
    if(is.data.frame(svinput)) pars <- svinput
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype,
                      mixed.design=mixed.design, method='MIXED', constrain=constrain, pars=pars,
                      draws=draws, ...)
    if(is(mod, 'MixedClass'))
        mod@Call <- Call
    return(mod)
}
