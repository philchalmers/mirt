#' Mixed effects modeling for MIRT models
#' 
#' \code{mixedmirt} fits MIRT models using FIML estimation to dichotomous and polytomous
#' IRT models conditional on fixed and random effect covariates. The method uses the MH-RM 
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
#' @param random a formula similar to the \code{nlme} random variable sepcifications for declaring
#' the random slope and intercept predictors. Not currently available, but will be soon (hopefully)....
#' @param itemtype same as itemtype in \code{\link{mirt}}
#' @param ... additinonal arguments to be passed to the MH-RM estimation engine. See \code{\link{confmirt}}
#' for more detail
#' 
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @export mixedmirt
#' @examples
#'  
#' \dontrun{
#' #make some data
#' set.seed(1234)
#' a <- matrix(rlnorm(10,.2,.5),10,1)
#' d <- matrix(rnorm(10), 10)
#' Theta <- matrix(sort(rnorm(750)))
#' pseudoIQ <- Theta * 5 + 100  + rnorm(750, 0 , 5)
#' group <- factor(rep(c('G1','G2','G3'), each = 250))
#' data <- simdata(a,d,750, itemtype = rep('dich',10), Theta=Theta)
#' covdata <- data.frame(group, pseudoIQ) 
#' 
#' #specify IRT model 
#' model <- confmirt.model()
#'      Theta = 1-10
#'      
#'      
#' #model with no person predictors
#' mod0 <- mirt(data, model, itemtype = 'Rasch')
#' #group as a fixed effect predictor     
#' mod1 <- mixedmirt(data, covdata, model, fixed = ~ group, itemtype = 'Rasch')   
#' anova(mod0, mod1)
#' summary(mod1)
#' coef(mod1) 
#' 
#' #same model as above in lme4
#' wide <- data.frame(id=1:750,data,covdata)      
#' long <- reshape2::melt(wide, id.vars = c('id', 'group', 'pseudoIQ'))
#' library(lme4)
#' lmod0 <- lmer(value ~ 0 + variable + (1|id), long, family = binomial)
#' lmod1 <- lmer(value ~ 0 + group + variable + (1|id), long, family = binomial)
#' anova(lmod0, lmod1)
#' 
#' #model using 2PL items instead of only Rasch
#' mod1b <- mixedmirt(data, covdata, model, fixed = ~ group)
#' anova(mod1, mod1b) #much better with 2PL models using all criteria (as expected, given simdata pars)
#'
#' #continuous predictor and interaction model with group 
#' mod2 <- mixedmirt(data, covdata, model, fixed = ~ group + pseudoIQ)
#' mod3 <- mixedmirt(data, covdata, model, fixed = ~ group * pseudoIQ)
#' summary(mod2)
#' anova(mod1b, mod2)
#' anova(mod2, mod3)  
#'      
#' }
mixedmirt <- function(data, covdata, model, fixed = ~ 1, random = NULL, itemtype = NULL, ...)
{       
    Call <- match.call() 
    for(i in 1:ncol(covdata))
        if(is(covdata[,i], 'numeric') || is(covdata[,i], 'integer'))
            covdata[,i] <- matrix(scale(covdata[,i], scale = FALSE))    
    if(fixed == ~ 1) {
        fixed.design <- NULL
    } else fixed.design <- model.matrix(fixed, covdata)[ ,-1, drop = FALSE]
    ### TEMPORARY    
    if(!is.null(random)) 
        stop('random effect covariates not yet supported')
    random <- ~ 1 
    if(random == ~ 1){
        random.design <- NULL
    } else random.design <- model.matrix(random, covdata)[ ,-1, drop = FALSE]    
    ###
    if(is.null(fixed.design) && is.null(random.design))
        stop('No fixed or random effects have been specified.')
        
    mixedlist <- list(fixed.design=fixed.design, betas=rep(0, ncol(fixed.design)), 
                      random.design=random.design)    
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype, 
                      D=1, mixedlist=mixedlist, method='MIXED', ...)
    if(is(mod, 'MixedClass'))
        mod@Call <- Call
    return(mod)    
}