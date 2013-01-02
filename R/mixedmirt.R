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
#' @param fixed a standard R formula for specifying the fixed effect predictors from \code{covdata} and 
#' \code{itemdesign}. By default constraints are not imposed, so the fixed person effects are 
#' not equal accross items, but this can be enabled using \code{fixed.constrain = TRUE} 
#' @param random a formula similar to the \code{nlme} random variable sepcifications for declaring
#' the random slope and intercept predictors. Not currently available, but will be soon (hopefully)....
#' @param itemtype same as itemtype in \code{\link{mirt}}
#' @param itemdesign a data.frame object used to create a design matrix for the items, where each 
#' \code{nrow(itemdesign) == nitems} and the number of columns is equal to the number of fixed effect 
#' predictors (i.e., item intercepts). For design based effects use a matrix of 1's, and apply the appropriate
#' item constraints (see examples below) 
#' @param fixed.constrain logical; constrain the fixed person effects to be equal accross items? Disable
#' this when modelling item level covariates and apply the constraints manually
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
#' N <- 750
#' a <- matrix(rlnorm(10,.2,.5),10,1)
#' d <- matrix(rnorm(10), 10)
#' Theta <- matrix(sort(rnorm(N)))
#' pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
#' group <- factor(rep(c('G1','G2','G3'), each = N/3))
#' data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
#' covdata <- data.frame(group, pseudoIQ) 
#' 
#' #specify IRT model 
#' model <- confmirt.model()
#'      Theta = 1-10
#'      
#'      
#' #model with no person predictors
#' mod0 <- mirt(data, model, itemtype = 'Rasch')
#' #group as a fixed effect predictor (aka, uniform dif) and equal effect for all items with 
#' #   fixed.constrain = TRUE
#' mod1 <- mixedmirt(data, covdata, model, fixed = ~ group, itemtype = 'Rasch', fixed.constrain = TRUE)   
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
#' mod1b <- mixedmirt(data, covdata, model, fixed = ~ group, itemtype = '2PL', fixed.constrain = TRUE)
#' anova(mod1, mod1b) #much better with 2PL models using all criteria (as expected, given simdata pars)
#' 
#' #global nonuniform dif (group interaction with latent variable)
#' (dif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, fixed.constrain = TRUE))
#' anova(mod1b, dif)
#' #free the interaction terms for each item to detect dif (less power for each item though)
#' sv <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, pars = 'values')
#' constrain <- list(sv$parnum[sv$name == 'groupG2'], sv$parnum[sv$name == 'groupG3']) # main effects
#' itemdif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, fixed.constrain = FALSE, 
#'      constrain=constrain)
#' anova(dif, itemdif)
#'                                       
#' #continuous predictor and interaction model with group 
#' mod2 <- mixedmirt(data, covdata, model, fixed = ~ group + pseudoIQ, fixed.constrain = TRUE)
#' mod3 <- mixedmirt(data, covdata, model, fixed = ~ group * pseudoIQ, fixed.constrain = TRUE)
#' summary(mod2)
#' anova(mod1b, mod2)
#' anova(mod2, mod3)  
#' 
#' ###############################################################
#' ###### Advanced use of mixedmirt: LLTM, and 2PL version of LLTM
#' ###############################################################
#'  
#' # flexible LLTM model by customizing the structure
#' data(SAT12)
#' data <- key2binary(SAT12,
#'                    key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' model <- confmirt.model()
#' Theta = 1-32
#' 
#' 
# #Suppose that the first 16 items were suspected to be easier than the last 16 items, and we wish
# #to test this item structure hypothesis. First, create a design matrix entirely of 1's for each
# #item design effect (only one here)
#' itemdesign <- data.frame(itemorder = matrix(1, 32))
#' 
# #Obtain starting values and identify parameters that are to be equal (first 16, last 16), and 
# #fix item intercepts to 0
#' sv <- mixedmirt(data, model = model, fixed = ~ itemorder, itemtype = 'Rasch', 
#'                 itemdesign = itemdesign, pars = 'values')
#' constrain <- list()
#' constrain[[1]] <- sv$parnum[sv$name == 'itemorder'][1:16]
#' constrain[[2]] <- sv$parnum[sv$name == 'itemorder'][-(1:16)]
#' sv$value[sv$name == 'd'] <- 0
#' sv$est[sv$name == 'd'] <- FALSE
#' #estimate the LLTM model using the new starting values and constraints
#' LLTM <- mixedmirt(data, model = model, fixed = ~ itemorder, itemtype = 'Rasch', 
#'         itemdesign = itemdesign, pars = sv, constrain = constrain)
#' coef(LLTM)
#' #compare to standard items with estimated slopes (2PL)?
#' sv$est[sv$name == 'a1'] <- TRUE
#' twoPL <- mixedmirt(data, model = model, fixed = ~ itemorder, 
#'         itemdesign = itemdesign, pars = sv, constrain = constrain)
#' coef(twoPL)
#' anova(twoPL, LLTM) 
#' #twoPL model better than LLTM, and don't draw the (spurious?) conclusion that the first 
#' #    half of the test is any easier/harder than the last
#' 
#' ### Similar example, but with simulated data instead
#' 
#' set.seed(1234)
#' N <- 750
#' a <- matrix(rep(1,10))
#' d <- matrix(c(rep(-1,5), rep(1,5)))    
#' Theta <- matrix(rnorm(N))
#' data <- simdata(a, d, N, itemtype = rep('dich',10), Theta=Theta, D=1)
#' itemdesign <- data.frame(itempred=rep(1, ncol(data)))
#' model <- confmirt.model()
#'    Theta = 1-10
#' 
#' 
#' sv <- mixedmirt(data, model = model, fixed = ~ itempred, pars = 'values', 
#'                  itemtype = 'Rasch', itemdesign = itemdesign)
#' sv$value[sv$name == 'd'] <- 0
#' sv$est[sv$name == 'd'] <- FALSE
#' 
#' #make design such that the first 5 items are systematically more difficult than the last 5
#' constrain <- list()
#' constrain[[1]] <- sv$parnum[sv$name == 'itempred'][1:5]
#' constrain[[2]] <- sv$parnum[sv$name == 'itempred'][-c(1:5)]
#' mod <- mixedmirt(data, model = model, fixed = ~ itempred, pars = sv, 
#'                  itemtype = 'Rasch', constrain = constrain, itemdesign = itemdesign)
#' coef(mod)                 
#' rasch <- mirt(data, 1, itemtype = 'Rasch', D=1)
#' anova(mod, rasch) #n.s., LLTM model a much better choice compared to Rasch
#' 
#' }
mixedmirt <- function(data, covdata = NULL, model, fixed = ~ 1, random = NULL, itemtype = NULL, 
                      itemdesign = NULL, fixed.constrain = FALSE, ...)
{       
    Call <- match.call() 
    if(is.null(covdata))
        covdata <- data.frame(InTeRnAlUsElESsNaMe = matrix(1, nrow(data)))
    if(is.null(itemdesign))
        itemdesign <- data.frame(InTeRnAlUsElESsNaMe2 = matrix(1, ncol(data)))
    for(i in 1:ncol(covdata))
        if(is(covdata[,i], 'numeric') || is(covdata[,i], 'integer'))            
            covdata[,i] <- matrix(scale(covdata[,i], scale = FALSE))    
    ### TEMPORARY    
    if(!is.null(random)) 
        stop('random effect covariates not yet supported')
    random <- ~ 1  
    ###
    if(fixed == ~ 1 && random == ~ 1)
        stop('No fixed or random effects have been specified.')
    Theta <- matrix(0, nrow(data), nrow(model$x))
    colnames(Theta) <- model$x[,1]
    fixed.identical <- FALSE
    if(all(itemdesign == 1)) fixed.identical <- TRUE
    fixed.design.list <- designMats(covdata=covdata, fixed=fixed, Thetas=Theta, nitems=ncol(data),
                                   itemdesign=itemdesign, fixed.identical=fixed.identical)    
    mixedlist <- list(fixed=fixed, random=random, covdata=covdata, factorNames=model$x[,1], 
                      FDL=fixed.design.list, itemdesign=itemdesign, fixed.constrain=fixed.constrain,
                      fixed.identical=fixed.identical)    
    mod <- ESTIMATION(data=data, model=model, group=rep('all', nrow(data)), itemtype=itemtype, 
                      D=1, mixedlist=mixedlist, method='MIXED', ...)
    if(is(mod, 'MixedClass'))
        mod@Call <- Call
    return(mod)    
}