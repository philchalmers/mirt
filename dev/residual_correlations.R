#' ### load SAT12 and compute bifactor model with 3 specific factors
#' data(SAT12)
#' data <- key2binary(SAT12,
#'   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
#' specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
#'
#' # alternative model definition via ?mirt.model syntax
#' specific2 <- "S1 = 7,9,10,11,13,15,17,18,21,22,24,27,31
#'               S2 = 1,3,6,8,16,29,32
#'               S3 = 2,4,5,12,14,19,20,23,25,26,28,30"
#' mod2 <- bfactor(data, specific2)
#'
#' #########
#' # Defining unconstrained residual correlations via item doublets
#' # Requires one specific factor per bivariate residual correlation with
#' # slopes constrained to be equal. NOTE: cannot be performed with bfactor(),
#' # but is included for comparison
#'
#' # SAT12 where second specific factor has unconstrained residual correlation
#' # structure (more complex and less causally specific, but more flexible)
#'
#' specific_model <- "G = 1-32
#'                    S1 = 7,9,10,11,13,15,17,18,21,22,24,27,31
#'                    S3 = 2,4,5,12,14,19,20,23,25,26,28,30"
#' ucor <- c(1,3,6,8,16,29,32)
#'
#' # there are this many bivariate correlations to define
#' r <- length(ucor)
#' r * (r-1)/2
#'
#' constrain <- NULL
#' anum <- 4    # in this model, residual factors start at slope a4
#' for(i in 1:length(ucor)){
#'    for(j in 1:length(ucor)){
#'       if(i < j){
#'           specific_model <- c(specific_model,
#'                           paste0('s', ucor[i], '_', ucor[j], ' = ', ucor[i], ',', ucor[j]))
#'           constrain <- c(constrain, sprintf("(%i,%i,a%i)", ucor[i], ucor[j], anum))
#'           anum <- anum + 1
#'       }
#'    }
#' }
#'
#' constrain <- paste0('\nCONSTRAIN = ', paste0(constrain, collapse=', '))
#' specific_model <- paste0(paste0(specific_model, collapse='\n'),
#'                                 constrain, collapse='')
#' cat(specific_model)
#'
#' # many more parameters to estimate than bifactor, and numerical
#' # integration more difficult (Monte Carlo EM used here with 500 draws)
#' resid_mod <- mirt(data, specific_model, method = 'MCEM')
#' coef(resid_mod, simplify=TRUE)$items
#' anova(mod2, resid_mod)