context('mixedmirt')

test_that('mixed dich', {
    set.seed(1234)
    N <- 750
    a <- matrix(rlnorm(10,.2,.5),10,1)
    d <- matrix(rnorm(10), 10)
    Theta <- matrix(sort(rnorm(N)))
    pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
    group <- factor(rep(c('G1','G2','G3'), each = N/3))
    data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
    covdata <- data.frame(group, pseudoIQ)
    mixedmirt1 <- 'Theta = 1-10'
    model <- confmirt.model(mixedmirt1, quiet = TRUE)  
    
    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ group, itemtype = 'Rasch', verbose = FALSE)
    expect_is(mod1, 'MixedClass')                  
    
    #model using 2PL items instead of only Rasch
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ group, verbose = FALSE)    
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1@df - mod1b@df, 10)
    suppressWarnings(dif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, verbose = FALSE))
    expect_is(dif, 'MixedClass')          
    sv <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, pars = 'values')
    constrain <- list(sv$parnum[sv$name == 'groupG2'], sv$parnum[sv$name == 'groupG3']) # main effects
    suppressWarnings(itemdif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, 
                         constrain=constrain, verbose = FALSE, fixed.constrain=FALSE))
    expect_is(itemdif, 'MixedClass')          

    #item covs
    set.seed(1234)
    N <- 750
    a <- matrix(rep(1,10))
    d <- matrix(c(rep(-1,5), rep(1,5)))    
    Theta <- matrix(rnorm(N))
    data <- simdata(a, d, N, itemtype = rep('dich',10), Theta=Theta, D=1)
    itemdesign <- data.frame(itempred=rep(1, ncol(data)))
    mixedmirt1 <- 'Theta = 1-10'
    model <- confmirt.model(mixedmirt1, quiet = TRUE)     
    sv <- mixedmirt(data, model = model, fixed = ~ itempred, pars = 'values', 
                      itemtype = 'Rasch', itemdesign = itemdesign, fixed.constrain=FALSE)
    sv$value[sv$name == 'd'] <- 0
    sv$est[sv$name == 'd'] <- FALSE
    
    #make design such that the first 5 items are systematically more difficult than the last 5
    constrain <- list()
    constrain[[1]] <- sv$parnum[sv$name == 'itempred'][1:5]
    constrain[[2]] <- sv$parnum[sv$name == 'itempred'][-c(1:5)]
    mod <- mixedmirt(data, model = model, fixed = ~ itempred, pars = sv, , fixed.constrain=FALSE,
                  itemtype = 'Rasch', constrain = constrain, itemdesign = itemdesign, verbose = FALSE)   
    expect_is(mod, 'MixedClass')                          
})   

test_that('item and group predictors', {    
    data <- key2binary(SAT12,
                       key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
    model <- confmirt.model('Theta = 1-32', quiet = TRUE)
    
    itemdesign <- data.frame(itemorder = factor(c(rep('easier', 16), rep('harder', 16))))
    fs <- scale(rowSums(data))
    covdata <- data.frame(gender=ifelse(fs > 1, 'M', 'F'))
    
    sv <- mixedmirt(data, covdata, model = model, fixed = ~ itemorder + gender, itemtype = 'Rasch', 
                    itemdesign = itemdesign, pars = 'values')
    expect_is(sv, 'data.frame')       
    LLTM <- mixedmirt(data, covdata, model = model, fixed = ~ itemorder + gender, itemtype = 'Rasch', 
                      itemdesign = itemdesign, verbose = FALSE)
    expect_is(LLTM, 'MixedClass')
    
    sv2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender, itemtype = 'Rasch', 
                     itemdesign = itemdesign, pars='values'))
    expect_is(sv2, 'data.frame')       
    LLTM2 <- suppressWarnings(mixedmirt(data, covdata, model = model, fixed = ~ itemorder * gender, itemtype = 'Rasch', 
                       itemdesign = itemdesign, verbose = FALSE))
    expect_is(LLTM2, 'MixedClass')
    expect_equal(LLTM@df - LLTM2@df, 1)    
    
}) 