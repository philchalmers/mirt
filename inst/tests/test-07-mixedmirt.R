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
    
    model <- confmirt.model('confmods/mixedmirt1', quiet = TRUE)  
    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ group, itemtype = 'Rasch', verbose = FALSE)
    expect_is(mod1, 'MixedClass')                  
    #model using 2PL items instead of only Rasch
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ group, verbose = FALSE)    
    expect_is(mod1b, 'MixedClass')                      
    dif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, verbose = FALSE)
    expect_is(dif, 'MixedClass')          
    sv <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, pars = 'values')
    constrain <- list(sv$parnum[sv$name == 'groupG2'], sv$parnum[sv$name == 'groupG3']) # main effects
    itemdif <- mixedmirt(data, covdata, model, fixed = ~ group * Theta, fixed.constrain = FALSE,
                         constrain=constrain, verbose = FALSE)
    expect_is(itemdif, 'MixedClass')          
})

test_that('mixed poly', {
    set.seed(1234)
    N <- 750
    a <- matrix(rlnorm(10,.2,.5),10,1)
    d <- matrix(rnorm(10), 10)
    Theta <- matrix(sort(rnorm(N)))
    pseudoIQ <- Theta * 5 + 100  + rnorm(N, 0 , 5)
    group <- factor(rep(c('G1','G2','G3'), each = N/3))
    data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
    covdata <- data.frame(group, pseudoIQ)
    
})

test_that('item covs', {
    set.seed(1234)
    N <- 750
    a <- matrix(rep(1,10))
    d <- matrix(c(rep(-1,5), rep(1,5)))    
    Theta <- matrix(rnorm(N))
    data <- simdata(a, d, N, itemtype = rep('dich',10), Theta=Theta, D=1)
    covdata <- data.frame(itempred=rep(1, nrow(data)))
    # d > 0
    model <- confmirt.model('confmods/mixedmirt1', quiet = TRUE)  
    sv <- mixedmirt(data, covdata, model, fixed= ~ itempred, pars = 'values', 
                               itemtype = 'Rasch')
    sv$value[sv$name == 'd'] <- 0
    sv$est[sv$name == 'd'] <- FALSE
    constrain <- list()
    constrain[[1]] <- sv$parnum[sv$name == 'itempred'][1:5]
    constrain[[2]] <- sv$parnum[sv$name == 'itempred'][-c(1:5)]
    mod <- mixedmirt(data, covdata, model, fixed= ~ itempred, pars = sv, 
                     itemtype = 'Rasch', constrain = constrain, fixed.constrain = FALSE, 
                     verbose = FALSE)    
    expect_is(mod, 'MixedClass')                          
})    