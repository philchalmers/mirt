context('bfactor')

test_that('dich data', {
    data <- key2binary(SAT12,
                       key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
    specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
    mod1 <- bfactor(data, specific)    
    expect_is(mod1, 'ConfirmatoryClass')                  
    fs <- fscores(mod1, verbose = FALSE)
    expect_is(fs, 'matrix')        
    cof <- coef(mod1, verbose = FALSE)
    expect_is(cof, 'list')
    sum <- summary(mod1, verbose = FALSE)
    expect_is(sum, 'list')    
    pfit1 <- personfit(mod1)
    expect_is(pfit1, 'data.frame')    
    ifit <- itemfit(mod1)
    expect_is(ifit, 'data.frame')

    #simulate data
    set.seed(1234)
    a <- matrix(c(
        1,0.5,NA,
        1,0.5,NA,
        1,0.5,NA,
        1,0.5,NA,
        1,0.5,NA,
        1,0.5,NA,
        1,0.5,NA,
        1,NA,0.5,
        1,NA,0.5,
        1,NA,0.5,
        1,NA,0.5,
        1,NA,0.5,
        1,NA,0.5,
        1,NA,0.5),ncol=3,byrow=TRUE)
    
    d <- matrix(c(
        -1.0,NA,NA,
        -1.5,NA,NA,
        1.5,NA,NA,
        0.0,NA,NA,
        0.0,-1.0,1.5,
        0.0,2.0,-0.5,
        3.0,2.0,-0.5,
        3.0,2.0,-0.5,
        2.5,1.0,-1,
        2.0,0.0,NA,
        -1.0,NA,NA,
        -1.5,NA,NA,
        1.5,NA,NA,
        0.0,NA,NA),ncol=3,byrow=TRUE)
    
    nominal <- matrix(NA, nrow(d), ncol(d))
    nominal[5, ] <- c(0,1.2,2)    
    sigma <- diag(3)
    set.seed(1234)
    items <- itemtype <- c(rep('dich', 4), 'nominal', 'gpcm', rep('graded',4),rep('dich', 4))
    dataset <- simdata(a,d,2000,itemtype, sigma=sigma, nominal=nominal)  
     
    specific <- c(rep(1,7),rep(2,7))
    items[items == 'dich'] <- '2PL'
    simmod <- bfactor(dataset, specific, itemtype = items)
    expect_is(simmod, 'ConfirmatoryClass')              
    fs <- fscores(simmod, verbose = FALSE)
    expect_is(fs, 'matrix')
    
    res <- residuals(simmod, verbose = FALSE)
    expect_is(res, 'matrix')
    fit <- fitted(simmod)
    expect_is(fit, 'matrix')  
    sum <- summary(simmod, verbose = FALSE)
    expect_is(sum, 'list')
})
