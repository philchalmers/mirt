context('confmirt')

test_that('all mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    explor <- confmirt(fulldata, 1, verbose = FALSE)
    expect_is(explor, 'confmirtClass')
    
    set.seed(1234)
    a <- matrix(c(
        1.5,NA,
        0.5,NA,
        1.0,NA,
        1.0,0.5,
        NA,1.5,
        NA,0.5,
        NA,1.0,
        NA,1.0),ncol=2,byrow=TRUE)
    
    d <- matrix(c(
        -1.0,NA,NA,
        -1.5,NA,NA,
        1.5,NA,NA,
        0.0,NA,NA,
        3.0,2.0,-0.5,
        2.5,1.0,-1,
        2.0,0.0,NA,
        1.0,NA,NA),ncol=3,byrow=TRUE)
    
    sigma <- diag(2)
    sigma[1,2] <- sigma[2,1] <- .4
    items <- c(rep('dich',4), rep('graded',3), 'dich')
    dataset <- simdata(a,d,2000,items,sigma)
    
    #analyses
    #CIFA for 2 factor crossed structure    
        
    model.1 <- confmirt.model('confmods/model1', quiet = TRUE)    
    model.quad <- confmirt.model('confmods/modelquad', quiet = TRUE)
    model.combo <- confmirt.model('confmods/modelcombo', quiet = TRUE)    
    
    mod1 <- confmirt(dataset,model.1, verbose = FALSE)    
    expect_is(mod1, 'confmirtClass')

    mod2 <- confmirt(dataset,model.1, constrain = list(c(1,5)), parprior = list(c(2, 'norm', 0, 1)),
                     verbose = FALSE)
    expect_is(mod2, 'confmirtClass')
    
    mod3 <- confmirt(dataset,model.1, itemtype = c(rep('2PL',3), '3PL', rep('graded',3), '2PL'), 
                     verbose = FALSE)
    expect_is(mod3, 'confmirtClass')

    mod.quad <- confmirt(dataset, model.quad, verbose = FALSE)
    expect_is(mod.quad, 'confmirtClass')
    
    mod.combo <- confmirt(dataset, model.combo, verbose = FALSE)
    expect_is(mod.combo, 'confmirtClass')
        
    fs1 <- fscores(mod1)
    expect_is(fs1, 'matrix')
    fs2 <- fscores(mod2)
    expect_is(fs2, 'matrix')
    fs3 <- fscores(mod.quad)
    expect_is(fs3, 'matrix')
    fs4 <- fscores(mod.combo)
    expect_is(fs4, 'matrix')    
})
 
