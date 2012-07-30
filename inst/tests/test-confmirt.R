context('confmirt')

test_that('all mods', {
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
    dataset <- simdata(a,d,2000,sigma)
    
    #analyses
    #CIFA for 2 factor crossed structure
        
    model.1 <- confmirt.model('confmods/model1', quiet = TRUE)
    model.2 <- confmirt.model('confmods/model2', quiet = TRUE)
    model.quad <- confmirt.model('confmods/modelquad', quiet = TRUE)
    model.combo <- confmirt.model('confmods/modelcombo', quiet = TRUE)
    model.part <- confmirt.model('confmods/noncomp3PL', quiet = TRUE)
    
    mod1 <- confmirt(dataset,model.1, verbose = FALSE)    
    expect_is(mod1, 'confmirtClass')

    mod2 <- confmirt(dataset,model.2, verbose = FALSE)
    expect_is(mod2, 'confmirtClass')

    mod.quad <- confmirt(dataset, model.quad, verbose = FALSE)
    expect_is(mod.quad, 'confmirtClass')
    
    mod.combo <- confmirt(dataset, model.combo, verbose = FALSE)
    expect_is(mod.combo, 'confmirtClass')
    
    mod.part <- confmirt(dataset, model.part, verbose = FALSE)
    expect_is(mod.part, 'confmirtClass')
    
    fs1 <- fscores(mod1)
    expect_is(fs1, 'matrix')
    fs2 <- fscores(mod2)
    expect_is(fs2, 'matrix')
    fs3 <- fscores(mod.quad)
    expect_is(fs3, 'matrix')
    fs4 <- fscores(mod.combo)
    expect_is(fs4, 'matrix')
    fs5 <- fscores(mod.part)
    expect_is(fs5, 'matrix')
})
 
