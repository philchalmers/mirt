context('confmirtTwo')

test_that('confirmatory mods', {
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
    model1 <- '
    F1 = 1-4
    F2 = 4-8
    COV = F1*F2'
    
    modelquad <- '
    F = 1-8
    (F*F) = 1-4
    '
    
    modelcombo <- '
    F1 = 1-4
    F2 = 5-8
    (F1*F2) = 1,5
    '    
        
    model.1 <- confmirt.model(model1, quiet = TRUE)    
    model.quad <- confmirt.model(modelquad, quiet = TRUE)
    model.combo <- confmirt.model(modelcombo, quiet = TRUE)    
    
    suppressWarnings(mod1 <- confmirt(dataset,model.1, verbose = FALSE))
    expect_is(mod1, 'ConfirmatoryClass')    
    
    suppressWarnings(mod3 <- confmirt(dataset,model.1, itemtype = c(rep('2PL',3), '3PL', rep('graded',3), '2PL'), 
                     verbose = FALSE))
    expect_is(mod3, 'ConfirmatoryClass')

    mod.quad <- confmirt(dataset, model.quad, verbose = FALSE)
    expect_is(mod.quad, 'ConfirmatoryClass')
    
    suppressWarnings(mod.combo <- confmirt(dataset, model.combo, verbose = FALSE))
    expect_is(mod.combo, 'ConfirmatoryClass')
        
    fs1 <- fscores(mod1, verbose = FALSE)
    expect_is(fs1, 'matrix')    
    fs3 <- fscores(mod.quad, full.scores=TRUE, verbose = FALSE)
    expect_is(fs3, 'matrix')
    fs4 <- fscores(mod.combo, verbose = FALSE)
    expect_is(fs4, 'matrix')
    
    TI <- plot(mod1)
    expect_is(TI, 'trellis')
    fit <- fitted(mod1)
    expect_is(fit, 'matrix')
    res <- residuals(mod1, verbose = FALSE)
    expect_is(res, 'matrix')
    IP <- itemplot(mod1, 1)
    expect_is(IP, 'trellis')
    
    TI <- plot(mod.quad)
    expect_is(TI, 'trellis')
    IP <- itemplot(mod.quad, 3, CE = TRUE)
    expect_is(IP, 'trellis')
    
    TI <- plot(mod.combo)
    expect_is(TI, 'trellis')
    IP <- itemplot(mod.combo, 1)
    expect_is(IP, 'trellis')    
})
 
