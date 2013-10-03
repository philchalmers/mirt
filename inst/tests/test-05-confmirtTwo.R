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
        
    model.1 <- mirt.model(model1, quiet = TRUE)    
    model.quad <- mirt.model(modelquad, quiet = TRUE)
    model.combo <- mirt.model(modelcombo, quiet = TRUE)    
    
    suppressWarnings(mod1 <- mirt(dataset,model.1, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod1, 'ConfirmatoryClass')
    expect_equal(mod1@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.3397, 0.2406, 0, NA, -1.0203, 0.0904, 0, NA, 1, NA, 0.4893, 0.0584, 0, NA, -1.4922, 0.0584, 0, NA, 1, NA, 1.5095, 0.1242, 0, NA, 1.7329, 0.025, 0, NA, 1, NA, 0.8927, 0.1551, 0.5822, 0.1019, 0.0222, 0.0747, 0, NA, 1, NA, 0, NA, 1.3949, 0.1412, 3.0503, 0.1187, 2.0627, 0.0965, -0.5311, 0.1134, 0, NA, 0.5654, 0.0808, 2.5973, 0.106, 1.0653, 0.071, -0.9324, 0.0571, 0, NA, 1.0509, 0.0851, 2.0222, 0.0914, -0.0196, 0.082, 0, NA, 1.0421, 0.166, 1.032, 0.1243, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0.3909, 0.0083, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod1b <- mirt(dataset,model.1, verbose = FALSE))
    expect_is(mod1b, 'ConfirmatoryClass')
    expect_equal(mod1b@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.3799, 0, -1.0412, 0, 1, 0.4931, 0, -1.4959, 0, 1, 1.4262, 0, 1.6824, 0, 1, 0.8877, 0.5976, 0.0155, 0, 1, 0, 1.3563, 3.012, 2.0348, -0.5284, 0, 0.5736, 2.5993, 1.0654, -0.9359, 0, 1.061, 2.0249, -0.0222, 0, 1.0538, 1.0328, 0, 1, 0, 0, 1, 0.3749, 1),
                 tollerance = 1e-2)
    
    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 586)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.72, 0.0804, 0.1646, 0.0797, -1.0007, 0.0263, 0, NA, 1, NA, 0.1791, 0.0718, 0.1442, 0.0605, -1.5885, 0.0435, 0, NA, 1, NA, 0.9014, 0.1922, 0.2022, 0.1483, 1.2499, 0.099, 0, NA, 1, NA, 1.198, 0.1188, 0.2552, 0.0451, -0.1514, 0.0636, 0, NA, 1, NA, 1.0922, 0.0374, 0, NA, 2.8028, 0.0533, 1.8784, 0.0443, -0.492, 0.0614, 0.5514, 0.0574, 0, NA, 2.5889, 0.088, 1.0584, 0.0557, -0.9317, 0.0531, 0.9165, 0.0893, 0, NA, 1.9461, 0.0753, -0.0239, 0.0553, 1.0117, 0.1118, 0, NA, 1.017, 0.075, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.4545, 0.4704, 0, NA, 0.3622, 0.266, -1.0869, 0.1665, 0, NA, 1, NA, 0.53, 0.0851, 0, NA, 0, NA, -1.5063, 0.0631, 0, NA, 1, NA, 1.4065, 0.1357, 0, NA, 0, NA, 1.6712, 0.0931, 0, NA, 1, NA, 1.0489, 0.2746, 0, NA, 0, NA, 0.0147, 0.0564, 0, NA, 1, NA, 0, NA, 1.4956, 0.1005, -0.1951, 0.1087, 3.1891, 0.1402, 2.1614, 0.1135, -0.5252, 0.091, 0, NA, 0.5436, 0.0778, 0, NA, 2.5914, 0.0882, 1.065, 0.0559, -0.9239, 0.0617, 0, NA, 1.0359, 0.087, 0, NA, 2.0201, 0.0979, -0.0138, 0.0714, 0, NA, 0.9892, 0.086, 0, NA, 1.023, 0.0816, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combob <- mirt(dataset, model.combo, verbose = FALSE))
    expect_is(mod.combob, 'ConfirmatoryClass')
    expect_equal(mod.combob@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combob, digits=4)))
    expect_equal(cfs, c(1.5249, 0, 0.4488, -1.1222, 0, 1, 0.5352, 0, 0, -1.5111, 0, 1, 1.3661, 0, 0, 1.6435, 0, 1, 1.0015, 0, 0, 0.0085, 0, 1, 0, 1.4237, -0.1974, 3.1168, 2.1096, -0.5172, 0, 0.5498, 0, 2.5916, 1.0638, -0.9272, 0, 1.0643, 0, 2.0327, -0.0166, 0, 1.0047, 0, 1.0242, 0, 1, 0, 0, 1, 0, 1),
                 tollerance = 1e-2)
        
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
 
