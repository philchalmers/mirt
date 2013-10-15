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
    expect_equal(cfs, c(1.3397, 0.8681, 1.8113, 0, NA, NA, -1.0203, -1.1974, -0.8431, 0, NA, NA, 1, NA, NA, 0.4893, 0.3749, 0.6037, 0, NA, NA, -1.4922, -1.6068, -1.3777, 0, NA, NA, 1, NA, NA, 1.5095, 1.2662, 1.7529, 0, NA, NA, 1.7329, 1.6839, 1.7819, 0, NA, NA, 1, NA, NA, 0.8927, 0.5886, 1.1967, 0.5822, 0.3825, 0.7818, 0.0222, -0.1243, 0.1687, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.3949, 1.1182, 1.6716, 3.0503, 2.8177, 3.2829, 2.0627, 1.8735, 2.2519, -0.5311, -0.7534, -0.3088, 0, NA, NA, 0.5654, 0.407, 0.7238, 2.5973, 2.3895, 2.8051, 1.0653, 0.9261, 1.2046, -0.9324, -1.0444, -0.8205, 0, NA, NA, 1.0509, 0.8841, 1.2177, 2.0222, 1.8431, 2.2014, -0.0196, -0.1803, 0.1412, 0, NA, NA, 1.0421, 0.7168, 1.3674, 1.032, 0.7884, 1.2756, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0.3909, 0.3745, 0.4072, 1, NA, NA),
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
    expect_equal(cfs, c(0.72, 0.5625, 0.8776, 0.1646, 0.0085, 0.3207, -1.0007, -1.0522, -0.9491, 0, NA, NA, 1, NA, NA, 0.1791, 0.0384, 0.3199, 0.1442, 0.0256, 0.2628, -1.5885, -1.6737, -1.5033, 0, NA, NA, 1, NA, NA, 0.9014, 0.5248, 1.278, 0.2022, -0.0886, 0.4929, 1.2499, 1.056, 1.4439, 0, NA, NA, 1, NA, NA, 1.198, 0.9651, 1.4308, 0.2552, 0.1668, 0.3436, -0.1514, -0.2761, -0.0268, 0, NA, NA, 1, NA, NA, 1.0922, 1.0188, 1.1656, 0, NA, NA, 2.8028, 2.6982, 2.9073, 1.8784, 1.7916, 1.9653, -0.492, -0.6124, -0.3716, 0.5514, 0.4388, 0.664, 0, NA, NA, 2.5889, 2.4165, 2.7612, 1.0584, 0.9492, 1.1675, -0.9317, -1.0357, -0.8277, 0.9165, 0.7414, 1.0916, 0, NA, NA, 1.9461, 1.7985, 2.0938, -0.0239, -0.1323, 0.0844, 1.0117, 0.7927, 1.2308, 0, NA, NA, 1.017, 0.87, 1.164, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.4545, 0.5325, 2.3765, 0, NA, NA, 0.3622, -0.1592, 0.8836, -1.0869, -1.4132, -0.7605, 0, NA, NA, 1, NA, NA, 0.53, 0.3631, 0.6968, 0, NA, NA, 0, NA, NA, -1.5063, -1.6301, -1.3826, 0, NA, NA, 1, NA, NA, 1.4065, 1.1406, 1.6724, 0, NA, NA, 0, NA, NA, 1.6712, 1.4886, 1.8537, 0, NA, NA, 1, NA, NA, 1.0489, 0.5108, 1.5871, 0, NA, NA, 0, NA, NA, 0.0147, -0.0959, 0.1252, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.4956, 1.2986, 1.6926, -0.1951, -0.4081, 0.018, 3.1891, 2.9144, 3.4639, 2.1614, 1.939, 2.3838, -0.5252, -0.7037, -0.3468, 0, NA, NA, 0.5436, 0.3912, 0.6961, 0, NA, NA, 2.5914, 2.4184, 2.7643, 1.065, 0.9554, 1.1745, -0.9239, -1.0449, -0.8029, 0, NA, NA, 1.0359, 0.8653, 1.2064, 0, NA, NA, 2.0201, 1.8283, 2.2119, -0.0138, -0.1537, 0.1261, 0, NA, NA, 0.9892, 0.8206, 1.1577, 0, NA, NA, 1.023, 0.863, 1.183, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
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
    
    TI <- plot(mod.combo)
    expect_is(TI, 'trellis')
    IP <- itemplot(mod.combo, 1)
    expect_is(IP, 'trellis')    
})
 
