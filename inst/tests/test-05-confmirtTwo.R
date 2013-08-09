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
    expect_equal(cfs, c(0.7133, 0.0306, 0.215, 0.0315, -1.0448, 0.0602, 0, NA, 1, NA, 0.1703, 0.0548, 0.1674, 0.0193, -1.6138, 0.0579, 0, NA, 1, NA, 0.9461, 0.111, 0.2572, 0.0507, 1.218, 0.0633, 0, NA, 1, NA, 1.2201, 0.118, 0.3083, 0.0796, -0.187, 0.0711, 0, NA, 1, NA, 1.0956, 0.0407, 0, NA, 2.8043, 0.1013, 1.878, 0.0781, -0.4944, 0.05, 0.5466, 0.0715, 0, NA, 2.5857, 0.0884, 1.0563, 0.0534, -0.931, 0.0522, 0.9073, 0.0645, 0, NA, 1.9394, 0.075, -0.0254, 0.0487, 1.0116, 0.105, 0, NA, 1.0154, 0.0615, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.3897, 0.1406, 0, NA, 0.323, 0.1207, -1.0661, 0.0656, 0, NA, 1, NA, 0.5282, 0.1996, 0, NA, 0, NA, -1.5064, 0.0772, 0, NA, 1, NA, 1.4684, 0.328, 0, NA, 0, NA, 1.7016, 0.1507, 0, NA, 1, NA, 1.0604, 0.1277, 0, NA, 0, NA, 0.0135, 0.0588, 0, NA, 1, NA, 0, NA, 1.4235, 0.4811, -0.18, 0.1285, 3.1183, 0.3407, 2.1112, 0.2518, -0.5154, 0.1253, 0, NA, 0.5449, 0.1713, 0, NA, 2.5908, 0.1023, 1.0641, 0.0606, -0.9248, 0.0679, 0, NA, 1.0619, 0.0961, 0, NA, 2.0342, 0.1171, -0.014, 0.0647, 0, NA, 1.0206, 0.1878, 0, NA, 1.0312, 0.1039, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0, NA, 1, NA),
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
 
