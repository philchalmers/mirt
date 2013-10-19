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
    expect_equal(cfs, c(1.3353, 1.0643, 1.6063, 0, NA, NA, -1.0209, -1.1738, -0.868, 0, NA, NA, 1, NA, NA, 0.4913, 0.3074, 0.6753, 0, NA, NA, -1.4936, -1.616, -1.3711, 0, NA, NA, 1, NA, NA, 1.5247, 1.4077, 1.6418, 0, NA, NA, 1.7391, 1.6165, 1.8616, 0, NA, NA, 1, NA, NA, 0.8886, 0.6148, 1.1624, 0.5785, 0.3752, 0.7819, 0.0199, -0.1018, 0.1416, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.3988, 1.2066, 1.591, 3.0458, 2.8186, 3.2731, 2.0588, 1.8794, 2.2383, -0.5336, -0.673, -0.3941, 0, NA, NA, 0.5679, 0.4383, 0.6975, 2.5963, 2.4208, 2.7718, 1.0646, 0.9543, 1.1748, -0.9336, -1.0414, -0.8257, 0, NA, NA, 1.0514, 0.8607, 1.2422, 2.0187, 1.8345, 2.203, -0.0214, -0.1373, 0.0946, 0, NA, NA, 1.0436, 0.8196, 1.2676, 1.0294, 0.8807, 1.1781, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0.3903, 0.3809, 0.3998, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod1b <- mirt(dataset,model.1, verbose = FALSE))
    expect_is(mod1b, 'ConfirmatoryClass')
    expect_equal(mod1b@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.3802, 0, -1.0413, 0, 1, 0.4931, 0, -1.4959, 0, 1, 1.4261, 0, 1.6823, 0, 1, 0.8878, 0.5975, 0.0155, 0, 1, 0, 1.3563, 3.012, 2.0348, -0.5284, 0, 0.5736, 2.5992, 1.0654, -0.9359, 0, 1.0611, 2.0249, -0.0222, 0, 1.0538, 1.0328, 0, 1, 0, 0, 1, 0.3749, 1),
                 tollerance = 1e-2)
    
    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 586)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.7082, 0.5326, 0.8837, 0.1967, 0.121, 0.2725, -1.0261, -1.1296, -0.9226, 0, NA, NA, 1, NA, NA, 0.1762, 0.0564, 0.296, 0.1639, 0.0035, 0.3244, -1.6103, -1.7454, -1.4752, 0, NA, NA, 1, NA, NA, 0.927, 0.6217, 1.2322, 0.2409, 0.0659, 0.416, 1.2269, 1.0634, 1.3904, 0, NA, NA, 1, NA, NA, 1.1917, 0.9639, 1.4195, 0.2804, 0.1532, 0.4075, -0.1673, -0.3013, -0.0334, 0, NA, NA, 1, NA, NA, 1.101, 0.9485, 1.2535, 0, NA, NA, 2.8099, 2.6092, 3.0107, 1.883, 1.7253, 2.0407, -0.4926, -0.6114, -0.3738, 0.5472, 0.4193, 0.675, 0, NA, NA, 2.5873, 2.4203, 2.7542, 1.0574, 0.9553, 1.1595, -0.9304, -1.0398, -0.821, 0.9143, 0.715, 1.1135, 0, NA, NA, 1.9453, 1.796, 2.0947, -0.023, -0.1333, 0.0873, 1.0223, 0.8199, 1.2248, 0, NA, NA, 1.0213, 0.8854, 1.1573, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.4575, 1.2148, 1.7002, 0, NA, NA, 0.3736, -0.0422, 0.7893, -1.0956, -1.2711, -0.9201, 0, NA, NA, 1, NA, NA, 0.5204, 0.3, 0.7407, 0, NA, NA, 0, NA, NA, -1.506, -1.6231, -1.389, 0, NA, NA, 1, NA, NA, 1.4134, 1.2336, 1.5932, 0, NA, NA, 0, NA, NA, 1.6664, 1.515, 1.8178, 0, NA, NA, 1, NA, NA, 1.0597, 0.8119, 1.3075, 0, NA, NA, 0, NA, NA, 0.0092, -0.0814, 0.0999, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.5029, 1.4057, 1.6002, -0.1896, -1.0428, 0.6635, 3.1876, 2.9085, 3.4667, 2.1608, 1.9445, 2.377, -0.5254, -0.6342, -0.4167, 0, NA, NA, 0.5493, 0.4088, 0.6897, 0, NA, NA, 2.5928, 2.4179, 2.7677, 1.0659, 0.9627, 1.169, -0.9246, -1.0198, -0.8293, 0, NA, NA, 1.0255, 0.5366, 1.5145, 0, NA, NA, 2.0125, 1.6273, 2.3977, -0.0136, -0.0879, 0.0607, 0, NA, NA, 1.0015, 0.7865, 1.2166, 0, NA, NA, 1.0253, 0.9418, 1.1088, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combob <- mirt(dataset, model.combo, verbose = FALSE))
    expect_is(mod.combob, 'ConfirmatoryClass')
    expect_equal(mod.combob@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combob, digits=4)))
    expect_equal(cfs, c(1.5244, 0, 0.4483, -1.1219, 0, 1, 0.5352, 0, 0, -1.5111, 0, 1, 1.3662, 0, 0, 1.6436, 0, 1, 1.0016, 0, 0, 0.0085, 0, 1, 0, 1.4241, -0.1975, 3.1172, 2.1099, -0.5173, 0, 0.5498, 0, 2.5916, 1.0638, -0.9273, 0, 1.0643, 0, 2.0327, -0.0166, 0, 1.0047, 0, 1.0242, 0, 1, 0, 0, 1, 0, 1),
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
 
