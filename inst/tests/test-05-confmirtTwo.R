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
    
    suppressWarnings(mod1 <- confmirt(dataset,model.1, verbose = FALSE, draws = 10))
    expect_is(mod1, 'ConfirmatoryClass')
    expect_equal(mod1@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.3154, 0.0715, 0, NA, -1.0109, 0.0642, 0, NA, 1, NA, 0.4721, 0.0591, 0,
                        NA, -1.4871, 0.0598, 0, NA, 1, NA, 1.5082, 0.1189, 0, NA, 1.733, 0.0492, 
                        0, NA, 1, NA, 0.9285, 0.0834, 0.5756, 0.057, 0.0239, 0.0576, 0, NA, 1, NA, 
                        0, NA, 1.3883, 0.0906, 3.0442, 0.1097, 2.0589, 0.0754, -0.5301, 0.0617, 0, 
                        NA, 0.5719, 0.0513, 2.5999, 0.0885, 1.0673, 0.0533, -0.9337, 0.0519, 0, NA, 
                        1.0401, 0.0603, 2.0151, 0.0686, -0.0196, 0.0528, 0, NA, 1.0507, 0.072, 
                        1.0352, 0.0582, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0.3926, 0.036, 1, NA),
                 tollerance = 1e-2)
    
    mod.quad <- confmirt(dataset, model.quad, verbose = FALSE, draws = 10)
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 586)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.7297, 0.0984, 0.115, 0.0459, -0.9579, 0.0785, 0, NA, 1, NA, 0.1942, 
                        0.0697, 0.1086, 0.14, -1.5497, 0.1617, 0, NA, 1, NA, 0.92, 0.0652, 0.1914, 
                        0.053, 1.2679, 0.0981, 0, NA, 1, NA, 1.1917, 0.0763, 0.1858, 0.1028, 
                        -0.1045, 0.1062, 0, NA, 1, NA, 1.1021, 0.0904, 0, NA, 2.8099, 0.0909, 
                        1.8838, 0.0643, -0.4929, 0.056, 0.5526, 0.0712, 0, NA, 2.5894, 0.0892, 
                        1.0583, 0.0544, -0.9314, 0.0514, 0.9057, 0.1425, 0, NA, 1.9411, 0.0862, 
                        -0.0213, 0.048, 0.9999, 0.1004, 0, NA, 1.0138, 0.0635, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combo <- confmirt(dataset, model.combo, verbose = FALSE, draws = 10))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.4545, 0.1909, 0, NA, 0.3825, 0.2138, -1.0907, 0.0429, 0, NA, 1, NA, 
                        0.5299, 0.1255, 0, NA, 0, NA, -1.5074, 0.0693, 0, NA, 1, NA, 1.3929, 
                        0.1927, 0, NA, 0, NA, 1.6616, 0.1127, 0, NA, 1, NA, 1.0435, 0.1415, 0, NA, 
                        0, NA, 0.0126, 0.0545, 0, NA, 1, NA, 0, NA, 1.4504, 0.0718, -0.1897, 0.0854,
                        3.1346, 0.0725, 2.1209, 0.0609, -0.5221, 0.07, 0, NA, 0.5433, 0.0709, 0, NA,
                        2.5888, 0.0892, 1.063, 0.0565, -0.9248, 0.0525, 0, NA, 1.0544, 0.0935, 0, NA,
                        2.0273, 0.0637, -0.0158, 0.0569, 0, NA, 1.0067, 0.1038, 0, NA, 1.0244, 0.0639,
                        0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0, NA, 1, NA),
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
 
