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
    expect_equal(cfs, c(1.3482, 0.4221, 0, NA, -1.0264, 0.1148, 0, NA, 1, NA, 0.484, 0.0344, 0, NA, -1.4921, 0.0569, 0, NA, 1, NA, 1.458, 0.2411, 0, NA, 1.7049, 0.0917, 0, NA, 1, NA, 0.8998, 0.0679, 0.5763, 0.1288, 0.0209, 0.0561, 0, NA, 1, NA, 0, NA, 1.3933, 0.2008, 3.0463, 0.0613, 2.06, 0.0504, -0.5314, 0.063, 0, NA, 0.5678, 0.0633, 2.598, 0.0886, 1.0659, 0.0545, -0.9328, 0.0488, 0, NA, 1.0444, 0.1057, 2.0169, 0.0848, -0.0203, 0.0395, 0, NA, 1.0579, 0.0831, 1.0363, 0.0344, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0.406, 0.1683, 1, NA),
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
    expect_equal(cfs, c(0.7184, 0.0949, 0.1423, 0.0129, -0.9815, 0.0699, 0, NA, 1, NA, 0.1853, 0.0742, 0.1244, 0.5222, -1.5669, 0.5897, 0, NA, 1, NA, 0.9165, 0.3787, 0.2013, 0.2966, 1.2562, 0.2814, 0, NA, 1, NA, 1.196, 0.2022, 0.2208, 0.3461, -0.1295, 0.1876, 0, NA, 1, NA, 1.101, 0.6199, 0, NA, 2.8071, 0.4663, 1.8802, 0.3403, -0.4946, 0.0406, 0.5501, 0.1051, 0, NA, 2.5876, 0.0822, 1.057, 0.0532, -0.9322, 0.0275, 0.9072, 0.0539, 0, NA, 1.9396, 0.0752, -0.0251, 0.03, 1.0053, 0.1328, 0, NA, 1.0139, 0.0602, 0, NA, 1, NA, 0, NA, 1, NA),
                 tollerance = 1e-2)
    
    suppressWarnings(mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 588)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.5184, 0.1011, 0, NA, 0.4103, 0.2855, -1.1117, 0.0928, 0, NA, 1, NA, 0.5425, 0.1393, 0, NA, 0, NA, -1.511, 0.0763, 0, NA, 1, NA, 1.4324, 0.1564, 0, NA, 0, NA, 1.6807, 0.1487, 0, NA, 1, NA, 0.9981, 0.1085, 0, NA, 0, NA, 0.0128, 0.0593, 0, NA, 1, NA, 0, NA, 1.4742, 0.277, -0.1972, 0.0609, 3.1639, 0.2312, 2.1423, 0.1626, -0.5258, 0.0927, 0, NA, 0.5444, 0.1152, 0, NA, 2.5897, 0.0933, 1.0632, 0.0566, -0.9262, 0.0586, 0, NA, 1.0499, 0.0949, 0, NA, 2.0246, 0.102, -0.0168, 0.0579, 0, NA, 0.9994, 0.1224, 0, NA, 1.0222, 0.0772, 0, NA, 1, NA, 0, NA, 0, NA, 1, NA, 0, NA, 1, NA),
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
 
