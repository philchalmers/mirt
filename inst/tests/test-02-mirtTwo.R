context('mirtTwo')

test_that('poly', {
    modp1 <- mirt(Science, 1, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')     
    expect_equal(modp1@df, 73)
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.864, 2.64, -1.466, 1.226, 2.924, 0.901, -2.267, 2.296, 5.239, 2.216, -1.965, 1.094, 3.347, 0.992, -1.688, 0, 1),
                 tollerance = 1e-2)
    vals <- mirt(Science, 1, large = TRUE, verbose=FALSE)
    modp1 <- mirt(Science, 1, large = vals, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')                  
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.864, 2.64, -1.466, 1.226, 2.924, 0.901, -2.267, 2.296, 5.239, 2.216, -1.965, 1.094, 3.347, 0.992, -1.688, 0, 1),
                 tollerance = 1e-2)    
    modp1 <- mirt(Science, 1, SE=TRUE, SE.type = 'SEM', verbose=FALSE, technical = list(TOL=1e-6))
    expect_is(modp1, 'ConfirmatoryClass')          
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 0.193, 4.863, 0.493, 2.639, 0.226, -1.465, 0.158, 1.226, 0.184, 2.924, 0.241, 0.901, 0.143, -2.266, 0.204, 2.301, 0.519, 5.245, 0.78, 2.219, 0.377, -1.967, 0.346, 1.094, 0.188, 3.347, 0.28, 0.991, 0.141, -1.688, 0.168, 0, NA, 1, NA),
                 tollerance = 1e-2)        
    suppressMessages(modp2 <- mirt(Science, 2, verbose=FALSE))
    expect_is(modp2, 'ExploratoryClass')
    expect_equal(modp2@df, 70)
    cfs <- as.numeric(do.call(c, coef(modp2, verbose = FALSE)))
    expect_equal(cfs, c(1.313, -0.052, 5.2, 2.863, -1.603, -0.128, -2.27, 3.869, 1.207, -3.04, 0.897, -1.211, 4.596, 1.926, -1.712, 1.758, 0.074, 3.984, 1.194, -2.044, 0, 0, 1, -0.494, 1),
                 tollerance = 1e-2)    
    modp3 <- mirt(Science, 1, constrain = list(c(1,5)), parprior = list(c(2,'norm',0,1)), verbose=FALSE)
    expect_is(modp3, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp3, verbose = FALSE)))
    expect_true(mirt:::closeEnough(cfs - c(1.090,  4.248,  2.550, -1.507,  1.090,  2.817,  0.853, -2.198,  2.269,  
                        5.176,  2.173, -1.978,  1.121,  3.357,  0.987, -1.714,  0.000,  1.000),
                 -1e-2, 1e-2))    
    modp4 <- suppressMessages(mirt(Science, 1, itemtype = c(rep('graded',3), 'nominal'), verbose=FALSE))
    expect_is(modp4, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp4, verbose = FALSE)))
    expect_equal(cfs, c(1.04, 4.861, 2.638, -1.466, 1.206, 2.908, 0.896, -2.254, 2.341, 5.302, 2.243, -1.991, 0.798, 0, 1.078, 1.775, 3, 0, 2.195, 2.962, 1.673, 0, 1),
                 tollerance = 1e-2)    
    modp5 <- mirt(Science, 1, itemtype = c(rep('graded',3), 'gpcm'), SE = TRUE, verbose=FALSE,
                  technical = list(TOL=1e-6))
    expect_is(modp5, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp5, verbose = FALSE)))
    expect_equal(cfs, c(1.057, 0.187, 4.876, 0.475, 2.649, 0.223, -1.472, 0.157, 1.219, 0.188, 2.918, 0.245, 0.9, 0.144, -2.263, 0.206, 2.255, 0.484, 5.178, 0.733, 2.19, 0.357, -1.943, 0.332, 0.77, 0.168, 0, NA, 2.159, 0.4, 2.973, 0.455, 1.767, 0.43, 0, NA, 1, NA),
                 tollerance = 1e-2)    
    
    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-2.7173474, -1.4189304, -0.7155405,
                                                     -0.4452374, -2.5339610, -1.2481305), -1e-4, 1e-4))
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
    expect_true(mirt:::closeEnough(as.numeric(fm2[1:6,c('F1','F2')]) - 
                                       c(-2.5585638, -1.8553231, -0.6358876, -1.1380305, -2.4068433, -0.7184635,  2.3547920,  0.7818212, -0.1434812, -0.3693318, 2.2966688, 1.5279756), -1e-4, 1e-4))
    fm3 <- fscores(modp3, rotate = 'oblimin', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modp4, verbose = FALSE)
    expect_is(fm4, 'matrix')
    expect_true(mirt:::closeEnough(fm4[1:6,'F1'] - c(-2.70139168482405, -1.44318291860916,
                                                     -0.792968643669662, -0.538984114117961, 
                                                     -2.52165495027628, -1.14464829050991), -1e-4, 1e-4))
    fm5 <- fscores(modp5, verbose = FALSE)
    expect_is(fm5, 'matrix')
    expect_true(mirt:::closeEnough(fm5[1:6,'F1'] - c(-2.6953386, -1.4445425, -0.7365539,
                                                     -0.5624948, -2.5085663, -1.1733173), -1e-4, 1e-4))
    
    cof1 <- coef(modp1)
    expect_is(cof1, 'list')
    cof2 <- coef(modp2, verbose = FALSE)
    expect_is(cof2, 'list')
    IP1 <- itemplot(modp1, 1)
    IP2 <- itemplot(modp2, 1)
    expect_is(IP1, 'trellis')
    expect_is(IP2, 'trellis')
    
    ##rating scale test
    set.seed(1234)
    a <- matrix(rep(1, 10))
    d <- matrix(c(1,0.5,-.5,-1), 10, 4, byrow = TRUE)
    cc <- seq(-1, 1, length.out=10)
    data <- simdata(a, d + cc, 2000, itemtype = rep('graded',10))
    sv <- mirt(data, 1, itemtype = 'grsm', pars = 'values', verbose=FALSE)
    sv[,'value'] <- c(as.vector(t(cbind(a,d,cc))),0,1)    
    grsm <- mirt(data, 1, itemtype = 'grsm', pars = sv, calcNull= FALSE, verbose=FALSE)
    rsm <- mirt(data, 1, itemtype = 'rsm', calcNull= FALSE, verbose=FALSE)
    expect_is(grsm, 'ConfirmatoryClass')    
    expect_is(rsm, 'ConfirmatoryClass') 
    cfs <- as.numeric(do.call(c, coef(grsm, verbose = FALSE)))
    expect_equal(cfs, c(0.998, 1.009, 0.504, -0.509, -1.003, -1, 0.862, 1.009, 0.504, -0.509, -1.003, -0.773, 1.038, 1.009, 0.504, -0.509, -1.003, -0.564, 0.949, 1.009, 0.504, -0.509, -1.003, -0.262, 0.995, 1.009, 0.504, -0.509, -1.003, -0.1, 0.905, 1.009, 0.504, -0.509, -1.003, 0.164, 0.948, 1.009, 0.504, -0.509, -1.003, 0.343, 1.091, 1.009, 0.504, -0.509, -1.003, 0.46, 0.893, 1.009, 0.504, -0.509, -1.003, 0.762, 0.967, 1.009, 0.504, -0.509, -1.003, 0.977, 0, 1),
                 tollerance = 1e-2)    
    cfs <- as.numeric(do.call(c, coef(rsm, verbose = FALSE)))
    expect_equal(cfs, c(1, 0, -2.037, -1.222, -2.047, -0.938, 0, 1, 0, -2.037, -1.222, -2.047, -0.938, 0.22, 1, 0, -2.037, -1.222, -2.047, -0.938, 0.397, 1, 0, -2.037, -1.222, -2.047, -0.938, 0.793, 1, 0, -2.037, -1.222, -2.047, -0.938, 0.938, 1, 0, -2.037, -1.222, -2.047, -0.938, 1.158, 1, 0, -2.037, -1.222, -2.047, -0.938, 1.343, 1, 0, -2.037, -1.222, -2.047, -0.938, 1.403, 1, 0, -2.037, -1.222, -2.047, -0.938, 1.758, 1, 0, -2.037, -1.222, -2.047, -0.938, 2.012, 0, 0.11),
                 tollerance = 1e-2)        
    expect_equal(rsm@df, 1945)
    expect_equal(grsm@df, 1936)
    
    #item and test info
    Theta <- matrix(seq(-4,4,.01))
    x <- extract.item(modp1, 1)
    iinfo <- iteminfo(x, Theta)
    expect_is(iinfo, 'matrix')    
    iinfo <- iteminfo(x, Theta, total.info=FALSE)
    expect_is(iinfo, 'matrix')
    tinfo <- testinfo(modp1, Theta)
    expect_is(tinfo, 'matrix')    
})