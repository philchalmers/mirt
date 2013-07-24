context('mirtTwo')

test_that('poly', {
    modp1 <- mirt(Science, 1, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')     
    expect_equal(modp1@df, 73)
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.863, 2.64, -1.466, 1.226, 2.924, 0.901, -2.267, 2.297, 5.239, 2.216, 
                        -1.965, 1.094, 3.347, 0.992, -1.688, 0, 1),
                 tollerance = 1e-2)
    vals <- mirt(Science, 1, large = TRUE, verbose=FALSE)
    modp1 <- mirt(Science, 1, large = vals, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')                  
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.863, 2.64, -1.466, 1.226, 2.924, 0.901, -2.267, 2.297, 5.239, 2.216, 
                        -1.965, 1.094, 3.347, 0.992, -1.688, 0, 1),
                 tollerance = 1e-2)    
    modp1 <- mirt(Science, 1, SE=TRUE, SE.type = 'SEM', verbose=FALSE, technical = list(TOL=1e-6))
    expect_is(modp1, 'ConfirmatoryClass')          
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 0.191, 4.863, 0.49, 2.639, 0.222, -1.465, 0.159, 1.226, 0.182, 2.924, 
                        0.24, 0.901, 0.143, -2.266, 0.203, 2.301, 0.498, 5.245, 0.757, 2.219, 0.37, 
                        -1.967, 0.331, 1.094, 0.187, 3.347, 0.277, 0.991, 0.14, -1.688, 0.169, 0, NA, 1, NA),
                 tollerance = 1e-2)        
    suppressMessages(modp2 <- mirt(Science, 2, verbose=FALSE))
    expect_is(modp2, 'ExploratoryClass')
    expect_equal(modp2@df, 70)
    cfs <- as.numeric(do.call(c, coef(modp2, verbose = FALSE)))
    expect_equal(cfs, c(1.313, -0.052, 5.2, 2.863, -1.603, -0.129, -2.269, 3.867, 1.206, -3.038, 
                        0.897, -1.211, 4.597, 1.927, -1.712, 1.758, 0.074, 3.984, 1.194, -2.044, 0, 0, 1, -0.494, 1),
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
    expect_equal(cfs, c(1.04, 4.861, 2.638, -1.466, 1.206, 2.908, 0.896, -2.254, 2.341, 5.302, 2.242,
                        -1.991, 0.798, 0, 1.078, 1.775, 3, 0, 2.195, 2.962, 1.673, 0, 1),
                 tollerance = 1e-2)    
    modp5 <- mirt(Science, 1, itemtype = c(rep('graded',3), 'gpcm'), SE = TRUE, verbose=FALSE,
                  technical = list(TOL=1e-6))
    expect_is(modp5, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp5, verbose = FALSE)))
    expect_equal(cfs, c(1.057, 0.198, 4.876, 0.493, 2.649, 0.225, -1.472, 0.162, 1.219, 0.183, 2.918, 0.241,
                        0.9, 0.144, -2.263, 0.203, 2.255, 0.5, 5.178, 0.758, 2.19, 0.373, -1.943, 0.326,
                        0.77, 0.162, 0, NA, 2.159, 0.315, 2.973, 0.353, 1.767, 0.325, 0, NA, 1, NA),
                 tollerance = 1e-2)    
    
    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-2.7173474, -1.4189304, -0.7155405,
                                                     -0.4452374, -2.5339610, -1.2481305), -1e-4, 1e-4))
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
    expect_true(mirt:::closeEnough(as.numeric(fm2[1:6,c('F1','F2')]) - 
                                       c(-2.5582887, -1.8554674, -0.6362017, -1.1384554, -2.4065656,
                                         -0.7180512,  2.3547802,  0.7824582, -0.1423948, -0.3684551,
                                         2.2966297, 1.5275234), -1e-4, 1e-4))
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
    expect_equal(cfs, c(0.997, 1.008, 0.503, -0.509, -1.003, -1, 0.862, 1.008, 0.503, -0.509, -1.003,
                        -0.773, 1.038, 1.008, 0.503, -0.509, -1.003, -0.564, 0.949, 1.008, 0.503, -0.509,
                        -1.003, -0.261, 0.995, 1.008, 0.503, -0.509, -1.003, -0.1, 0.905, 1.008, 0.503, 
                        -0.509, -1.003, 0.164, 0.948, 1.008, 0.503, -0.509, -1.003, 0.344, 1.091, 1.008, 
                        0.503, -0.509, -1.003, 0.461, 0.893, 1.008, 0.503, -0.509, -1.003, 0.763, 0.967, 
                        1.008, 0.503, -0.509, -1.003, 0.978, 0, 1),
                 tollerance = 1e-2)    
    cfs <- as.numeric(do.call(c, coef(rsm, verbose = FALSE)))
    expect_equal(cfs, c(1, 0, -2.036, -1.222, -2.047, -0.938, 0, 1, 0, -2.036, -1.222, -2.047, -0.938,
                        0.219, 1, 0, -2.036, -1.222, -2.047, -0.938, 0.396, 1, 0, -2.036, -1.222, 
                        -2.047, -0.938, 0.792, 1, 0, -2.036, -1.222, -2.047, -0.938, 0.937, 1, 0,
                        -2.036, -1.222, -2.047, -0.938, 1.157, 1, 0, -2.036, -1.222, -2.047, -0.938,
                        1.343, 1, 0, -2.036, -1.222, -2.047, -0.938, 1.402, 1, 0, -2.036, -1.222, 
                        -2.047, -0.938, 1.757, 1, 0, -2.036, -1.222, -2.047, -0.938, 2.012, 0, 0.11),
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