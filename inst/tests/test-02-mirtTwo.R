context('mirtTwo')

test_that('poly', {
    modp1 <- mirt(Science, 1, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')     
    expect_equal(modp1@df, 73)
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.042,  4.864,  2.640, -1.466,  1.226,  2.924,  0.901, -2.267,  2.294,
                        5.235,  2.214, -1.964,  1.095,  3.348,  0.992, -1.688,  0.000,  1.000),
                 tollerance = 1e-3)
    vals <- mirt(Science, 1, large = TRUE, verbose=FALSE)
    modp1 <- mirt(Science, 1, large = vals, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')                  
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.042,  4.864,  2.640, -1.466,  1.226,  2.924,  0.901, -2.267,  2.294,
                        5.235,  2.214, -1.964,  1.095,  3.348,  0.992, -1.688,  0.000,  1.000),
                 tollerance = 1e-3)    
    modp1 <- mirt(Science, 1, SE=TRUE, SE.type = 'SEM', verbose=FALSE, technical = list(TOL=1e-6))
    expect_is(modp1, 'ConfirmatoryClass')          
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041,  0.191,  4.863,  0.491,  2.639,  0.223, -1.465,  0.159,  1.226,
                        0.181,  2.924,  0.239,  0.901,  0.143, -2.266,  0.204,  2.301,  0.501,
                        5.245,  0.764,  2.219,  0.374, -1.967,  0.332,  1.094,  0.187,  3.347,
                        0.277,  0.991,  0.140, -1.688,  0.170,  0.000,     NA,  1.000,     NA),
                 tollerance = 1e-3)        
    suppressMessages(modp2 <- mirt(Science, 2, verbose=FALSE))
    expect_is(modp2, 'ExploratoryClass')
    expect_equal(modp2@df, 70)
    cfs <- as.numeric(do.call(c, coef(modp2, verbose = FALSE)))
    expect_equal(cfs, c(1.320, -0.042,  5.198,  2.864, -1.604, -0.115, -2.375,  3.989,  1.246,
                        -3.137,  0.930, -1.158,  4.560,  1.910, -1.698,  1.760,  0.083,
                        3.981,  1.193, -2.043,  0.000,  0.000,  1.000, -0.485,  1.000),
                 tollerance = 1e-3)    
    modp3 <- mirt(Science, 1, constrain = list(c(1,5)), parprior = list(c(2,'norm',0,1)), verbose=FALSE)
    expect_is(modp3, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp3, verbose = FALSE)))
    expect_equal(cfs, c(1.090,  4.248,  2.550, -1.507,  1.090,  2.817,  0.853, -2.198,  2.269,  
                        5.176,  2.173, -1.978,  1.121,  3.357,  0.987, -1.714,  0.000,  1.000),
                 tollerance = 1e-3)    
    modp4 <- suppressMessages(mirt(Science, 1, itemtype = c(rep('graded',3), 'nominal'), verbose=FALSE))
    expect_is(modp4, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp4, verbose = FALSE)))
    expect_equal(cfs, c(1.040,  4.862,  2.638, -1.466,  1.206,  2.908,  0.896, -2.254,  2.339,
                        5.297,  2.241, -1.989,  0.798,  0.000,  1.079,  1.776,  3.000,
                        0.000,  2.197,  2.964,  1.675,  0.000,  1.000),
                 tollerance = 1e-3)    
    modp5 <- mirt(Science, 1, itemtype = c(rep('graded',3), 'gpcm'), SE = TRUE, verbose=FALSE,
                  technical = list(TOL=1e-6))
    expect_is(modp5, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp5, verbose = FALSE)))
    expect_equal(cfs, c(1.057,  0.199,  4.876,  0.494,  2.649,  0.226, -1.472,  0.162,  1.219,
                        0.183,  2.918,  0.241,  0.900,  0.144, -2.263,  0.204,  2.255,
                        0.499,  5.178,  0.758,  2.190,  0.373, -1.943,  0.327,  0.770,  0.163,
                        0.000,     NA,  2.159,  0.316,  2.973,  0.356,  1.767,  0.327,
                        0.000,     NA,  1.000,     NA),
                 tollerance = 1e-3)    
    
    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-2.7173474, -1.4189304, -0.7155405,
                                                     -0.4452374, -2.5339610, -1.2481305), -1e-4, 1e-4))
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
    expect_true(mirt:::closeEnough(as.numeric(fm2[1:6,c('F1','F2')]) - 
                                       c(-2.5675049, -1.8503899, -0.6211492, -1.1228684, -2.4158367,
                                         -0.7326870,  2.3463451,  0.7306712, -0.2198164, -0.4388248,
                                         2.2909837, 1.5594119), -1e-4, 1e-4))
    fm3 <- fscores(modp3, rotate = 'oblimin', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modp4, verbose = FALSE)
    expect_is(fm4, 'matrix')
    expect_true(mirt:::closeEnough(fm4[1:6,'F1'] - c(-2.7016115, -1.4438297, -0.7924922,
                                                     -0.5401814, -2.5217015, -1.1443791), -1e-4, 1e-4))
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
    expect_equal(cfs, c(0.979,  1.003,  0.498, -0.514, -1.009, -1.000,  0.843,  1.003,  
                        0.498, -0.514, -1.009, -0.767,  1.014,  1.003,  0.498, -0.514, -1.009,
                       -0.558,  0.928,  1.003,  0.498, -0.514, -1.009, -0.256,  0.972,  
                        1.003,  0.498, -0.514, -1.009, -0.095,  0.885,  1.003,  0.498, -0.514,
                       -1.009,  0.170,  0.926,  1.003,  0.498, -0.514, -1.009,  0.349,  
                        1.065,  1.003,  0.498, -0.514, -1.009,  0.466,  0.873,  1.003,  0.498,
                       -0.514, -1.009,  0.768,  0.946,  1.003,  0.498, -0.514, -1.009,  0.983,  0.000,  1.045),
                 tollerance = 1e-3)    
    cfs <- as.numeric(do.call(c, coef(rsm, verbose = FALSE)))
    expect_equal(cfs, c(1.000,  0.000, -2.037, -1.223, -2.048, -0.939,  0.000,  1.000,  0.000, 
                        -2.037, -1.223, -2.048, -0.939,  0.220,  1.000,  0.000, -2.037,
                       -1.223, -2.048, -0.939,  0.397,  1.000,  0.000, -2.037, -1.223, -2.048, 
                        -0.939,  0.793,  1.000,  0.000, -2.037, -1.223, -2.048, -0.939,
                        0.939,  1.000,  0.000, -2.037, -1.223, -2.048, -0.939,  1.159,  1.000,  
                        0.000, -2.037, -1.223, -2.048, -0.939,  1.344,  1.000,  0.000,
                       -2.037, -1.223, -2.048, -0.939,  1.404,  1.000,  0.000, -2.037, -1.223, 
                        -2.048, -0.939,  1.759,  1.000,  0.000, -2.037, -1.223, -2.048,
                       -0.939,  2.013,  0.000,  0.110),
                 tollerance = 1e-3)        
    expect_equal(rsm@df, 1945)
    expect_equal(grsm@df, 1935)
    
    #item and test info
    Theta <- matrix(seq(-4,4,.01))
    x <- extract.item(modp1, 1)
    iinfo <- iteminfo(x, Theta)
    expect_is(iinfo, 'matrix')    
    tinfo <- testinfo(modp1, Theta)
    expect_is(tinfo, 'matrix')    
})