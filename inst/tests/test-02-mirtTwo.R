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
    cfs <- as.numeric(do.call(c, coef(modp2, digits=4, verbose=FALSE)))
    expect_equal(abs(cfs), abs(c(1.31, -0.0563, 5.199, 2.863, -1.6029, -0.1337, -2.2349, 3.8279, 1.1934, -3.0068, 0.8848, -1.2312, 4.6096, 1.9326, -1.7172, 1.7579, 0.0709, 3.9851, 1.1943, -2.0445, 0, 0, 1, -0.4981, 1)),
                 tollerance = 1e-2)    
    modp3 <- mirt(Science, 1, constrain = list(c(1,5)), parprior = list(c(2,'norm',0,1)), verbose=FALSE)
    expect_is(modp3, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp3, verbose = FALSE)))
    expect_true(mirt:::closeEnough(cfs - c(1.090,  4.248,  2.550, -1.507,  1.090,  2.817,  0.853, -2.198,  2.269,  
                        5.176,  2.173, -1.978,  1.121,  3.357,  0.987, -1.714,  0.000,  1.000),
                 -1e-2, 1e-2))    
    modp4 <- suppressMessages(mirt(Science, 1, itemtype = c(rep('graded',3), 'nominal'), verbose=FALSE))
    expect_is(modp4, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp4, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(1.0399, 4.8613, 2.6381, -1.4661, 1.2064, 2.9082, 0.8957, -2.2541, 2.3415, 5.3016, 2.2425, -1.9909, 0.7975, 0, 1.0775, 1.7747, 3, 0, 2.1944, 2.9615, 1.6723, 0, 1),
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
    expect_true(mirt:::closeEnough(abs(as.numeric(fm2[1:6,c('F1','F2')])) - 
                                       abs(c(-2.5558577, -1.8568859, -0.6405448, -1.1427467, -2.4041356, -0.7137254,  2.3594542,  0.8015557, -0.1156903,
                                         -0.3435933,  2.3002348,  1.5171179)),
                                   -1e-4, 1e-4))
    fm3 <- fscores(modp3, rotate = 'oblimin', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modp4, verbose = FALSE)
    expect_is(fm4, 'matrix')
    expect_true(mirt:::closeEnough(fm4[1:6,'F1'] - c(-2.7012919, -1.4429877, -0.7929602, -0.5387171, -2.5215738, -1.1446546),
                                   -1e-4, 1e-4))
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
    cfs <- as.numeric(do.call(c, coef(grsm, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(0.9976, 1.0086, 0.5036, -0.5086, -1.0029, -1, 0.8623, 1.0086, 0.5036, -0.5086, -1.0029, -0.773, 1.0377, 1.0086, 0.5036, -0.5086, -1.0029, -0.5642, 0.9493, 1.0086, 0.5036, -0.5086, -1.0029, -0.2619, 0.9949, 1.0086, 0.5036, -0.5086, -1.0029, -0.1004, 0.9054, 1.0086, 0.5036, -0.5086, -1.0029, 0.164, 0.9482, 1.0086, 0.5036, -0.5086, -1.0029, 0.3433, 1.0908, 1.0086, 0.5036, -0.5086, -1.0029, 0.4602, 0.8928, 1.0086, 0.5036, -0.5086, -1.0029, 0.7623, 0.9668, 1.0086, 0.5036, -0.5086, -1.0029, 0.9775, 0, 1),
                 tollerance = 1e-2)    
    cfs <- as.numeric(do.call(c, coef(rsm, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 0, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 0.2165, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 0.3938, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 0.7899, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 0.9351, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 1.155, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 1.3404, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 1.4, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 1.7549, 1, 0, -2.0339, -1.2193, -2.0444, -0.9355, 2.0091, 0, 0.1103),
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