context('mirtTwo')

test_that('poly', {
    modp1 <- mirt(Science, 1, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    expect_equal(modp1@df, 239)
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.864, 2.64, -1.466, 1.226, 2.924, 0.901, -2.266, 2.296, 5.238, 2.216, -1.965, 1.095, 3.348, 0.992, -1.688, 0, 1),
                 tolerance = 1e-2)
    modLouis <- mirt(Science, 1, SE=T, SE.type='Louis', verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modLouis, digits = 5, printSE=TRUE)))
    expect_equal(cfs, c(1.04242, 0.18838, 4.86516, 0.49082, 2.64042, 0.22267, -1.46627, 0.15868, 1.2257, 0.18191, 2.9241, 0.2393, 0.90113, 0.14289, -2.26671, 0.2031, 2.28964, 0.48218, 5.2285, 0.72744, 2.21148, 0.3561, -1.96184, 0.32184, 1.0956, 0.18336, 3.34852, 0.27659, 0.99191, 0.14053, -1.68858, 0.16865, 0, NA, 1, NA),
                 tolerance = 1e-4)
    expect_equal(modLouis@condnum, 97.04485, tolerance = 1e-2)
    modsandwich <- mirt(Science, 1, SE=T, SE.type='sandwich', verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modsandwich, digits = 5, printSE=TRUE)))
    expect_equal(cfs, c(1.04242, 0.23839, 4.86516, 0.4678, 2.64042, 0.24657, -1.46627, 0.17162, 1.2257, 0.19225, 2.9241, 0.24659, 0.90113, 0.14593, -2.26671, 0.19902, 2.28964, 0.51948, 5.2285, 0.80622, 2.21148, 0.37255, -1.96184, 0.33646, 1.0956, 0.22698, 3.34852, 0.29202, 0.99191, 0.14491, -1.68858, 0.18015, 0, NA, 1, NA),
                 tolerance = 1e-4)
    expect_equal(modsandwich@condnum, 140.0871, tolerance = 1e-2)
    modp1 <- mirt(Science, 1, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    expect_equal(modp1@df, 239)
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.864, 2.64, -1.466, 1.226, 2.924, 0.901, -2.266, 2.296, 5.238, 2.216, -1.965, 1.095, 3.348, 0.992, -1.688, 0, 1),
                 tolerance = 1e-2)
    vals <- mirt(Science, 1, large = TRUE, verbose=FALSE)
    modp1 <- mirt(Science, 1, large = vals, verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 4.864, 2.64, -1.466, 1.226, 2.924, 0.901, -2.266, 2.296, 5.238, 2.216, -1.965, 1.095, 3.348, 0.992, -1.688, 0, 1),
                 tolerance = 1e-2)
    modp1 <- mirt(Science, 1, SE=TRUE, SE.type = 'SEM', verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 0.657, 1.424, 4.863, 3.343, 6.382, 2.639, 2.18, 3.099, -1.466, -1.77, -1.161, 1.226, 0.884, 1.567, 2.924, 2.454, 3.394, 0.901, 0.614, 1.189, -2.266, -2.646, -1.887, 2.3, 1.34, 3.259, 5.244, 3.837, 6.651, 2.218, 1.518, 2.919, -1.967, -2.602, -1.333, 1.094, 0.732, 1.456, 3.347, 2.808, 3.886, 0.991, 0.723, 1.26, -1.688, -2.017, -1.358, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)
    suppressMessages(modp2 <- mirt(Science, 2, verbose=FALSE))
    expect_is(modp2, 'ExploratoryClass')
    expect_equal(modp2@df, 236)
    cfs <- as.numeric(do.call(c, coef(modp2, digits=4, verbose=FALSE)))
    expect_equal(abs(cfs), abs(c(1.3248, -0.0373, 5.2047, 2.8665, -1.605, -0.1116, -2.4002, 4.02, 1.2556, -3.1606, 0.9372, -1.1447, 4.5539, 1.9064, -1.6944, 1.7547, 0.0843, 3.9762, 1.1913, -2.0395, 0, 0, 1, -0.4825, 1)),
                 tolerance = 1e-2)
    modp3 <- mirt(Science, 1, constrain = list(c(1,5)), parprior = list(c(2,'norm',0,1)), verbose=FALSE)
    expect_is(modp3, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp3, verbose = FALSE)))
    expect_true(mirt:::closeEnough(cfs - c(1.090,  4.248,  2.550, -1.507,  1.090,  2.817,  0.853, -2.198,  2.269,
                        5.176,  2.173, -1.978,  1.121,  3.357,  0.987, -1.714,  0.000,  1.000),
                 -1e-2, 1e-2))
    newmodel <- mirt.model('F = 1-4
                           CONSTRAIN = (1-2,a1)
                           PRIOR = (1, d1, norm, 0, 1)')
    modp3 <- mirt(Science, newmodel, verbose=FALSE)
    expect_is(modp3, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp3, verbose = FALSE)))
    expect_true(mirt:::closeEnough(cfs - c(1.090,  4.248,  2.550, -1.507,  1.090,  2.817,  0.853, -2.198,  2.269,
                                           5.176,  2.173, -1.978,  1.121,  3.357,  0.987, -1.714,  0.000,  1.000),
                                   -1e-2, 1e-2))

    modp4 <- suppressMessages(mirt(Science, 1, itemtype = c(rep('graded',3), 'nominal'), verbose=FALSE))
    expect_is(modp4, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp4, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(1.0408, 4.862, 2.6387, -1.4664, 1.2063, 2.9083, 0.8958, -2.254, 2.3376, 5.2972, 2.2404, -1.9886, 0.7986, 0, 1.0782, 1.7756, 3, 0, 2.1964, 2.9637, 1.6742, 0, 1),
                 tolerance = 1e-2)
    modp5 <- mirt(Science, 1, itemtype = c(rep('graded',3), 'gpcm'), SE = TRUE, SE.type = 'SEM', verbose=FALSE)
    expect_is(modp5, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp5, verbose = FALSE)))
    expect_equal(cfs, c(1.057, 0.676, 1.438, 4.876, 3.912, 5.84, 2.65, 2.21, 3.09, -1.472, -1.786, -1.159, 1.219, 0.839, 1.599, 2.918, 2.425, 3.41, 0.9, 0.615, 1.185, -2.263, -2.679, -1.848, 2.254, 1.322, 3.187, 5.177, 3.765, 6.588, 2.19, 1.49, 2.89, -1.942, -2.567, -1.318, 0.771, 0.464, 1.077, 0, NA, NA, 1, NA, NA, 2, NA, NA, 3, NA, NA, 0, NA, NA, 2.16, 1.552, 2.767, 2.973, 2.29, 3.657, 1.767, 1.137, 2.397, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)
    modp6 <- mirt(Science, 1, empiricalhist=TRUE, verbose = FALSE)
    expect_is(modp6, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp6, verbose = FALSE)))
    expect_equal(cfs, c(0.789, 5.158, 2.624, -1.344, 1.047, 2.943, 0.962, -2.184, 2.681, 5.726, 2.782, -1.808, 0.931, 3.423, 1.047, -1.592, 0, 1),
                 tolerance = 1e-2)
    
    fm0 <- fscores(modp1, method='EAP', response.pattern = c(1,2,3,4))
    expect_equal(as.numeric(fm0[,c('F1','SE_F1')]), c(-0.3494903, 0.6004922), tolerance=1e-4)
    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-2.7173474, -1.4189304, -0.7155405,
                                                     -0.4452374, -2.5339610, -1.2481305), -1e-2, 1e-2))
    fm1b <- fscores(modp1, verbose = FALSE, full.scores=TRUE)
    expect_equal(cor(fm1b, rowSums(Science))[1], .969, tolerance = .02)
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
    expect_true(mirt:::closeEnough(abs(as.numeric(fm2[1:6,c('F1','F2')])) -
                                       abs(c(2.5707445, 1.8477329, 0.6198135, 1.1152277, 2.4184494,
                                             0.7455961, 2.3420774, 0.7085009, 0.2521220, 0.4697231, 2.2881277, 1.5731773)),
                                   -1e-2, 1e-2))
    fm3 <- fscores(modp3, rotate = 'oblimin', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modp4, verbose = FALSE)
    expect_is(fm4, 'matrix')
    expect_true(mirt:::closeEnough(fm4[1:6,'F1'] - c(-2.7025522, -1.4458413, -0.7908526, -0.5438003, -2.5218037, -1.1433440), -1e-4, 1e-4))
    fm5 <- fscores(modp5, verbose = FALSE)
    expect_is(fm5, 'matrix')
    expect_true(mirt:::closeEnough(fm5[1:6,'F1'] - c(-2.6953386, -1.4445425, -0.7365539,
                                                     -0.5624948, -2.5085663, -1.1733173), -1e-2, 1e-2))

    resmat <- residuals(modp3, type = 'Q3', Theta = fm3[,'F1'], verbose = FALSE)
    expect_equal(as.numeric(resmat), c(1, -0.167, -0.144, 0.084, -0.167, 1, -0.055, -0.218, -0.144, -0.055, 1, -0.447, 0.084, -0.218, -0.447, 1))
    resmatLD <- residuals(modp3, type = 'LD', verbose = FALSE)
    expect_equal(as.numeric(resmatLD), c(NA, -20.136, -13.303, 19.573, -0.131, NA, 10.895, -21.592, -0.106, 0.096, NA, -17.521, 0.129, -0.136, -0.122, NA))
    resmatG2 <- residuals(modp3, type = 'LDG2', verbose = FALSE)
    expect_equal(as.numeric(resmatG2), c(NA, -22.906, -12.609, 23.21, -0.14, NA, 10.753, -17.254, -0.104, 0.096, NA, -18.509, 0.14, -0.121, -0.125, NA))
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
    rsm <- mirt(data, 1, itemtype = 'rsm', calcNull= FALSE, verbose=FALSE, TOL = 1e-3)
    expect_is(grsm, 'ConfirmatoryClass')
    expect_is(rsm, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(grsm, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(0.997, 1.008, 0.5031, -0.5092, -1.0034, -1, 0.8616, 1.008, 0.5031, -0.5092, -1.0034, -0.7723, 1.0369, 1.008, 0.5031, -0.5092, -1.0034, -0.5635, 0.9486, 1.008, 0.5031, -0.5092, -1.0034, -0.2613, 0.9941, 1.008, 0.5031, -0.5092, -1.0034, -0.0998, 0.9047, 1.008, 0.5031, -0.5092, -1.0034, 0.1646, 0.9475, 1.008, 0.5031, -0.5092, -1.0034, 0.3439, 1.09, 1.008, 0.5031, -0.5092, -1.0034, 0.4608, 0.8921, 1.008, 0.5031, -0.5092, -1.0034, 0.7628, 0.966, 1.008, 0.5031, -0.5092, -1.0034, 0.9779, 0, 1),
                 tolerance = 1e-2)
    cfs <- as.numeric(do.call(c, coef(rsm, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 0, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 0.2204, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 0.3977, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 0.7937, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 0.9389, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 1.1589, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 1.3442, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 1.4039, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 1.7589, 1, 0, 1, 2, 3, 4, 0, -2.0381, -1.224, -2.05, -0.9416, 2.0131, 0, 0.1104),
                 tolerance = 1e-2)
    expect_equal(rsm@df, 9765610)
    expect_equal(grsm@df, 9765601)

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
