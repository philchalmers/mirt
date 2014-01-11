context('mirtTwo')

test_that('poly', {
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
    modp1 <- mirt(Science, 1, SE=TRUE, SE.type = 'SEM', verbose=FALSE, technical = list(TOL=1e-6))
    expect_is(modp1, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp1)))
    expect_equal(cfs, c(1.041, 0.669, 1.412, 4.863, 3.904, 5.822, 2.639, 2.203, 3.075, -1.465, -1.775, -1.156, 1.226, 0.866, 1.586, 2.924, 2.442, 3.406, 0.901, 0.617, 1.185, -2.266, -2.663, -1.869, 2.301, 1.334, 3.267, 5.245, 3.782, 6.708, 2.219, 1.506, 2.932, -1.967, -2.613, -1.322, 1.094, 0.731, 1.456, 3.347, 2.804, 3.89, 0.991, 0.716, 1.267, -1.688, -2.016, -1.359, 0, NA, NA, 1, NA, NA),
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
    modp5 <- mirt(Science, 1, itemtype = c(rep('graded',3), 'gpcm'), SE = TRUE, verbose=FALSE,
                  technical = list(TOL=1e-6))
    expect_is(modp5, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp5, verbose = FALSE)))
    expect_equal(cfs, c(1.057, 0.674, 1.439, 4.876, 3.913, 5.84, 2.649, 2.21, 3.089, -1.472, -1.789, -1.156, 1.219, 0.856, 1.582, 2.918, 2.445, 3.39, 0.9, 0.619, 1.18, -2.263, -2.662, -1.864, 2.255, 1.289, 3.221, 5.178, 3.712, 6.643, 2.19, 1.47, 2.911, -1.943, -2.575, -1.311, 0.77, 0.464, 1.076, 0, NA, NA, 2.159, 1.552, 2.767, 2.973, 2.294, 3.653, 1.767, 1.136, 2.398, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)
    modp6 <- mirt(Science, 1, empiricalhist=TRUE, verbose = FALSE)
    expect_is(modp6, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp6, verbose = FALSE)))
    expect_equal(cfs, c(0.764, 5.211, 2.662, -1.281, 1.031, 2.999, 1.025, -2.095, 3.172, 6.42, 3.44, -1.639, 0.924, 3.479, 1.109, -1.518, 0, 1),
                 tolerance = 1e-2)

    fm1 <- fscores(modp1, verbose = FALSE)
    expect_is(fm1, 'matrix')
    expect_true(mirt:::closeEnough(fm1[1:6,'F1'] - c(-2.7173474, -1.4189304, -0.7155405,
                                                     -0.4452374, -2.5339610, -1.2481305), -1e-4, 1e-4))
    fm2 <- fscores(modp2, rotate = 'oblimin', verbose = FALSE)
    expect_is(fm2, 'matrix')
    expect_true(mirt:::closeEnough(abs(as.numeric(fm2[1:6,c('F1','F2')])) -
                                       abs(c(-2.5689536, -1.8493210, -0.6202252, -1.1195381, -2.4169928, -0.7386485,  2.3431580,  0.7169414, -0.2393583, -0.4576444,  2.2886770,  1.5675075)),
                                   -1e-2, 1e-2))
    fm3 <- fscores(modp3, rotate = 'oblimin', full.scores = TRUE, verbose = FALSE)
    expect_is(fm3, 'matrix')
    fm4 <- fscores(modp4, verbose = FALSE)
    expect_is(fm4, 'matrix')
    expect_true(mirt:::closeEnough(fm4[1:6,'F1'] - c(-2.7014033, -1.4432191, -0.7929087, -0.5390761, -2.5216331, -1.1445815), -1e-4, 1e-4))
    fm5 <- fscores(modp5, verbose = FALSE)
    expect_is(fm5, 'matrix')
    expect_true(mirt:::closeEnough(fm5[1:6,'F1'] - c(-2.6953386, -1.4445425, -0.7365539,
                                                     -0.5624948, -2.5085663, -1.1733173), -1e-4, 1e-4))

    resmat <- residuals(modp3, restype = 'Q3', Theta = fm3[,'F1'], verbose = FALSE)
    expect_equal(as.numeric(resmat), c(1, -0.167, -0.144, 0.085, -0.167, 1, -0.055, -0.217, -0.144, -0.055, 1, -0.448, 0.085, -0.217, -0.448, 1))
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
    expect_equal(cfs, c(0.997, 1.008, 0.5031, -0.5092, -1.0034, -1, 0.8616, 1.008, 0.5031, -0.5092, -1.0034, -0.7723, 1.0369, 1.008, 0.5031, -0.5092, -1.0034, -0.5635, 0.9486, 1.008, 0.5031, -0.5092, -1.0034, -0.2613, 0.9941, 1.008, 0.5031, -0.5092, -1.0034, -0.0998, 0.9047, 1.008, 0.5031, -0.5092, -1.0034, 0.1646, 0.9475, 1.008, 0.5031, -0.5092, -1.0034, 0.3439, 1.09, 1.008, 0.5031, -0.5092, -1.0034, 0.4608, 0.8921, 1.008, 0.5031, -0.5092, -1.0034, 0.7628, 0.966, 1.008, 0.5031, -0.5092, -1.0034, 0.9779, 0, 1),
                 tolerance = 1e-2)
    cfs <- as.numeric(do.call(c, coef(rsm, verbose = FALSE, digits=4)))
    expect_equal(cfs, c(1, 0, -2.0346, -1.22, -2.045, -0.9358, 0, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 0.2175, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 0.3947, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 0.7905, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 0.9357, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 1.1556, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 1.3409, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 1.4005, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 1.7554, 1, 0, -2.0346, -1.22, -2.045, -0.9358, 2.0096, 0, 0.11),
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