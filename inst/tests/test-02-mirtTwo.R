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
    expect_equal(cfs, c(1.0416, 0.1883, 4.86371, 0.49058, 2.63978, 0.22257, -1.46592, 0.15863, 1.22583, 0.1818, 2.92402, 0.23925, 0.90114, 0.14289, -2.26652, 0.20304, 2.29511, 0.48523, 5.23747, 0.73208, 2.21511, 0.35808, -1.96461, 0.32343, 1.09488, 0.1833, 3.34766, 0.27648, 0.99163, 0.14048, -1.6881, 0.16859, 0, NA, 1, NA),
                 tolerance = 1e-4)
    expect_equal(modLouis@condnum, 98.36997, tolerance = 1e-2)
    modsandwich <- mirt(Science, 1, SE=T, SE.type='sandwich', verbose=FALSE)
    expect_is(modp1, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modsandwich, digits = 5, printSE=TRUE)))
    expect_equal(cfs, c(1.0416, 0.23843, 4.86371, 0.46757, 2.63978, 0.24649, -1.46592, 0.17158, 1.22583, 0.19197, 2.92402, 0.24647, 0.90114, 0.14591, -2.26652, 0.19888, 2.29511, 0.52407, 5.23747, 0.81335, 2.21511, 0.37544, -1.96461, 0.33878, 1.09488, 0.22703, 3.34766, 0.29195, 0.99163, 0.14486, -1.6881, 0.18012, 0, NA, 1, NA),
                 tolerance = 1e-4)
    expect_equal(modsandwich@condnum, 142.5913, tolerance = 1e-2)
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
    modp5 <- mirt(Science, 1, itemtype = c(rep('graded',3), 'gpcm'), SE = TRUE, verbose=FALSE)
    expect_is(modp5, 'ConfirmatoryClass')
    cfs <- as.numeric(do.call(c, coef(modp5, verbose = FALSE)))
    expect_equal(cfs, c(1.057, 0.671, 1.442, 4.876, 3.909, 5.843, 2.649, 2.209, 3.09, -1.472, -1.789, -1.155, 1.219, 0.858, 1.58, 2.918, 2.448, 3.387, 0.9, 0.621, 1.179, -2.263, -2.663, -1.863, 2.255, 1.297, 3.213, 5.178, 3.707, 6.648, 2.19, 1.466, 2.914, -1.943, -2.574, -1.311, 0.77, 0.458, 1.083, 0, NA, NA, 1, NA, NA, 2, NA, NA, 3, NA, NA, 0, NA, NA, 2.159, 1.544, 2.775, 2.973, 2.284, 3.663, 1.767, 1.131, 2.403, 0, NA, NA, 1, NA, NA),
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

    resmat <- residuals(modp3, type = 'Q3', Theta = fm3[,'F1'], verbose = FALSE)
    expect_equal(as.numeric(resmat), c(1, -0.167, -0.144, 0.085, -0.167, 1, -0.055, -0.217, -0.144, -0.055, 1, -0.448, 0.085, -0.217, -0.448, 1))
    resmatLD <- residuals(modp3, type = 'LD', verbose = FALSE)
    expect_equal(as.numeric(resmatLD), c(NA, -20.127, -13.328, 19.607, -0.131, NA, 10.872, -21.561, -0.106, 0.096, NA, -17.541, 0.129, -0.135, -0.122, NA))
    resmatG2 <- residuals(modp3, type = 'LDG2', verbose = FALSE)
    expect_equal(as.numeric(resmatG2), c(NA, -22.903, -12.621, 23.24, -0.14, NA, 10.728, -17.242, -0.104, 0.096, NA, -18.519, 0.141, -0.121, -0.125, NA))
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