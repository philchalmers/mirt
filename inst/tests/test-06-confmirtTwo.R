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
    expect_equal(mod1@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.3353, 1.0286, 1.6421, 0, NA, NA, -1.0209, -1.18, -0.8617, 0, NA, NA, 1, NA, NA, 0.4913, 0.2731, 0.7095, 0, NA, NA, -1.4936, -1.6233, -1.3638, 0, NA, NA, 1, NA, NA, 1.5247, 0.9073, 2.1422, 0, NA, NA, 1.7391, 1.4144, 2.0637, 0, NA, NA, 1, NA, NA, 0.8886, 0.5608, 1.2164, 0.5785, 0.3267, 0.8304, 0.0199, -0.1018, 0.1416, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.3988, 1.1487, 1.6489, 3.0458, 2.7906, 3.3011, 2.0588, 1.8635, 2.2542, -0.5336, -0.6773, -0.3898, 0, NA, NA, 0.5679, 0.4354, 0.7004, 2.5963, 2.4186, 2.774, 1.0646, 0.9525, 1.1767, -0.9336, -1.039, -0.8281, 0, NA, NA, 1.0514, 0.8387, 1.2641, 2.0187, 1.8429, 2.1945, -0.0214, -0.1379, 0.0951, 0, NA, NA, 1.0436, 0.8277, 1.2595, 1.0294, 0.8876, 1.1712, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0.3903, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    suppressWarnings(mod1b <- mirt(dataset,model.1, verbose = FALSE))
    expect_is(mod1b, 'ConfirmatoryClass')
    expect_equal(mod1b@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.3794, 0, -1.041, 0, 1, 0.4917, 0, -1.4955, 0, 1, 1.4234, 0, 1.6808, 0, 1, 0.8878, 0.5941, 0.0156, 0, 1, 0, 1.356, 3.0118, 2.0346, -0.5283, 0, 0.5739, 2.5993, 1.0655, -0.9359, 0, 1.0602, 2.0243, -0.0221, 0, 1.0539, 1.0328, 0, 1, 0, 0, 1, 0.3815, 1),
                 tolerance = 1e-2)

    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 1510)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.7082, 0.5163, 0.9001, 0.1967, 0.1189, 0.2746, -1.0261, -1.1473, -0.905, 0, NA, NA, 1, NA, NA, 0.1762, 0.0107, 0.3416, 0.1639, -0.1033, 0.4311, -1.6103, -1.9317, -1.2889, 0, NA, NA, 1, NA, NA, 0.927, 0.598, 1.256, 0.2409, 0.0554, 0.4265, 1.2269, 1.0521, 1.4018, 0, NA, NA, 1, NA, NA, 1.1917, 0.922, 1.4614, 0.2804, -0.0373, 0.5981, -0.1673, -0.351, 0.0163, 0, NA, NA, 1, NA, NA, 1.101, 0.944, 1.258, 0, NA, NA, 2.8099, 2.6072, 3.0127, 1.883, 1.7182, 2.0479, -0.4926, -0.6289, -0.3562, 0.5472, 0.419, 0.6754, 0, NA, NA, 2.5873, 2.4183, 2.7563, 1.0574, 0.9505, 1.1643, -0.9304, -1.0448, -0.8161, 0.9143, 0.6788, 1.1497, 0, NA, NA, 1.9453, 1.796, 2.0947, -0.023, -0.1477, 0.1017, 1.0223, 0.8185, 1.2262, 0, NA, NA, 1.0213, 0.8745, 1.1682, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    suppressWarnings(mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM'))
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.4575, 0.8709, 2.044, 0, NA, NA, 0.3736, -0.3325, 1.0796, -1.0956, -1.4419, -0.7493, 0, NA, NA, 1, NA, NA, 0.5204, 0.1871, 0.8537, 0, NA, NA, 0, NA, NA, -1.506, -1.6404, -1.3717, 0, NA, NA, 1, NA, NA, 1.4134, 0.525, 2.3018, 0, NA, NA, 0, NA, NA, 1.6664, 1.1831, 2.1498, 0, NA, NA, 1, NA, NA, 1.0597, 0.8049, 1.3146, 0, NA, NA, 0, NA, NA, 0.0092, -0.116, 0.1345, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.5029, 1.1664, 1.8395, -0.1896, -1.1032, 0.7239, 3.1876, 2.8886, 3.4865, 2.1608, 1.9401, 2.3815, -0.5254, -0.7062, -0.3447, 0, NA, NA, 0.5493, 0.3881, 0.7104, 0, NA, NA, 2.5928, 2.4156, 2.77, 1.0659, 0.954, 1.1777, -0.9246, -1.0328, -0.8163, 0, NA, NA, 1.0255, 0.3063, 1.7448, 0, NA, NA, 2.0125, 1.5081, 2.5168, -0.0136, -0.1318, 0.1045, 0, NA, NA, 1.0015, 0.7831, 1.2199, 0, NA, NA, 1.0253, 0.88, 1.1706, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    suppressWarnings(mod.combob <- mirt(dataset, model.combo, verbose = FALSE))
    expect_is(mod.combob, 'ConfirmatoryClass')
    expect_equal(mod.combob@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combob, digits=4)))
    expect_equal(cfs, c(1.5244, 0, 0.4483, -1.1219, 0, 1, 0.5352, 0, 0, -1.5111, 0, 1, 1.3662, 0, 0, 1.6436, 0, 1, 1.0016, 0, 0, 0.0085, 0, 1, 0, 1.4241, -0.1975, 3.1172, 2.1099, -0.5173, 0, 0.5498, 0, 2.5916, 1.0638, -0.9273, 0, 1.0643, 0, 2.0327, -0.0166, 0, 1.0047, 0, 1.0242, 0, 1, 0, 0, 1, 0, 1),
                 tolerance = 1e-2)

    fs1 <- fscores(mod1, verbose = FALSE)
    expect_is(fs1, 'matrix')
    fs3 <- fscores(mod.quad, full.scores=TRUE, verbose = FALSE)
    expect_is(fs3, 'matrix')
    fs4 <- fscores(mod.combo, verbose = FALSE)
    expect_is(fs4, 'matrix')

    TI <- plot(mod1)
    expect_is(TI, 'trellis')
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

