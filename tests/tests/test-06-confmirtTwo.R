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

    mod1 <- mirt(dataset,model.1, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod1, 'ConfirmatoryClass')
    expect_equal(mod1@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.4404, 0.8408, 2.0399, 0, NA, NA, -0.928, -1.0857, -0.7702, 0, NA, NA, 1, NA, NA, 0.5379, 0.3835, 0.6923, 0, NA, NA, -1.4507, -1.5751, -1.3264, 0, NA, NA, 1, NA, NA, 0.7851, 0.64, 0.9303, 0, NA, NA, 1.3396, 1.1932, 1.486, 0, NA, NA, 1, NA, NA, 1.3087, 0.0061, 2.6113, 0.5335, 0.1948, 0.8722, 0.0461, 0.0105, 0.0818, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.4021, 1.1905, 1.6137, 2.8563, 2.6652, 3.0474, 1.8527, 1.7067, 1.9987, -0.4618, -0.6038, -0.3198, 0, NA, NA, 0.4311, 0.3212, 0.5409, 2.5296, 2.3682, 2.6909, 0.9916, 0.9033, 1.0799, -1.0149, -1.0999, -0.93, 0, NA, NA, 0.8428, 0.5195, 1.1662, 1.8889, 1.6365, 2.1413, 0.079, 0.062, 0.096, 0, NA, NA, 0.9901, 0.7501, 1.23, 1.0208, 0.921, 1.1206, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0.44, 0.3234, 0.5566, 1, NA, NA),
                 tolerance = 1e-2)

    mod1b <- mirt(dataset,model.1, verbose = FALSE)
    expect_is(mod1b, 'ConfirmatoryClass')
    expect_equal(mod1b@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.5526, 0, -0.9582, 0, 1, 0.5616, 0, -1.4569, 0, 1, 0.7901, 0, 1.3428, 0, 1, 1.1662, 0.6007, 0.0474, 0, 1, 0, 1.4272, 2.8804, 1.8688, -0.4648, 0, 0.4275, 2.5285, 0.991, -1.0141, 0, 0.8374, 1.8871, 0.0792, 0, 0.9758, 1.0169, 0, 1, 0, 0, 1, 0.4246, 1),
                 tolerance = 1e-2)

    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 1510)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.871, 0.7172, 1.0248, 0.2957, 0.1718, 0.4196, -1.0208, -1.1598, -0.8818, 0, NA, NA, 1, NA, NA, 0.4826, 0.3796, 0.5856, -0.0556, -0.3647, 0.2534, -1.3805, -1.6113, -1.1496, 0, NA, NA, 1, NA, NA, 0.7484, 0.5265, 0.9704, 0.2752, -0.0324, 0.5828, 1.0777, 0.8278, 1.3276, 0, NA, NA, 1, NA, NA, 1.5858, 1.3878, 1.7839, 0.4895, 0.3593, 0.6197, -0.2305, -0.3194, -0.1417, 0, NA, NA, 1, NA, NA, 1.0549, 0.4532, 1.6567, 0, NA, NA, 2.5958, 2.1809, 3.0107, 1.6652, 1.3772, 1.9532, -0.4206, -0.5231, -0.3181, 0.3504, 0.2936, 0.4073, 0, NA, NA, 2.5036, 2.3433, 2.6639, 0.9778, 0.8808, 1.0748, -1.0013, -1.1031, -0.8995, 0.7009, 0.5795, 0.8223, 0, NA, NA, 1.8259, 1.7006, 1.9513, 0.0738, -0.0159, 0.1634, 0.9391, 0.5567, 1.3216, 0, NA, NA, 1.0044, 0.9215, 1.0873, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.7616, 1.3354, 2.1878, 0, NA, NA, 0.4541, 0.2335, 0.6747, -1.0423, -1.189, -0.8956, 0, NA, NA, 1, NA, NA, 0.533, 0.3471, 0.7188, 0, NA, NA, 0, NA, NA, -1.4486, -1.5776, -1.3197, 0, NA, NA, 1, NA, NA, 0.7591, 0.6681, 0.8501, 0, NA, NA, 0, NA, NA, 1.3316, 1.2102, 1.453, 0, NA, NA, 1, NA, NA, 1.2975, 1.0146, 1.5804, 0, NA, NA, 0, NA, NA, 0.0434, -0.0693, 0.1561, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1.4591, 1.1663, 1.7519, -0.2748, -0.4434, -0.1063, 2.9675, 2.7869, 3.1482, 1.9248, 1.8079, 2.0418, -0.44, -0.5637, -0.3162, 0, NA, NA, 0.4784, 0.3611, 0.5956, 0, NA, NA, 2.5502, 2.3795, 2.7209, 1.0036, 0.8983, 1.1089, -1.0209, -1.1248, -0.9169, 0, NA, NA, 0.8376, 0.6631, 1.012, 0, NA, NA, 1.892, 1.7321, 2.0519, 0.0855, -0.0178, 0.1888, 0, NA, NA, 0.905, 0.8179, 0.9922, 0, NA, NA, 1.0025, 0.8866, 1.1183, 0, NA, NA, 1, NA, NA, 0, NA, NA, 0, NA, NA, 1, NA, NA, 0, NA, NA, 1, NA, NA),
                 tolerance = 1e-2)

    mod.combob <- mirt(dataset, model.combo, verbose = FALSE)
    expect_is(mod.combob, 'ConfirmatoryClass')
    expect_equal(mod.combob@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combob, digits=4)))
    expect_equal(cfs, c(1.6765, 0, 0.4179, -1.0214, 0, 1, 0.542, 0, 0, -1.4527, 0, 1, 0.7514, 0, 0, 1.327, 0, 1, 1.3204, 0, 0, 0.0407, 0, 1, 0, 1.4954, -0.2818, 3.0102, 1.9579, -0.4358, 0, 0.4644, 0, 2.5463, 1.0025, -1.0165, 0, 0.837, 0, 1.8952, 0.089, 0, 0.892, 0, 1.0029, 0, 1, 0, 0, 1, 0, 1),
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

