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
    expect_equal(cfs, c(1.4404,0,-0.928,0,1,0.5379,0,-1.4507,0,1,0.7851,0,1.3396,0,1,1.3087,0.5335,0.0461,0,1,0,1.4021,2.8563,1.8527,-0.4618,0,0.4311,2.5296,0.9916,-1.0149,0,0.8428,1.8889,0.079,0,0.9901,1.0208,0,1,0,0,1,0.44,1),
                 tolerance = 1e-2)

    mod1b <- mirt(dataset,model.1, verbose = FALSE)
    expect_is(mod1b, 'ConfirmatoryClass')
    expect_equal(mod1b@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.5517,0,-0.9579,0,1,0.5616,0,-1.4568,0,1,0.7901,0,1.3428,0,1,1.1667,0.6004,0.0475,0,1,0,1.427,2.8803,1.8686,-0.4648,0,0.4275,2.5285,0.9911,-1.0142,0,0.8374,1.8871,0.0792,0,0.9759,1.0169,0,1,0,0,1,0.4248,1),
                 tolerance = 1e-2)

    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 1510)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.871,0.2957,-1.0208,0,1,0.4826,-0.0556,-1.3805,0,1,0.7484,0.2752,1.0777,0,1,1.5858,0.4895,-0.2305,0,1,1.0549,0,2.5958,1.6652,-0.4206,0.3504,0,2.5036,0.9778,-1.0013,0.7009,0,1.8259,0.0738,0.9391,0,1.0044,0,1,0,1),
                 tolerance = 1e-2)

    mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.7616,0,0.4541,-1.0423,0,1,0.533,0,0,-1.4486,0,1,0.7591,0,0,1.3316,0,1,1.2975,0,0,0.0434,0,1,0,1.4591,-0.2748,2.9675,1.9248,-0.44,0,0.4784,0,2.5502,1.0036,-1.0209,0,0.8376,0,1.892,0.0855,0,0.905,0,1.0025,0,1,0,0,1,0,1),
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

