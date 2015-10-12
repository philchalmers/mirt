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
    expect_is(mod1, 'SingleGroupClass')
    expect_equal(mod1@Fit$df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.5565,0,-0.9627,0,1,0.5455,0,-1.4531,0,1,0.7978,0,1.3458,0,1,1.1692,0.5865,0.0477,0,1,0,1.4023,2.8587,1.8546,-0.4592,0,0.4415,2.5338,0.9943,-1.0163,0,0.8477,1.8935,0.0816,0,0.9704,1.0164,0,1,0,0,1,0.4309,1),
                 tolerance = 1e-2)

    mod1b <- mirt(dataset,model.1, verbose = FALSE)
    expect_is(mod1b, 'SingleGroupClass')
    expect_equal(mod1b@Fit$df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.5493,0,-0.9572,0,1,0.5616,0,-1.4568,0,1,0.7901,0,1.3428,0,1,1.1682,0.5995,0.0475,0,1,0,1.4268,2.88,1.8685,-0.4648,0,0.4276,2.5285,0.9911,-1.0141,0,0.8374,1.8871,0.0792,0,0.9759,1.017,0,1,0,0,1,0.4253,1),
                 tolerance = 1e-2)

    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'SingleGroupClass')
    expect_equal(mod.quad@Fit$df, 1510)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.8852,0.1199,-0.8806,0,1,0.499,-0.0729,-1.3696,0,1,0.7136,0.1747,1.1533,0,1,1.4444,0.2611,-0.1055,0,1,1.055,0,2.5917,1.6627,-0.4225,0.3559,0,2.5045,0.9777,-1.0034,0.7043,0,1.8252,0.0726,0.9228,0,0.9983,0,1,0,1),
                 tolerance = 1e-2)

    mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.combo, 'SingleGroupClass')
    expect_equal(mod.combo@Fit$df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.5661,0,0.346,-0.9944,0,1,0.5206,0,0,-1.4488,0,1,0.7675,0,0,1.3309,0,1,1.3356,0,0,0.0361,0,1,0,1.454,-0.2666,2.9673,1.9292,-0.431,0,0.481,0,2.5526,1.0059,-1.0202,0,0.8513,0,1.902,0.0896,0,0.9081,0,1.0071,0,1,0,0,1,0,1),
                 tolerance = 1e-2)

    mod.combob <- mirt(dataset, model.combo, verbose = FALSE)
    expect_is(mod.combob, 'SingleGroupClass')
    expect_equal(mod.combob@Fit$df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combob, digits=4)))
    expect_equal(cfs, c(1.6765, 0, 0.4179, -1.0214, 0, 1, 0.542, 0, 0, -1.4527, 0, 1, 0.7514, 0, 0, 1.327, 0, 1, 1.3204, 0, 0, 0.0407, 0, 1, 0, 1.4954, -0.2818, 3.0102, 1.9579, -0.4358, 0, 0.4644, 0, 2.5463, 1.0025, -1.0165, 0, 0.837, 0, 1.8952, 0.089, 0, 0.892, 0, 1.0029, 0, 1, 0, 0, 1, 0, 1),
                 tolerance = 1e-2)

    fs1 <- fscores(mod1, verbose = FALSE, full.scores=FALSE)
    expect_is(fs1, 'matrix')
    fs3 <- fscores(mod.quad, full.scores=TRUE, verbose = FALSE)
    expect_is(fs3, 'matrix')
    fs4 <- fscores(mod.combo, verbose = FALSE, full.scores=FALSE)
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

