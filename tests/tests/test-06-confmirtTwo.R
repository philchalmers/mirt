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
    expect_equal(cfs, c(1.5534,0,-0.9591,0,1,0.5688,0,-1.4589,0,1,0.7977,0,1.3456,0,1,1.1897,0.5863,0.0463,0,1,0,1.3994,2.8545,1.85,-0.4642,0,0.4358,2.5308,0.9919,-1.0165,0,0.8481,1.8906,0.0775,0,0.9863,1.0179,0,1,0,0,1,0.4256,1),
                 tolerance = 1e-2)

    mod1b <- mirt(dataset,model.1, verbose = FALSE)
    expect_is(mod1b, 'ConfirmatoryClass')
    expect_equal(mod1b@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.5493,0,-0.9572,0,1,0.5616,0,-1.4568,0,1,0.7901,0,1.3428,0,1,1.1682,0.5995,0.0475,0,1,0,1.4268,2.88,1.8685,-0.4648,0,0.4276,2.5285,0.9911,-1.0141,0,0.8374,1.8871,0.0792,0,0.9759,1.017,0,1,0,0,1,0.4253,1),
                 tolerance = 1e-2)

    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'ConfirmatoryClass')
    expect_equal(mod.quad@df, 1510)
    cfs <- as.numeric(do.call(c, coef(mod.quad, digits=4)))
    expect_equal(cfs, c(0.8923,0.1328,-0.8929,0,1,0.5117,-0.088,-1.3597,0,1,0.7162,0.1825,1.1458,0,1,1.5009,0.3353,-0.1473,0,1,1.0598,0,2.5935,1.6636,-0.4238,0.3515,0,2.5026,0.9769,-1.0026,0.7102,0,1.8271,0.0722,0.9147,0,0.9944,0,1,0,1),
                 tolerance = 1e-2)

    mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.combo, 'ConfirmatoryClass')
    expect_equal(mod.combo@df, 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo, digits=4)))
    expect_equal(cfs, c(1.6228,0,0.3441,-1.0014,0,1,0.5473,0,0,-1.4535,0,1,0.7777,0,0,1.3372,0,1,1.339,0,0,0.0431,0,1,0,1.4432,-0.2535,2.9491,1.9152,-0.437,0,0.4848,0,2.5517,1.0041,-1.0227,0,0.8448,0,1.8946,0.0856,0,0.9293,0,1.0089,0,1,0,0,1,0,1),
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

