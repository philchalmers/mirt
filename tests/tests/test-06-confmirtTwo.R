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
    expect_equal(extract.mirt(mod1, 'df'), 1512)
    cfs <- as.numeric(do.call(c, coef(mod1)))
    expect_equal(cfs, c(1.6592,0,-1.0015,0,1,0.4588,0,-1.4308,0,1,0.9201,0,1.4623,0,1,0.8413,0.4269,0.038,0,1,0,1.4071,2.9023,1.9254,-0.4004,0,0.6263,2.4973,1.0223,-1.0313,0,0.8677,1.998,0.0892,0,0.9166,1.1144,0,1,0,0,1,0.4774,1),
                 tolerance = 1e-2)
    Theta <- expand.grid(-4:4, -4:4)
    info <- testinfo(mod1, Theta, degrees = c(30,40))
    expect_equal(info[1:4], c(0.3018910, 0.3689341, 0.4608809, 0.6010538), tolerance = 1e-4)

    mod1b <- mirt(dataset,model.1, verbose = FALSE)
    expect_is(mod1b, 'SingleGroupClass')
    expect_equal(extract.mirt(mod1b, 'df'), 1512)
    cfs <- as.numeric(do.call(c, coef(mod1b)))
    expect_equal(cfs, c(1.6146,0,-0.9786,0,1,0.455,0,-1.4272,0,1,0.9052,0,1.46,0,1,0.8899,0.4345,0.0463,0,1,0,1.3875,2.8938,1.9211,-0.3892,0,0.6358,2.5058,1.0294,-1.028,0,0.8851,2.0124,0.0958,0,0.9291,1.1242,0,1,0,0,1,0.4679,1),
                 tolerance = 1e-2)

    mod.quad <- mirt(dataset, model.quad, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.quad, 'SingleGroupClass')
    expect_equal(extract.mirt(mod.quad, 'df'), 1510)
    cfs <- as.numeric(do.call(c, coef(mod.quad)))
    expect_equal(cfs, c(0.9216,0.4482,-1.1633,0,1,0.3867,0.0833,-1.4982,0,1,0.7688,0.2348,1.1864,0,1,1.0083,0.0735,-0.0144,0,1,1.1327,0,2.6885,1.771,-0.3681,0.5721,0,2.4748,1.0122,-1.0161,0.7917,0,1.9597,0.0861,0.9188,0,1.1163,0,1,0,1),
                 tolerance = 1e-2)

    mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.combo, 'SingleGroupClass')
    expect_equal(extract.mirt(mod.combo, 'df'), 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo)))
    expect_equal(cfs, c(2.1777,0,1.0919,-1.2019,0,1,0.3737,0,0,-1.413,0,1,0.8558,0,0,1.4316,0,1,0.9581,0,0,0.0316,0,1,0,1.3907,-0.1952,2.9383,1.9537,-0.3665,0,0.6466,0,2.5157,1.0365,-1.0288,0,0.9055,0,2.033,0.1057,0,0.8483,0,1.1054,0,1,0,0,1,0,1),
                 tolerance = 1e-2)

    mod.combob <- mirt(dataset, model.combo, verbose = FALSE)
    expect_is(mod.combob, 'SingleGroupClass')
    expect_equal(extract.mirt(mod.combob, 'df'), 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combob)))
    expect_equal(cfs, c(2.0294,0,0.9383,-1.1511,0,1,0.4035,0,0,-1.418,0,1,0.8781,0,0,1.4422,0,1,0.9923,0,0,0.0352,0,1,0,1.3989,-0.2031,2.9465,1.9576,-0.3648,0,0.6501,0,2.5172,1.0377,-1.0276,0,0.9077,0,2.0328,0.1067,0,0.8806,0,1.1162,0,1,0,0,1,0,1),
                 tolerance = 1e-2)

    fs1 <- fscores(mod1, verbose = FALSE, full.scores=FALSE)
    expect_is(fs1, 'matrix')
    fs2 <- fscores(mod1, method = 'WLE', response.pattern = dataset[1:2,])
    expect_equal(as.vector(fs2[,c('F1', 'F2')]), c(-0.3940940, -0.8773481, -2.1421944,  0.7922160),
                 tolerance=1e-4)
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

