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
    expect_equal(cfs, c(1.725311,0,-1.003746,0,1,0.4296214,0,-1.418875,0,1,0.8736704,0,1.451735,0,1,0.8736917,0.4298316,0.05164668,0,1,0,1.380402,2.886302,1.919996,-0.3822423,0,0.6367169,2.507395,1.030947,-1.026989,0,0.8973311,2.019848,0.09872695,0,0.923406,1.123934,0,1,0,0,1,0.4793236,1),
                 tolerance = 1e-2)
    Theta <- expand.grid(-4:4, -4:4)
    info <- testinfo(mod1, Theta, degrees = c(30,40))
    expect_equal(info[1:4], c(0.3079868,0.3687883,0.4502187,0.593109), tolerance = 1e-4)

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
    expect_equal(cfs, c(0.9008854,0.4617326,-1.174335,0,1,0.3751464,0.07415233,-1.487041,0,1,0.7719363,0.255426,1.169278,0,1,1.037831,0.1488222,-0.06758965,0,1,1.136647,0,2.685626,1.767618,-0.3708233,0.5774083,0,2.475542,1.012218,-1.018043,0.7840536,0,1.954257,0.0851495,0.9447677,0,1.122107,0,1,0,1),
                 tolerance = 1e-2)

    mod.combo <- mirt(dataset, model.combo, verbose = FALSE, draws = 10, method = 'MHRM')
    expect_is(mod.combo, 'SingleGroupClass')
    expect_equal(extract.mirt(mod.combo, 'df'), 1512)
    cfs <- as.numeric(do.call(c, coef(mod.combo)))
    expect_equal(cfs, c(2.175565,0,1.124398,-1.196377,0,1,0.3772646,0,0,-1.411797,0,1,0.8485753,0,0,1.431696,0,1,0.9730396,0,0,0.03585494,0,1,0,1.414372,-0.2296381,2.96445,1.970084,-0.3634654,0,0.6340928,0,2.509455,1.035055,-1.022543,0,0.9094451,0,2.034747,0.1080139,0,0.877999,0,1.115907,0,1,0,0,1,0,1),
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
    expect_equal(as.vector(fs2[,c('F1', 'F2')]), c(-0.3647738,-0.8722322,-2.168684,0.7796325),
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

