context('mixedmirt')

test_that('mixed dich', {
    set.seed(1234)
    N <- 750
    a <- matrix(rlnorm(10,.2,.5),10,1)
    d <- matrix(rnorm(10), 10)
    Theta <- matrix(sort(rnorm(N)))
    pseudoIQ <- scale(Theta * 5 + 100  + rnorm(N, 0 , 5))
    group <- factor(rep(c('G1','G2','G3'), each = N/3))
    data <- simdata(a,d,N, itemtype = rep('dich',10), Theta=Theta)
    covdata <- data.frame(group, pseudoIQ)
    mixedmirt1 <- 'Theta = 1-10'
    model <- mirt.model(mixedmirt1, quiet = TRUE)

    #simple latent regression
    mod0 <- mirt(data, 1, 'Rasch', covdata=covdata, formula = ~ group + pseudoIQ, verbose=FALSE)
    expect_equal(mod0@logLik, -4058.968, tolerance = 1e-2)
    cfs <- coef(mod0, digits = 10)
    expect_equal(as.numeric(cfs$lr.betas), c(0.0000000, 0.9548921, 1.9383165, 0.1877870), tolerance=1e-4)
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    set.seed(1)
    bs <- boot.mirt(mod0, R = 3)
    expect_is(bs, 'boot')
    fs <- fscores(mod0, full.scores.SE=TRUE)
    expect_equal(as.numeric(head(fs)), c(-0.2821604,-0.5666264,-0.4924516,-0.3397658,
                                         -0.486331,-0.4766829,0.2692515,0.2709036,
                                         0.2704817,0.2695925,0.2704465,0.270391), tolerance=1e-4)

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 10)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.6944,-1.886,-1.5027,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-2.1067,-2.3065,-1.9069,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.7012,-1.893,-1.5094,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.0454,-1.2298,-0.861,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-0.3108,-0.4959,-0.1256,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.2415,-1.4274,-1.0556,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.6807,-1.8721,-1.4893,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-2.1216,-2.3218,-1.9215,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-2.0255,-2.2235,-1.8276,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,1.6168,1.3461,1.8875,0,NA,NA,1,NA,NA,0,NA,NA,0.1057,NaN,NaN),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(ncol(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], 2.262686, tolerance = 1e-4)

    mod1a <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, SE=FALSE,
                      verbose = FALSE, draws = 10, internal_constraints = FALSE)
    cfs <- as.numeric(do.call(c, coef(mod1a, digits=4)))
    expect_equal(cfs, c(0.7872,1.6548,1,-1.3497,0,1,1.3846,2.7506,1,-2.4455,0,1,2.4523,4.3594,1,-3.1411,0,1,0.3555,0.6675,1,-0.2773,0,1,1.5055,2.6693,1,-0.5364,0,1,1.4066,3.1572,1,-1.6433,0,1,0.9762,1.9109,1,-1.5035,0,1,1.4161,2.5814,1,-2.394,0,1,0.8393,1.9342,1,-1.7969,0,1,0.878,2.0502,1,1.7316,0,1,0,0.1445),
                 tolerance = 1e-2)

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 10)
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.0419,-0.268,0.3519,-1.7627,-1.9567,-1.5686,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.2383,-0.0529,0.5295,-2.1902,-2.3965,-1.9838,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.1176,-0.1265,0.3617,-1.772,-1.9667,-1.5773,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,-2.2564,-2.9486,-1.5642,-1.0794,-1.3454,-0.8135,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.0714,-0.1633,0.3061,-0.3678,-0.5514,-0.1842,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.0955,-0.1772,0.3683,-1.3073,-1.4944,-1.1202,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.0462,-0.2475,0.3399,-1.7489,-1.9427,-1.5552,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,0.3284,0.0327,0.6241,-2.2154,-2.4242,-2.0065,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,-0.0514,-0.3201,0.2173,-2.0918,-2.2933,-1.8904,0,NA,NA,1,NA,NA,1.1851,1.0395,1.3306,2.3645,2.2091,2.5198,-0.0592,-0.4572,0.3389,1.5601,1.2901,1.83,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 10, verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1,NA,NA,-0.5298,-0.6795,-0.3801,0,NA,NA,1,NA,NA,1,NA,NA,-0.9425,-1.0999,-0.7851,0,NA,NA,1,NA,NA,1,NA,NA,-0.5366,-0.6864,-0.3868,0,NA,NA,1,NA,NA,1,NA,NA,0.1151,-0.0306,0.2607,0,NA,NA,1,NA,NA,1,NA,NA,0.8483,0.6959,1.0006,0,NA,NA,1,NA,NA,1,NA,NA,-0.0796,-0.2255,0.0663,0,NA,NA,1,NA,NA,1,NA,NA,-0.5161,-0.6657,-0.3666,0,NA,NA,1,NA,NA,1,NA,NA,-0.9575,-1.1152,-0.7998,0,NA,NA,1,NA,NA,1,NA,NA,-0.8609,-1.0165,-0.7054,0,NA,NA,1,NA,NA,1,NA,NA,2.8436,2.5879,3.0994,0,NA,NA,1,NA,NA,0,NA,NA,0.3625,NaN,NaN,0.5556,NaN,NaN),
                 tolerance = 1e-2)
})

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model, SE=FALSE,
                     fixed = ~ 0 + group, verbose = FALSE, draws = 10)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod, digits=4))))
    expect_equal(cfs, c(-0.0436,1,0,1,2,3,0,3.0499,5.6391,4.273,-0.0436,1,0,1,2,3,0,1.8812,2.7964,0.9689,-0.0436,1,0,1,2,3,0,2.6243,4.0441,2.9327,-0.0436,1,0,1,2,3,0,2.4307,3.3324,2.002,0,0.955),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(mod@df - mod2@df, 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2, digits=4))))
    expect_equal(cfs, c(-0.1617,-0.3441,0.0206,0.8287,0.6219,1.0356,0,1,2,3,0,2.8314,1.7696,3.8931,5.3694,4.2385,6.5002,4.1331,2.9835,5.2826,-0.1617,-0.3441,0.0206,0.8375,0.5892,1.0858,0,1,2,3,0,1.7767,1.2525,2.301,2.7154,2.1016,3.3292,1.0623,0.4093,1.7153,-0.1617,-0.3441,0.0206,2.5641,0,1,2,3,0,5.2528,7.6976,5.6962,-0.1617,-0.3441,0.0206,0.6946,0.5047,0.8846,0,1,2,3,0,2.1363,1.5783,2.6943,3.0002,2.3646,3.6358,1.9153,1.2339,2.5967,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3, digits=4))))
    expect_equal(cfs, c(-0.198,-0.5895,0.1935,0.9961,0.6142,1.378,4.9121,3.9526,5.8715,2.7023,2.2476,3.157,-1.3542,-1.7526,-0.9557,-0.198,-0.5895,0.1935,1.2157,0.8752,1.5562,3.0075,2.4889,3.5261,0.9901,0.6291,1.3511,-2.1652,-2.619,-1.7115,-0.198,-0.5895,0.1935,2.52,1.3941,3.6459,5.6312,3.6812,7.5813,2.4423,1.4088,3.4759,-1.9982,-2.6822,-1.3142,-0.198,-0.5895,0.1935,1.0548,0.6935,1.4161,3.4081,2.8626,3.9535,1.0723,0.7388,1.4057,-1.575,-1.9876,-1.1624,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1, digits=4))))
    expect_equal(cfs, c(0.9974,0.5716,1.4232,4.8226,3.8584,5.7869,2.6142,2.1828,3.0456,-1.4414,-1.7829,-1.0998,1.238,0.871,1.6049,2.938,2.4713,3.4046,0.9094,0.628,1.1907,-2.2668,-2.6928,-1.8409,2.5964,0.5357,4.6572,5.6735,2.5245,8.8225,2.418,0.7886,4.0475,-2.1134,-3.1905,-1.0363,1.0398,0.5392,1.5404,3.3063,2.7464,3.8661,0.982,0.7228,1.2412,-1.654,-2.0642,-1.2438,0,1,1e-04),
                 tolerance = 1e-2)

    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)

    ## latent regression
    set.seed(1234)
    n <- 250
    Theta <- matrix(c(rnorm(n, -1, sd = sqrt(1/3)),
                      rnorm(n,0, sd = sqrt(1/3)),
                      rnorm(n, 1, sd = sqrt(1/3))))
    dat <- simdata(matrix(rlnorm(10)), matrix(rnorm(10)), N=n*3, Theta=Theta, itemtype = 'dich')
    covdata <- data.frame(group=rep(c('g1', 'g2', 'g3'), each=n))

    mod <- mixedmirt(dat, covdata = covdata, 1, itemtype = '2PL', fixed = ~ 0 + items,
                     lr.fixed = ~ group, verbose=FALSE)
    cfs <- coef(mod, digits=10, printSE=TRUE)
    expect_equal(as.numeric(cfs$lr.betas)[-c(1:2)], c(1.36699167, 0.09689971, 2.86660204, 0.11157517),
                 tolerance=1e-4)
    expect_equal(mod@logLik, -4360.583, tolerance = 1e-2)
    expect_equal(mod@df, 1001)

})
