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

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 10)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.6944,-1.886,-1.5027,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-2.1067,-2.3065,-1.9069,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.7012,-1.893,-1.5094,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.0454,-1.2298,-0.861,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-0.3108,-0.4959,-0.1256,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.2415,-1.4274,-1.0556,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-1.6807,-1.8721,-1.4893,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-2.1216,-2.3218,-1.9215,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,-2.0255,-2.2235,-1.8276,0,NA,NA,1,NA,NA,1.112,0.9634,1.2606,2.2478,2.1006,2.395,1,NA,NA,1.6168,1.3461,1.8875,0,NA,NA,1,NA,NA,0,NA,NA,0.1057,0.0301,0.1814),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(ncol(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], 2.18281, tolerance = 1e-4)

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 10)
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.0419,-0.2248,0.3087,-1.7627,-1.9587,-1.5666,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.2383,-0.0433,0.52,-2.1902,-2.3971,-1.9833,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.1176,-0.1239,0.3591,-1.772,-1.9681,-1.5759,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,-2.2564,-3.0821,-1.4307,-1.0794,-1.3626,-0.7963,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.0714,-0.2124,0.3552,-0.3678,-0.5526,-0.183,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.0955,-0.1553,0.3463,-1.3073,-1.4959,-1.1187,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.0462,-0.3907,0.4832,-1.7489,-1.9437,-1.5541,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,0.3284,-0.0943,0.751,-2.2154,-2.4264,-2.0044,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,-0.0514,-0.3761,0.2733,-2.0918,-2.2948,-1.8889,0,NA,NA,1,NA,NA,1.1851,1.0369,1.3333,2.3645,2.2039,2.525,-0.0592,-0.5277,0.4093,1.5601,1.2897,1.8304,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 10, verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs, c(1,NA,NA,-0.5298,-0.6786,-0.3809,0,NA,NA,1,NA,NA,1,NA,NA,-0.9425,-1.0882,-0.7968,0,NA,NA,1,NA,NA,1,NA,NA,-0.5366,-0.6854,-0.3878,0,NA,NA,1,NA,NA,1,NA,NA,0.1151,-0.0255,0.2556,0,NA,NA,1,NA,NA,1,NA,NA,0.8483,0.7457,0.9508,0,NA,NA,1,NA,NA,1,NA,NA,-0.0796,-0.2247,0.0655,0,NA,NA,1,NA,NA,1,NA,NA,-0.5161,-0.665,-0.3673,0,NA,NA,1,NA,NA,1,NA,NA,-0.9575,-1.1031,-0.812,0,NA,NA,1,NA,NA,1,NA,NA,-0.8609,-1.0076,-0.7142,0,NA,NA,1,NA,NA,1,NA,NA,2.8436,2.7543,2.9329,0,NA,NA,1,NA,NA,0,NA,NA,0.3625,0.2327,0.4923,0.5556,0.1066,1.0046),
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
    expect_equal(cfs, c(-0.1617,-0.4868,0.1633,0.8287,0.3302,1.3273,0,1,2,3,0,2.8314,1.538,4.1247,5.3694,3.8254,6.9134,4.1331,2.7104,5.5558,-0.1617,-0.4868,0.1633,0.8375,0.4665,1.2084,0,1,2,3,0,1.7767,1.1748,2.3787,2.7154,1.997,3.4337,1.0623,0.3876,1.737,-0.1617,-0.4868,0.1633,2.5641,-3.3715,8.4998,0,1,2,3,0,5.2528,-4.3739,14.8796,7.6976,-5.5417,20.937,5.6962,-3.6021,14.9945,-0.1617,-0.4868,0.1633,0.6946,0.4359,0.9534,0,1,2,3,0,2.1363,1.5628,2.7099,3.0002,2.3482,3.6522,1.9153,1.2361,2.5944,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3, digits=4))))
    expect_equal(cfs, c(-0.198,-0.6323,0.2363,0.9961,0.6498,1.3424,4.9121,3.9477,5.8765,2.7023,2.2253,3.1792,-1.3542,-1.7989,-0.9094,-0.198,-0.6323,0.2363,1.2157,0.8604,1.5709,3.0075,2.4386,3.5764,0.9901,0.5704,1.4098,-2.1652,-2.6575,-1.673,-0.198,-0.6323,0.2363,2.52,1.2818,3.7581,5.6312,3.4718,7.7907,2.4423,1.2687,3.616,-1.9982,-2.706,-1.2904,-0.198,-0.6323,0.2363,1.0548,0.7337,1.3759,3.4081,2.8427,3.9734,1.0723,0.6923,1.4523,-1.575,-2.0336,-1.1164,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1, digits=4))))
    expect_equal(cfs, c(0.9974,0.4726,1.5222,4.8226,4.0091,5.6361,2.6142,2.3219,2.9065,-1.4414,-1.5174,-1.3653,1.238,0.8545,1.6214,2.938,2.4727,3.4032,0.9094,0.6518,1.167,-2.2668,-2.6777,-1.856,2.5964,-0.3535,5.5463,5.6735,1.0999,10.2471,2.418,0.2389,4.5971,-2.1134,-3.6756,-0.5512,1.0398,0.4362,1.6434,3.3063,3.0207,3.5918,0.982,0.7405,1.2235,-1.654,-1.8281,-1.4798,0,1,1e-04,0,3e-04),
                 tolerance = 1e-2)

    re <- randef(rmod1, ndraws=100)
    expect_is(re, 'list')
    expect_equal(length(re), 2)

})
