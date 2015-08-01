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
    set.seed(1234)
    plaus <- fscores(mod0, plausible.draws = 2)
    expect_equal(plaus[[1]][1:4], c(-0.6074285, -0.4914628, -0.1991007, -0.9724681),
                 tolerance = 1e-4)
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    set.seed(1)
    bs <- boot.mirt(mod0, R = 3)
    expect_is(bs, 'boot')
    fs <- fscores(mod0, full.scores.SE=TRUE)
    expect_equal(as.numeric(head(fs)), c(-0.2823035,-0.5666475,-0.4925317,-0.3398507,-0.4864681,-0.4767635,0.2693515,0.2710047,0.2705827,0.2696926,0.2705478,0.270492), tolerance=1e-4)

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 1)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.704,-1.8948,-1.5132,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-2.1167,-2.3158,-1.9176,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.7108,-1.9018,-1.5199,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.0544,-1.238,-0.8709,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-0.319,-0.5036,-0.1344,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.2508,-1.4358,-1.0657,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.6903,-1.8809,-1.4997,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-2.1316,-2.3311,-1.9322,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-2.0354,-2.2327,-1.8382,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,1.6107,1.3396,1.8818,0,NA,NA,1,NA,NA,0,NA,NA,0.1111,NaN,NaN),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(ncol(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], 2.087929, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod1, R=2)
    expect_is(bs, 'boot')

    mod1a <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, SE=FALSE,
                      verbose = FALSE, draws = 1, internal_constraints = FALSE)
    cfs <- as.numeric(do.call(c, coef(mod1a, digits=4)))
    expect_equal(cfs, c(0.7945,1.6548,1,-1.3507,0,1,1.3901,2.7467,1,-2.4426,0,1,2.4551,4.3481,1,-3.137,0,1,0.3674,0.6741,1,-0.2854,0,1,1.5117,2.664,1,-0.5424,0,1,1.4113,3.148,1,-1.643,0,1,0.9827,1.9094,1,-1.5037,0,1,1.4213,2.5786,1,-2.3912,0,1,0.8462,1.9336,1,-1.7959,0,1,0.9109,2.0724,1,1.6948,0,1,0,0.1185),
                 tolerance = 1e-2)

    mod_items <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items,
                       verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items)[['GroupPars']], coef(mod_items)[['items']])
    expect_equal(cfs[1:3], c(0, .981, 1.114), tolerance = 1e-3)

    mod_items.group <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items:group,
                           verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items.group)[['GroupPars']], coef(mod_items.group)[['items:group']])
    expect_equal(cfs[1:3], c(0, .121, 2.198), tolerance = 1e-3)
    set.seed(1)
    bs <- boot.mirt(mod_items.group, R=2)
    expect_is(bs, 'boot')

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 1)
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.2094,1.061,1.3579,2.4111,2.2539,2.5682,0.0077,-0.25,0.2654,-1.7898,-1.9856,-1.594,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,0.1293,-0.1173,0.376,-2.2116,-2.4178,-2.0053,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,0.0435,-0.1835,0.2705,-1.7976,-1.9937,-1.6016,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,-3.0561,-4.2477,-1.8644,-1.0574,-1.4011,-0.7137,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,0.0312,-0.2154,0.2777,-0.388,-0.5721,-0.2038,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,0.0627,-0.1564,0.2817,-1.3321,-1.5207,-1.1436,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,-0.0471,-0.2961,0.202,-1.7749,-1.9704,-1.5795,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,0.2739,0.0146,0.5332,-2.2391,-2.448,-2.0301,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,-0.1397,-0.3932,0.1139,-2.1214,-2.3248,-1.9179,0,NA,NA,1,NA,NA,1.2094,1.061,1.3579,2.4111,2.2539,2.5682,-0.1228,-0.5105,0.265,1.5526,1.2819,1.8234,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 1, verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs[124:129], c(0.0581, 0.0181, 0.0982, 1.0808, 0.6386, 1.5231), tolerance = 1e-2)
})

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model, SE=FALSE,
                     fixed = ~ 0 + group, verbose = FALSE, draws = 1)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod, digits=4))))
    expect_equal(cfs, c(-0.0322,1,0,1,2,3,0,3.0385,5.6204,4.2602,-0.0322,1,0,1,2,3,0,1.8703,2.7828,0.9653,-0.0322,1,0,1,2,3,0,2.6123,4.0266,2.9216,-0.0322,1,0,1,2,3,0,2.4194,3.3174,1.9951,0,0.9342),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(mod@df - mod2@df, 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2, digits=4))))
    expect_equal(cfs, c(-0.1641,-0.4105,0.0824,0.8267,0.5379,1.1154,0,1,2,3,0,2.8341,1.7311,3.937,5.3795,4.1792,6.5798,4.151,2.9605,5.3414,-0.1641,-0.4105,0.0824,0.8448,0.5349,1.1546,0,1,2,3,0,1.7939,1.246,2.3418,2.7434,2.1063,3.3804,1.0905,0.42,1.761,-0.1641,-0.4105,0.0824,2.4661,0.5397,4.3925,0,1,2,3,0,5.1082,1.8438,8.3726,7.5164,2.9267,12.106,5.5924,2.0101,9.1747,-0.1641,-0.4105,0.0824,0.7003,0.4439,0.9566,0,1,2,3,0,2.1503,1.571,2.7296,3.0236,2.3635,3.6838,1.9408,1.2451,2.6366,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3, digits=4))))
    expect_equal(cfs, c(-0.2067,-0.5235,0.1102,1.0103,0.6619,1.3586,4.9372,3.9762,5.8982,2.7222,2.2735,3.1709,-1.3536,-1.6933,-1.0138,-0.2067,-0.5235,0.1102,1.2074,0.8647,1.5502,3.0128,2.5332,3.4923,0.9974,0.684,1.3109,-2.1541,-2.5795,-1.7287,-0.2067,-0.5235,0.1102,2.4936,1.3573,3.6299,5.6153,3.8747,7.356,2.4413,1.5736,3.3089,-1.976,-2.7274,-1.2245,-0.2067,-0.5235,0.1102,1.0676,0.6829,1.4523,3.4295,2.8692,3.9898,1.0854,0.7761,1.3947,-1.5752,-1.9484,-1.2019,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE, SE=FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1, digits=4))))
    expect_equal(cfs, c(1.0458,4.8686,2.6404,-1.4691,1.2377,2.9332,0.9035,-2.2763,2.2356,5.1481,2.1737,-1.9378,1.1222,3.3691,0.9976,-1.7041,0,1,0.003),
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
                     lr.fixed = ~ group, verbose=FALSE, SE=TRUE)
    cfs <- coef(mod, digits=10, printSE=TRUE)
    expect_equal(as.numeric(cfs$lr.betas)[-c(1:2)], c(1.4341209, 0.1320500, 2.9319028, 0.1505968),
                 tolerance=1e-4)
    expect_equal(mod@logLik, -4359.836, tolerance = 1e-2)
    expect_equal(mod@df, 1001)

    mod2 <- mixedmirt(dat, covdata = covdata, 1, itemtype = 'Rasch', fixed = ~ 0 + items,
                     lr.fixed = ~ group, verbose=FALSE, draws=1, SE=FALSE)
    so <- summary(mod2, verbose=FALSE)
    expect_equal(as.numeric(c(so$random$Theta, so$lr.out[,1])),
                 c(0.1958198, 0.0000000, 0.7081315, 1.4867710), tolerance=1e-4)
    expect_equal(mod2@logLik, -4685.077, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod2, R = 3)
    expect_is(bs, 'boot')

    #uncorrelated random slope
    covdata$theta <- Theta
    covdata$cut <- cut(Theta, breaks=2)
    mod <- mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items, SE=FALSE,
                     random = ~ -1 + cut + theta|group, verbose=FALSE, draws=1)
    so <- summary(mod, verbose=FALSE)
    expect_equal(as.numeric(diag(so$random$group)), c(0.04038883, 0.10988249, 0.57688499), tolerance = 1e-4)

})
