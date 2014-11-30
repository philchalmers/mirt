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
    expect_equal(plaus[[1]][1:4], c(-0.6071646, -0.4914698, -0.1991302, -0.9721484),
                 tolerance = 1e-4)
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
    expect_equal(cfs, c(1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.704,-1.8948,-1.5132,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-2.1167,-2.3158,-1.9176,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.7108,-1.9018,-1.5199,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.0544,-1.238,-0.8709,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-0.319,-0.5036,-0.1344,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.2508,-1.4358,-1.0657,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-1.6903,-1.8809,-1.4997,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-2.1316,-2.3311,-1.9322,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,-2.0354,-2.2327,-1.8382,0,NA,NA,1,NA,NA,1.1102,0.9689,1.2515,2.2606,2.1118,2.4094,1,NA,NA,1.6107,1.3396,1.8818,0,NA,NA,1,NA,NA,0,NA,NA,0.1111,NaN,NaN),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(ncol(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], 2.3378, tolerance = 1e-4)

    mod1a <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, SE=FALSE,
                      verbose = FALSE, draws = 10, internal_constraints = FALSE)
    cfs <- as.numeric(do.call(c, coef(mod1a, digits=4)))
    expect_equal(cfs, c(0.7945,1.6548,1,-1.3507,0,1,1.3901,2.7467,1,-2.4426,0,1,2.4551,4.3481,1,-3.137,0,1,0.3674,0.6741,1,-0.2854,0,1,1.5117,2.664,1,-0.5424,0,1,1.4113,3.148,1,-1.643,0,1,0.9827,1.9094,1,-1.5037,0,1,1.4213,2.5786,1,-2.3912,0,1,0.8462,1.9336,1,-1.7959,0,1,0.9109,2.0724,1,1.6948,0,1,0,0.1185),
                 tolerance = 1e-2)

    mod_items <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items,
                       verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items)[['GroupPars']], coef(mod_items)[['items']])
    expect_equal(cfs[1:3], c(0, 1.028, 1.094), tolerance = 1e-3)

    mod_items.group <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items:group,
                           verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items.group)[['GroupPars']], coef(mod_items.group)[['items:group']])
    expect_equal(cfs[1:3], c(0, .127, 2.166), tolerance = 1e-3)

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 10)
    expect_is(mod1b, 'MixedClass')
    expect_equal(mod1b@df, 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,0.0077,-0.3089,0.3243,-1.7898,-1.9915,-1.5881,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,0.1293,-0.1766,0.4353,-2.2116,-2.4254,-1.9978,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,0.0435,-0.181,0.2679,-1.7976,-1.9988,-1.5964,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,-3.0561,-5.5602,-0.5519,-1.0574,-1.4042,-0.7106,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,0.0312,-0.2424,0.3048,-0.388,-0.5755,-0.2004,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,0.0627,-0.2234,0.3487,-1.3321,-1.5253,-1.139,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,-0.0471,-0.3544,0.2602,-1.7749,-1.9759,-1.5739,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,0.2739,-0.0269,0.5746,-2.2391,-2.4521,-2.026,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,-0.1397,-0.5281,0.2487,-2.1214,-2.3301,-1.9127,0,NA,NA,1,NA,NA,1.2094,1.0533,1.3656,2.4111,2.2399,2.5822,-0.1228,-0.5219,0.2763,1.5526,1.2816,1.8237,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 10, verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs[124:129], c(0.0581,0.0466,0.0696,1.0838,0.6372,1.5305), tolerance = 1e-2)
})

test_that('polytomous', {
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model, SE=FALSE,
                     fixed = ~ 0 + group, verbose = FALSE, draws = 10)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod, digits=4))))
    expect_equal(cfs, c(-0.0341,1,0,1,2,3,0,3.0959,5.707,4.3338,-0.0341,1,0,1,2,3,0,1.9118,2.8388,0.9956,-0.0341,1,0,1,2,3,0,2.6607,4.0983,2.9785,-0.0341,1,0,1,2,3,0,2.4637,3.3796,2.0369,0,1.0008),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(mod@df - mod2@df, 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2, digits=4))))
    expect_equal(cfs, c(-0.1584,-0.3861,0.0693,0.8302,0.5955,1.0649,0,1,2,3,0,2.8404,1.7808,3.8999,5.3862,4.2391,6.5334,4.1557,2.9773,5.334,-0.1584,-0.3861,0.0693,0.8382,0.6411,1.0352,0,1,2,3,0,1.7847,1.271,2.2984,2.7302,2.1151,3.3453,1.0829,0.4113,1.7546,-0.1584,-0.3861,0.0693,2.5422,0,1,2,3,0,5.2397,7.6962,5.7294,-0.1584,-0.3861,0.0693,0.6934,0.4951,0.8918,0,1,2,3,0,2.1404,1.587,2.6938,3.0092,2.3793,3.6392,1.93,1.2478,2.6122,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 10,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(mod3@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3, digits=4))))
    expect_equal(cfs, c(-0.2067,-0.5178,0.1045,1.0103,0.5874,1.4332,4.9372,3.9357,5.9387,2.7222,2.2317,3.2127,-1.3536,-1.7094,-0.9977,-0.2067,-0.5178,0.1045,1.2074,0.8815,1.5334,3.0128,2.534,3.4916,0.9974,0.6821,1.3127,-2.1541,-2.5639,-1.7442,-0.2067,-0.5178,0.1045,2.4936,0.8239,4.1633,5.6153,3.4068,7.8238,2.4413,1.3762,3.5063,-1.976,-2.973,-0.9789,-0.2067,-0.5178,0.1045,1.0676,0.6611,1.4741,3.4295,2.8326,4.0264,1.0854,0.7645,1.4062,-1.5752,-1.9492,-1.2011,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(rmod1@df, 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1, digits=4))))
    expect_equal(cfs, c(1.0349,0.6996,1.3702,4.8749,3.9217,5.8281,2.6497,2.2231,3.0764,-1.4491,-1.7604,-1.1379,1.2446,0.8882,1.601,2.955,2.4756,3.4344,0.9207,0.6299,1.2116,-2.2647,-2.6724,-1.8571,2.3348,1.6211,3.0486,5.3065,4.0966,6.5164,2.2587,1.6466,2.8708,-1.9636,-2.502,-1.4252,1.0892,0.7292,1.4492,3.3603,2.8233,3.8974,1.0041,0.7246,1.2837,-1.6724,-2.0159,-1.3289,0,1,0.0044),
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
    expect_equal(as.numeric(cfs$lr.betas)[-c(1:2)], c(1.4272120, 0.2088452, 2.9126819, 0.1428500),
                 tolerance=1e-4)
    expect_equal(mod@logLik, -4359.447, tolerance = 1e-2)
    expect_equal(mod@df, 1001)

    #uncorrelated random slope
    covdata$theta <- Theta
    mod <- mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items,
                     random = ~ -1 + theta|group, verbose=FALSE)
    cfs <- coef(mod, digits=10)
    expect_equal(as.numeric(cfs$group[1,]), c(0.06698951, 0.00000000, 0.42557580), tolerance = 1e-4)

})
