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
    expect_equal(mod0@Fit$logLik, -4058.968, tolerance = 1e-2)
    cfs <- coef(mod0, digits = 10)
    expect_equal(as.numeric(cfs$lr.betas), c(0.0000000, 0.9548921, 1.9383165, 0.1877870), tolerance=1e-4)
    set.seed(1234)
    plaus <- fscores(mod0, plausible.draws = 2)
    expect_equal(plaus[[1]][1:4], c(-0.6070500, -0.2729041, -0.3763773, -0.4946551),
                 tolerance = 1e-4)
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    set.seed(1)
    bs <- boot.mirt(mod0, R = 3)
    expect_is(bs, 'boot')
    fs <- fscores(mod0, full.scores.SE=TRUE, full.scores=FALSE)
    expect_equal(as.numeric(head(fs)), c(-0.2820982,-0.5666356,-0.4924289,-0.3397344,-0.4862802,-0.4766593,0.2692081,0.2708598,0.2704379,0.2695491,0.2704026,0.2703473), tolerance=1e-4)

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 1)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1, digits=4)))
    expect_equal(cfs, c(1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-1.6979,-1.8905,-1.5053,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-2.1101,-2.3119,-1.9083,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-1.7048,-1.8975,-1.512,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-1.0493,-1.2334,-0.8652,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-0.3151,-0.4994,-0.1309,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-1.2453,-1.4313,-1.0594,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-1.6843,-1.8766,-1.4919,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-2.125,-2.3272,-1.9229,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,-2.029,-2.2287,-1.8292,0,NA,NA,1,NA,NA,1.111,0.9628,1.2593,2.2422,2.0891,2.3954,1,NA,NA,1.6114,1.3404,1.8824,0,NA,NA,1,NA,NA,0,NA,NA,0.1035,0.0541,0.1529),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(length(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], 2.156415, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod1, R=2)
    expect_is(bs, 'boot')

    mod1a <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, SE=FALSE,
                      verbose = FALSE, draws = 1, internal_constraints = FALSE)
    cfs <- as.numeric(do.call(c, coef(mod1a, digits=4)))
    expect_equal(cfs, c(0.7882,1.6403,1,-1.3455,0,1,1.3836,2.731,1,-2.4362,0,1,2.4475,4.3307,1,-3.1301,0,1,0.3616,0.6616,1,-0.2824,0,1,1.5031,2.649,1,-0.5388,0,1,1.4038,3.1316,1,-1.6373,0,1,0.9761,1.8945,1,-1.4983,0,1,1.4147,2.563,1,-2.3848,0,1,0.8402,1.9187,1,-1.7901,0,1,0.9049,2.062,1,1.6933,0,1,0,0.1088),
                 tolerance = 1e-2)

    mod_items <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items,
                       verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items)[['GroupPars']], coef(mod_items)[['items']])
    expect_equal(cfs[1:3], c(0, 1.005, 1.143), tolerance = 1e-3)

    mod_items.group <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items:group,
                           verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items.group)[['GroupPars']], coef(mod_items.group)[['items:group']])
    expect_equal(cfs[1:3], c(0, .115, 2.254), tolerance = 1e-3)
    set.seed(1)
    bs <- boot.mirt(mod_items.group, R=2)
    expect_is(bs, 'boot')

    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 1)
    expect_is(mod1b, 'MixedClass')
    expect_equal(extract.mirt(mod1b, 'df'), 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b, digits=4)))
    expect_equal(cfs, c(1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,-6e-04,-0.2538,0.2526,-1.7908,-1.9861,-1.5955,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,0.133,-0.1237,0.3896,-2.2139,-2.4201,-2.0077,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,0.056,-0.1629,0.2749,-1.7997,-1.9955,-1.6039,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,-3.132,-4.2471,-2.0169,-1.05,-1.4039,-0.6962,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,0.0428,-0.1931,0.2787,-0.3885,-0.5725,-0.2045,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,0.0864,-0.1338,0.3067,-1.3339,-1.5223,-1.1454,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,-0.0749,-0.3338,0.1841,-1.7758,-1.9707,-1.5809,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,0.2577,1e-04,0.5153,-2.2396,-2.4487,-2.0305,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,-0.1304,-0.3719,0.1111,-2.1229,-2.3255,-1.9204,0,NA,NA,1,NA,NA,1.2111,1.0631,1.3591,2.4128,2.2574,2.5681,-0.1141,-0.5071,0.279,1.5521,1.2816,1.8226,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 1, verbose = FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(extract.mirt(rmod1, 'df'), 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1, digits=4)))
    expect_equal(cfs[124:129], c(0.0536, 0.0397, 0.0676, 1.0902, 0.6417, 1.5387), tolerance = 1e-2)

    #polytomous
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model, SE=FALSE,
                     fixed = ~ 0 + group, verbose = FALSE, draws = 1)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod, digits=4))))
    expect_equal(cfs, c(-0.0515,1,0,1,2,3,0,3.0327,5.6104,4.2493,-0.0515,1,0,1,2,3,0,1.8659,2.7745,0.958,-0.0515,1,0,1,2,3,0,2.6075,4.0176,2.9117,-0.0515,1,0,1,2,3,0,2.4149,3.3088,1.9863,0,0.9265),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(extract.mirt(mod, 'df') - extract.mirt(mod2, 'df'), 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2, digits=4))))
    expect_equal(cfs, c(-0.1658,-0.3961,0.0645,0.8385,0.4843,1.1928,0,1,2,3,0,2.8555,1.7005,4.0105,5.4041,4.1164,6.6918,4.1651,2.9216,5.4087,-0.1658,-0.3961,0.0645,0.822,0.5568,1.0872,0,1,2,3,0,1.7642,1.2394,2.2889,2.7013,2.086,3.3166,1.0618,0.396,1.7277,-0.1658,-0.3961,0.0645,2.398,0.9417,3.8543,0,1,2,3,0,4.9963,2.4428,7.5497,7.3389,3.7964,10.8813,5.4399,2.6573,8.2226,-0.1658,-0.3961,0.0645,0.7278,0.4345,1.0211,0,1,2,3,0,2.1803,1.5836,2.777,3.057,2.3766,3.7374,1.9533,1.2594,2.6473,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(extract.mirt(mod3, 'df'), 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3, digits=4))))
    expect_equal(cfs, c(-0.2005,-0.5234,0.1223,1.0151,0.6339,1.3964,4.937,3.968,5.9061,2.7168,2.2619,3.1717,-1.3634,-1.713,-1.0138,-0.2005,-0.5234,0.1223,1.1613,0.8347,1.4878,2.969,2.4977,3.4403,0.9771,0.666,1.2882,-2.1311,-2.5415,-1.7206,-0.2005,-0.5234,0.1223,2.5549,1.2906,3.8192,5.695,3.7818,7.6083,2.4632,1.5444,3.382,-2.0238,-2.8748,-1.1729,-0.2005,-0.5234,0.1223,1.0885,0.719,1.4579,3.4395,2.8808,3.9982,1.0818,0.7716,1.392,-1.5923,-1.9637,-1.2208,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE, SE=FALSE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(extract.mirt(rmod1, 'df'), 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1, digits=4))))
    expect_equal(cfs, c(1.0388,4.8464,2.6252,-1.4723,1.2431,2.9268,0.8951,-2.2858,2.2871,5.1966,2.1899,-1.9791,1.0995,3.3382,0.9803,-1.6994,0,1,0.0016),
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
    expect_equal(as.numeric(cfs$lr.betas)[-c(1:2)], c(1.4712503, 0.2312698, 2.9707061, 0.2363475),
                 tolerance=1e-4)
    expect_equal(extract.mirt(mod, 'logLik'), -4360.41, tolerance = 1e-2)
    expect_equal(extract.mirt(mod, 'df'), 1001)

    mod2 <- mixedmirt(dat, covdata = covdata, 1, itemtype = 'Rasch', fixed = ~ 0 + items,
                     lr.fixed = ~ group, verbose=FALSE, draws=1, SE=FALSE)
    so <- summary(mod2, verbose=FALSE)
    expect_equal(as.numeric(c(so$random$Theta, so$lr.out[,1])),
                 c(0.1961649, 0.0000000, 0.7167553, 1.5104188), tolerance=1e-4)
    expect_equal(extract.mirt(mod2, 'logLik'), -4688.784, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod2, R = 3)
    expect_is(bs, 'boot')

    #uncorrelated random slope
    covdata$theta <- Theta
    covdata$cut <- factor(cut(Theta, breaks=2))
    mod <- mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items, SE=FALSE,
                     random = ~ -1 + cut|group, verbose=FALSE, draws=1)
    so <- summary(mod, verbose=FALSE)
    expect_equal(as.numeric(diag(so$random$group)), c(0.6211201, 0.6262019), tolerance = 1e-4)

    mod <- mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items, SE=FALSE,
                     random = ~ -1 + theta|group, verbose=FALSE, draws=1)
    so <- summary(mod, verbose=FALSE)
    expect_equal(as.numeric(diag(so$random$group)), c(0.05447179, 0.51292915), tolerance = 1e-4)

})
