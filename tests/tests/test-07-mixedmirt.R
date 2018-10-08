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
    cfs <- coef(mod0)
    expect_equal(as.numeric(cfs$lr.betas), c(0.0000000, 0.8916977, 1.9757340, 0.2168971), tolerance=1e-4)
    set.seed(1234)
    plaus <- fscores(mod0, plausible.draws = 2)
    expect_equal(plaus[[1]][1:4], c(-0.6875802, -0.3252686, -0.3528408, -0.4759154),
                 tolerance = 1e-4)
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    set.seed(1)
    bs <- boot.mirt(mod0, R = 3)
    expect_is(bs, 'boot')
    fs <- fscores(mod0, full.scores.SE=TRUE, full.scores=FALSE)
    expect_equal(as.numeric(head(fs)), c(-0.3326607,-0.6466616,-0.4794936,-0.3070341,-0.4769525,-0.5488866,0.2940349,0.2963674,0.2951423,0.2938395,0.2951234,0.2956562),
                 tolerance=1e-4)

    #group as a fixed effect predictor (aka, uniform dif)
    mod1 <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                       verbose = FALSE, draws = 1)
    expect_is(mod1, 'MixedClass')
    cfs <- as.numeric(do.call(c, coef(mod1)))
    expect_equal(cfs, c(1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-1.7432,-1.9407,-1.5458,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-2.0608,-2.2651,-1.8565,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-1.7014,-1.898,-1.5047,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-1.0408,-1.2289,-0.8527,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-0.2658,-0.4541,-0.0775,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-1.3741,-1.5656,-1.1825,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-1.6805,-1.8768,-1.4842,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-2.121,-2.3268,-1.9152,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,-2.0015,-2.2043,-1.7986,0,NA,NA,1,NA,NA,1.075,0.9182,1.2318,2.3263,2.1645,2.488,1,NA,NA,1.5431,1.2778,1.8084,0,NA,NA,1,NA,NA,0,NA,NA,0.1381,0.1015,0.1747),
                 tolerance = 1e-2)
    names <- wald(mod1)
    L <- matrix(c(1, numeric(length(names) - 1L)), 1L)
    wld <- wald(mod1, L, C=as.numeric(L))
    expect_equal(wld$W[1], .8784299, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod1, R=2)
    expect_is(bs, 'boot')

    mod1a <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group, SE=FALSE,
                      verbose = FALSE, draws = 1, internal_constraints = FALSE)
    cfs <- as.numeric(do.call(c, coef(mod1a)))
    expect_equal(cfs, c(0.8023,1.5944,1,-1.3592,0,1,1.3988,2.9293,1,-2.4552,0,1,2.4305,4.4838,1,-3.1509,0,1,0.307,0.7794,1,-0.2878,0,1,1.4845,3.0533,1,-0.5464,0,1,1.2012,3.0349,1,-1.6528,0,1,0.9915,1.9697,1,-1.513,0,1,1.4304,2.6556,1,-2.4036,0,1,0.8535,2.0584,1,-1.8064,0,1,0.8374,2.3511,1,1.6148,0,1,0,0.1489),
                 tolerance = 1e-2)

    mod_items <- mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items,
                       verbose = FALSE, draws = 1)
    cfs <- c(coef(mod_items)[['GroupPars']], coef(mod_items)[['items']])
    expect_equal(cfs[1:3], c(0.000, 1.083, 1.125), tolerance = 1e-2)

    mod_items.group <- try(mixedmirt(data, covdata, model, fixed = ~ 1, SE=FALSE, random = ~ 1|items:group,
                           verbose = FALSE, draws = 1), TRUE)
    cfs <- c(coef(mod_items.group)[['GroupPars']], coef(mod_items.group)[['items:group']])
    expect_equal(cfs[1:3], c(0.000, 0.1431, 2.2536), tolerance = 1e-3)
    set.seed(1)
    bs <- boot.mirt(mod_items.group, R=2)
    expect_is(bs, 'boot')


    #model using 2PL items instead of only Rasch, and with missing data
    data[1,1] <- covdata[1,2] <- NA
    mod1b <- mixedmirt(data, covdata, model, fixed = ~ 0 + items + group,
                                        itemtype = '2PL', verbose = FALSE, draws = 1)
    expect_is(mod1b, 'MixedClass')
    expect_equal(extract.mirt(mod1b, 'df'), 1001)
    cfs <- as.numeric(do.call(c, coef(mod1b)))
    expect_equal(cfs, c(1.0652,NaN,NaN,2.3363,NaN,NaN,0.6912,NaN,NaN,-1.7727,NaN,NaN,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.3748,NaN,NaN,-2.064,-2.2495,-1.8785,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.1637,NaN,NaN,-1.6932,NaN,NaN,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,-0.0519,NaN,NaN,-1.0393,NaN,NaN,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.3548,NaN,NaN,-0.2618,NaN,NaN,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.0425,NaN,NaN,-1.3679,NaN,NaN,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.3055,-0.113,0.724,-1.6781,-1.7586,-1.5976,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.1485,NaN,NaN,-2.1073,-2.1676,-2.047,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,1.0666,NaN,NaN,-2.1148,NaN,NaN,0,NA,NA,1,NA,NA,1.0652,NaN,NaN,2.3363,NaN,NaN,0.831,NaN,NaN,1.7523,1.4181,2.0865,0,NA,NA,1,NA,NA,0,NA,NA,1,NA,NA),
                 tolerance = 1e-2)

    covdata$group <- factor(rep(paste0('G',1:50), each = N/50))
    rmod1 <- try(mixedmirt(data, covdata, 1, fixed = ~ 0 + items, random = ~ 1|group,
                                        draws = 1, verbose = FALSE), TRUE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(extract.mirt(rmod1, 'df'), 1011)
    cfs <- as.numeric(do.call(c, coef(rmod1)))
    expect_equal(cfs[124:129], c(0.06756121, 0.05018307, 0.08493935, 1.13460966, 0.66202405, 1.60719526),
                 tolerance = 1e-2)

    #polytomous
    covdat <- data.frame(group = rep(c('m', 'f'), nrow(Science)/2))
    model <- mirt.model('F1 = 1-4', quiet = TRUE)
    mod <- mixedmirt(Science, covdat, model=model, SE=FALSE,
                     fixed = ~ 0 + group, verbose = FALSE, draws = 1)
    expect_is(mod, 'MixedClass')
    cfs <- as.numeric(na.omit(do.call(c, coef(mod))))
    expect_equal(cfs, c(-0.0515,1,0,1,2,3,0,3.0327,5.6104,4.2493,-0.0515,1,0,1,2,3,0,1.8659,2.7745,0.958,-0.0515,1,0,1,2,3,0,2.6075,4.0176,2.9117,-0.0515,1,0,1,2,3,0,2.4149,3.3088,1.9863,0,0.9265),
                 tolerance = 1e-2)

    mod2 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'gpcm', verbose = FALSE)
    expect_is(mod2, 'MixedClass')
    expect_equal(extract.mirt(mod, 'df') - extract.mirt(mod2, 'df'), 3)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod2))))
    expect_equal(cfs, c(-0.1658,-0.3961,0.0645,0.8385,0.4843,1.1928,0,1,2,3,0,2.8555,1.7005,4.0105,5.4041,4.1164,6.6918,4.1651,2.9216,5.4087,-0.1658,-0.3961,0.0645,0.822,0.5568,1.0872,0,1,2,3,0,1.7642,1.2394,2.2889,2.7013,2.086,3.3166,1.0618,0.396,1.7277,-0.1658,-0.3961,0.0645,2.398,0.9417,3.8543,0,1,2,3,0,4.9963,2.4428,7.5497,7.3389,3.7964,10.8813,5.4399,2.6573,8.2226,-0.1658,-0.3961,0.0645,0.7278,0.4345,1.0211,0,1,2,3,0,2.1803,1.5836,2.777,3.057,2.3766,3.7374,1.9533,1.2594,2.6473,0,1),
                 tolerance = 1e-2)

    mod3 <- mixedmirt(Science, covdat, model=model, draws = 1,
                                       fixed = ~ 0 + group, itemtype = 'graded', verbose = FALSE)
    expect_is(mod3, 'MixedClass')
    expect_equal(extract.mirt(mod3, 'df'), 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(mod3))))
    expect_equal(cfs, c(-0.2005,-0.5234,0.1223,1.0151,0.6339,1.3964,4.937,3.968,5.9061,2.7168,2.2619,3.1717,-1.3634,-1.713,-1.0138,-0.2005,-0.5234,0.1223,1.1613,0.8347,1.4878,2.969,2.4977,3.4403,0.9771,0.666,1.2882,-2.1311,-2.5415,-1.7206,-0.2005,-0.5234,0.1223,2.5549,1.2906,3.8192,5.695,3.7818,7.6083,2.4632,1.5444,3.382,-2.0238,-2.8748,-1.1729,-0.2005,-0.5234,0.1223,1.0885,0.719,1.4579,3.4395,2.8808,3.9982,1.0818,0.7716,1.392,-1.5923,-1.9637,-1.2208,0,1),
                 tolerance = 1e-2)

    covdat$group <- factor(rep(paste0('G',1:20), length.out = nrow(Science)))
    rmod1 <- try(mixedmirt(Science, covdat, model=model, draws=10, random = ~ 1|group,
                       itemtype = 'graded', verbose = FALSE, SE=FALSE), TRUE)
    expect_is(rmod1, 'MixedClass')
    expect_equal(extract.mirt(rmod1, 'df'), 238)
    cfs <- as.numeric(na.omit(do.call(c, coef(rmod1))))
    expect_equal(cfs, c(1.062578,4.895614,2.661662,-1.479457,1.195101,2.907836,0.897937,-2.253751,2.178545,5.083139,2.151222,-1.918151,1.08017,3.345869,0.9912499,-1.68683,0,1,0.006179096),
                 tolerance = 1e-4)
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
    cfs <- coef(mod, printSE=TRUE)
    expect_equal(as.numeric(cfs$lr.betas)[-c(1:2)], c(1.4192996, 0.1366249, 3.1089497, 0.1664278),
                 tolerance=1e-4)
    expect_equal(extract.mirt(mod, 'logLik'), -4314.043, tolerance = 1e-2)
    expect_equal(extract.mirt(mod, 'df'), 1001)

    mod2 <- mixedmirt(dat, covdata = covdata, 1, itemtype = 'Rasch', fixed = ~ 0 + items,
                     lr.fixed = ~ group, verbose=FALSE, draws=1, SE=FALSE)
    so <- summary(mod2, verbose=FALSE)
    expect_equal(as.numeric(c(so$random$Theta, so$lr.out[,1])),
                 c(0.1978728, 0.0000000, 0.7181347, 1.7006670), tolerance=1e-3)
    expect_equal(extract.mirt(mod2, 'logLik'), -4649.466, tolerance = 1e-4)
    set.seed(1)
    bs <- boot.mirt(mod2, R = 3)
    expect_is(bs, 'boot')

    #uncorrelated random slope
    covdata$theta <- Theta
    covdata$cut <- factor(cut(Theta, breaks=2))
    mod <- try(mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items, SE=FALSE,
                     random = ~ -1 + cut|group, verbose=FALSE, draws=1), TRUE)
    so <- summary(mod, verbose=FALSE)
    expect_equal(as.numeric(diag(so$random$group)), c(0.6155525, 0.4089267), tolerance = 1e-4)

    mod <- try(mixedmirt(dat, covdata = covdata, 1, fixed = ~ 0 + items, SE=FALSE,
                     random = ~ -1 + theta|group, verbose=FALSE, draws=1), TRUE)
    so <- summary(mod, verbose=FALSE)
    expect_equal(as.numeric(diag(so$random$group)), c(0.2255360, 0.4568764), tolerance = 1e-4)

})
