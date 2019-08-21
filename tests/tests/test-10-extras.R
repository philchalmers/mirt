context('extras')

test_that('extras', {
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    data <- expand.table(LSAT7)
    data <- rbind(data, data)
    mod1 <- mirt(data, 1, verbose=FALSE, SE=TRUE)

    fun <- function(Thetas, MIN, MAX, ...) as.numeric(dunif(Thetas, min=MIN, max=MAX))
    fs1 <- fscores(mod1, verbose = FALSE, custom_den = fun, MIN = -3, MAX = 3, full.scores=FALSE)
    fs2 <- suppressWarnings(fscores(mod1, custom_den = fun, MIN = -3, MAX = 3, verbose = FALSE, method = 'MAP',
                                    full.scores=FALSE))
    expect_equal(as.numeric(fs1[1,c('F1', 'SE_F1')]), c(-2.4766, .4988), tolerance = 1e-3)
    expect_equal(as.numeric(fs2[5,c('F1', 'SE_F1')]), c(-1.9121, .9381), tolerance = 1e-3)

    set.seed(12345)
    a1 <- a2 <- matrix(abs(rnorm(10,1,.3)), ncol=1)
    d1 <- d2 <- matrix(rnorm(10,0,.7),ncol=1)
    a2[1:2, ] <- a1[1:2, ]/3
    d1[c(1,3), ] <- d2[c(1,3), ]/4
    itemtype <- rep('dich', nrow(a1))
    N <- 1000
    dataset1 <- simdata(a1, d1, N, itemtype)
    dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
    dat <- rbind(dataset1, dataset2)
    group <- c(rep('D1', N), rep('D2', N))
    model1a <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, SE.type = 'Louis')
    model1b <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, SE.type = 'Richardson',
                             pars = mod2values(model1a), technical = list(warn=FALSE))
    expect_equal(as.numeric(model1a@vcov - model1b@vcov), numeric(ncol(model1a@vcov)^2),
                 tolerance = 1e-4)
    model2 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE,
                            invariance = c('slopes', 'intercepts', 'free_means', 'free_var'))
    modideal <- mirt(dataset1, model = mirt.model('F1 = 1-6
                                                  F2 = 5-10'), 'ideal', verbose = FALSE)
    cfs <- as.numeric(coef(modideal, verbose=FALSE)[[5]])
    expect_equal(extract.mirt(modideal, 'logLik'), -6408.54, tolerance = 1e-3)
    expect_equal(cfs, c(0.3044798,0.4071647,-1.432097), tolerance = 1e-3)

    acov <- fscores(mod1, return.acov=TRUE, full.scores=FALSE)
    expect_equal(acov[[1]][1], 0.4799239, tolerance=1e-3)
    acov <- fscores(mod1, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], 0.4987313, tolerance=1e-3)
    acov <- fscores(model1a, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], .325383, tolerance=1e-3)
    boot <- boot.mirt(model1a, R=5)
    cfs <- boot$t0
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(as.numeric(cfs), c(1.200941,0.03964278,1.202968,1.253258,1.162712,0.1054011,0.9359284,0.410968,1.20204,-0.5220437,0.5644495,0.696662,1.219704,-0.7151586,0.9018765,-0.1465849,0.8550952,0.7757281,0.7616766,0.1305874,0.7119052,-0.09370059,0.6546074,1.274378,1.36438,0.38146,1.207508,0.5272688,1.370619,-0.426826,0.6340348,0.57111,1.651928,-0.5331716,1.068324,-0.1782871,1.121677,0.8790566,0.7972691,0.2748593),
                 tolerance=1e-3)
    PLCI <- PLCI.mirt(mirt(Science, 1, verbose=FALSE), parnum=c(1,2))
    expect_equal(c(PLCI$lower_2.5, PLCI$upper_97.5), c(0.7008617, 4.0112247, 1.4530232, 5.9667792),
                  tolerance=1e-3)

    extr.2 <- extract.item(mod1, 2)
    Theta <- matrix(seq(-6,6, length.out=200))
    expected <- expected.item(extr.2, Theta)
    expect_equal(expected[1:3], c(0.003411029, 0.003639925, 0.003884122), tolerance=1e-4)
    expected <- expected.test(mod1, Theta=Theta)
    expect_equal(expected[1:3], c(0.1083150, 0.1133415, 0.1185956), tolerance=1e-4)
    data[1,1] <- NA
    data <- imputeMissing(mod1, Theta=fscores(mod1, full.scores=TRUE), warn=FALSE)
    expect_is(data, 'matrix')
    expect_true(!all(is.na(data)))
    pb <- probtrace(extr.2, Theta)
    expect_equal(pb[1:3, 1], c(0.9965890, 0.9963601, 0.9961159), tolerance = 1e-5)
    ti <- testinfo(mod1, Theta=Theta)
    expect_equal(ti[1:3], c(0.06607895, 0.06918330, 0.07242870), tolerance = 1e-5)
    set.seed(1234)
    fs <- fscores(mod1, MI=20, verbose=FALSE, full.scores=FALSE)
    expect_is(fs, 'matrix')
    expect_equal(fs[1:3,'F1'], c(-1.853914, -1.509722, -1.514913), tolerance=1e-3)

    set.seed(1)
    dat <- data.frame(Science, Science + sample(c(-1,0,1), prod(dim(Science)), TRUE))
    mats <- mats2 <- list()
    mats[1:4] <- mats2[1:4] <- list(matrix(c(0:3, 0:3), 4))
    mats[5:8] <- list(matrix(c(0:5, 1,1,0,0,0,0), 6))
    mats2[5:8] <- list(matrix(c(0:5, 0:5), 6))
    mod1 <- mirt(dat, 2, 'gpcm', TOL = 5e-2, verbose=FALSE)
    mod2 <- mirt(dat, 2, 'gpcm', gpcm_mats = mats2, TOL = 5e-2, verbose=FALSE)
    s1 <- coef(mod1, simplify=TRUE)$items
    s2 <- coef(mod2, simplify=TRUE)$items
    pick <- c('a1', 'a2', 'd1', 'd2', 'd3')
    expect_true(sum(abs(s1[,pick] - s2[,pick])) < 1e-8)
    mod3 <- mirt(dat, 2, 'gpcm', gpcm_mats = mats, TOL = 1e-2, verbose=FALSE)
    expect_equal(extract.mirt(mod3, 'logLik'), -3721.461, tolerance = 1e-4)
    cfs <- as.vector(coef(mod3, simplify=TRUE)$items)
    expect_equal(cfs, c(0.7835446,2.088519,4.489773,0.6727285,0.2357904,0.3148114,1.736133,0.3861383,0.07062774,6.563708,1.01623,0.05366885,0.2764068,-5.030398,-0.1089204,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,2.737057,9.954572,8.696347,2.051471,1.746521,1.493303,4.871186,2.184036,5.174404,13.00528,12.46809,2.830379,4.335373,6.261805,7.718796,3.259992,3.88189,5.255669,9.140148,1.683119,4.413986,6.18042,8.669859,3.579751,NA,NA,NA,NA,4,4,4,4,NA,NA,NA,NA,5,5,5,5,NA,NA,NA,NA,0,0,0,0,NA,NA,NA,NA,0,0,0,0,NA,NA,NA,NA,4.178023,5.776866,7.750279,3.300487,NA,NA,NA,NA,3.043197,4.313529,5.190609,1.550788), tolerance=1e-4)

})

