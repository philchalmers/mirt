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
    expect_equal(cfs, c(0.3044798, 0.4071647, -1.4320971), tolerance = 1e-3)

    acov <- fscores(mod1, return.acov=TRUE, full.scores=FALSE)
    expect_equal(acov[[1]][1], 0.4799239, tolerance=1e-3)
    acov <- fscores(mod1, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], 0.4987313, tolerance=1e-3)
    acov <- fscores(model1a, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], .325383, tolerance=1e-3)
    boot <- boot.mirt(model1a, R=5)
    cfs <- boot$t0
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(as.numeric(cfs), c(1.200932,0.03966039,1.202983,1.253283,1.162713,0.1054175,0.935932,0.4109794,1.202011,-0.5220202,0.5644482,0.6966761,1.219665,-0.7151297,0.9018733,-0.1465739,0.855102,0.7757409,0.7616776,0.1306053,0.3678973,-0.1117068,0.4865936,1.258095,1.270467,0.3109331,1.126136,0.3280035,1.441583,-0.3823625,0.5785964,0.6331239,1.531089,-0.5576562,1.289808,-0.1500941,1.141839,0.94723,0.9748633,0.2646415),
                 tolerance=1e-3)
    PLCI <- PLCI.mirt(mirt(Science, 1, verbose=FALSE), parnum=c(1,2))
    expect_equal(c(PLCI$lower_2.5, PLCI$upper_97.5), c(0.7008617, 4.0112247, 1.4530232, 5.9667792),
                  tolerance=1e-3)
    DIFF <- suppressMessages(DIF(model1a, which.par='d', items2test = c(1,3)))
    expect_is(DIFF, 'list')
    expect_equal(DIFF[[1L]][2,'logLik'], -12519.45, tolerance = 1e-3)
    DIFF2 <- suppressMessages(DIF(model2, which.par=c('a1', 'd'), items2test = c(1,3), scheme='drop'))
    expect_is(DIFF2, 'list')
    expect_equal(DIFF2[[1L]][,'logLik'], c(-12562.00, -12539.75), tolerance = 1e-3)

    WALD <- suppressMessages(DIF(model1a, which.par='d', items2test = 1:3, Wald=TRUE))
    expect_is(WALD, 'list')
    expect_equal(WALD[[1]]$W[1], 2.101652, tolerance = 1e-3)
    expect_equal(WALD[[1]]$p[1], .1471401, tolerance = 1e-3)
    WALD2 <- suppressMessages(DIF(model1a, which.par=c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr'))
    expect_equal(as.numeric(WALD2$adj_pvals), c(2.027518e-06,0.0001081836,0.2437397,0.370737,0.2437397,0.8028216,0.2010698,0.2010698,0.2437397,0.2437397), tolerance = 1e-3)
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
    expect_true(sum(abs(s1[,pick] - s2[,pick])) < 1e-10)
    mod3 <- mirt(dat, 2, 'gpcm', gpcm_mats = mats, TOL = 1e-2, verbose=FALSE)
    expect_equal(extract.mirt(mod3, 'logLik'), -3708.216, tolerance = 1e-4)
    cfs <- as.vector(coef(mod3, simplify=TRUE)$items)
    expect_equal(cfs, c(-1.23177,-2.76946,-1.54897,-1.34947,-0.41973,-0.45958,-0.71595,-0.58334,-0.29564,4.58694,0.36752,-0.31281,0.4241,-2.90723,-0.16516,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,3.69493,7.87579,3.59178,3.00087,2.60806,2.7068,2.83078,2.07429,6.50601,10.34317,5.36289,4.0315,5.42778,5.57857,4.69492,3.43515,4.9848,4.17678,3.89073,2.40695,5.74,5.87555,5.00542,3.86719,NA,NA,NA,NA,4,4,4,4,NA,NA,NA,NA,5,5,5,5,NA,NA,NA,NA,0,0,0,0,NA,NA,NA,NA,0,0,0,0,NA,NA,NA,NA,5.43095,5.52149,4.54781,3.38757,NA,NA,NA,NA,3.98583,3.27075,2.90466,1.19029), tolerance=1e-4)

})

