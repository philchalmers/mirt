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
    model1b <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, SE.type = 'BL',
                             pars = mod2values(model1a), technical = list(warn=FALSE))
    expect_equal(as.numeric(model1a@vcov - model1b@vcov), numeric(ncol(model1a@vcov)^2),
                 tolerance = 1e-4)
    model2 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE,
                            invariance = c('slopes', 'intercepts', 'free_means', 'free_var'))
    modideal <- mirt(dataset1, model = mirt.model('F1 = 1-6
                                                  F2 = 5-10'), 'ideal', verbose = FALSE)
    cfs <- as.numeric(coef(modideal, digits=5, verbose=FALSE)[[5]])
    expect_equal(extract.mirt(modideal, 'logLik'), -6435.647, tolerance = 1e-3)
    expect_equal(cfs, c(0.73123, 1.26069, -1.41785), tolerance = 1e-3)

    acov <- fscores(mod1, return.acov=TRUE, full.scores=FALSE)
    expect_equal(acov[[1]][1], 0.4799239, tolerance=1e-3)
    acov <- fscores(mod1, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], 0.4987313, tolerance=1e-3)
    acov <- fscores(model1a, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], .3452783, tolerance=1e-3)
    boot <- boot.mirt(model1a, R=5)
    cfs <- boot$t0
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(as.numeric(cfs), c(1.03027, -0.03172, 1.18018, 1.36523, 1.1824, 0.10923, 0.94459, 0.36701, 1.09781, -0.53252, 0.38686, 0.48484, 1.4409, -0.6983, 0.80128, -0.22106, 0.80456, 0.83075, 0.7102, 0.24493, 0.67924, 0.10731, 0.47428, 1.27867, 1.42737, 0.41311, 1.29134, 0.51946, 1.36572, -0.39269, 0.64571, 0.60175, 1.44829, -0.4938, 1.2331, -0.20888, 1.08338, 0.95549, 0.98608, 0.1816),
                 tolerance=1e-3)
    PLCI <- PLCI.mirt(mod1, parnum=c(1,2))
    expect_equal(c(PLCI$lower_2.5, PLCI$upper_97.5), c(0.7580446, 1.6843469, 1.2564076, 2.0529152),
                 tolerance=1e-3)
    DIFF <- suppressMessages(DIF(model1a, which.par='d', items2test = 1:3))
    expect_is(DIFF, 'list')
    expect_equal(DIFF[[1L]][2,'logLik'], -12508.15, tolerance = 1e-3)
    DIFF2 <- suppressMessages(DIF(model2, which.par=c('a1', 'd'), items2test = 1:3, scheme='drop'))
    expect_is(DIFF2, 'list')
    expect_equal(DIFF2[[1L]][,'logLik'], c(-12534.13, -12528.83), tolerance = 1e-3)

    WALD <- suppressMessages(DIF(model1a, which.par='d', items2test = 1:3, Wald=TRUE))
    expect_is(WALD, 'list')
    expect_equal(WALD[[1]]$W[1], 1.779187, tolerance = 1e-3)
    expect_equal(WALD[[1]]$p[1], .1822492, tolerance = 1e-3)
    WALD2 <- suppressMessages(DIF(model1a, which.par=c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr'))
    expect_equal(as.numeric(WALD2$adj_pvals), c(0.0777, 0.0010, 0.0777, 0.1611, 0.1652, 0.1611,
                                                0.2662, 0.0777, 0.2344, 0.1652), tolerance = 1e-3)
    extr.2 <- extract.item(mod1, 2)
    Theta <- matrix(seq(-6,6, length.out=200))
    expected <- expected.item(extr.2, Theta)
    expect_equal(expected[1:3], c(0.003411029, 0.003639925, 0.003884122), tolerance=1e-4)
    expected <- expected.test(mod1, Theta=Theta)
    expect_equal(expected[1:3], c(0.1083150, 0.1133415, 0.1185956), tolerance=1e-4)
    data[1,1] <- NA
    data <- imputeMissing(mod1, Theta=fscores(mod1, full.scores=TRUE))
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
    dat <- cbind(Science, Science + sample(c(-1,0,1), prod(dim(Science)), TRUE))
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
    expect_equal(extract.mirt(mod3, 'logLik'), -3708.0, tolerance = 1e-4)
    cfs <- as.vector(coef(mod3, simplify=TRUE, digits = 5)$items)
    expect_equal(cfs, c(-1.23177,-2.76946,-1.54897,-1.34947,-0.41973,-0.45958,-0.71595,-0.58334,-0.29564,4.58694,0.36752,-0.31281,0.4241,-2.90723,-0.16516,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,0,0,0,0,3,3,3,3,0,0,0,0,0,0,0,0,0,0,0,0,3.69493,7.87579,3.59178,3.00087,2.60806,2.7068,2.83078,2.07429,6.50601,10.34317,5.36289,4.0315,5.42778,5.57857,4.69492,3.43515,4.9848,4.17678,3.89073,2.40695,5.74,5.87555,5.00542,3.86719,NA,NA,NA,NA,4,4,4,4,NA,NA,NA,NA,5,5,5,5,NA,NA,NA,NA,0,0,0,0,NA,NA,NA,NA,0,0,0,0,NA,NA,NA,NA,5.43095,5.52149,4.54781,3.38757,NA,NA,NA,NA,3.98583,3.27075,2.90466,1.19029), tolerance=1e-4)

})

