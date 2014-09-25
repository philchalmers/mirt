context('extras')

test_that('extras', {
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    data <- expand.table(LSAT7)
    data <- rbind(data, data)
    mod1 <- mirt(data, 1, verbose=FALSE, SE=TRUE)
    
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
    expect_equal(as.numeric(model1a@information - model1b@information), numeric(ncol(model1a@information)^2),
                 tolerance = 1e-3)
    model2 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, 
                            invariance = c('slopes', 'intercepts', 'free_means', 'free_var'))
    modideal <- mirt(dataset1, model = mirt.model('F1 = 1-6
                                                  F2 = 5-10'), 'ideal', verbose = FALSE)
    cfs <- as.numeric(coef(modideal, digits=5, verbose=FALSE)[[5]])
    expect_equal(modideal@logLik, -6435.657, tolerance = 1e-3)
    expect_equal(cfs, c(0.73246, 1.26111, -1.41768), tolerance = 1e-3)
    
    acov <- fscores(mod1, return.acov=TRUE)
    expect_equal(acov[[1]][1], 0.4799239, tolerance=1e-3)
    acov <- fscores(mod1, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], 0.4987313, tolerance=1e-3)
    acov <- fscores(model1a, return.acov=TRUE, full.scores=TRUE)
    expect_equal(acov[[500]][1], .3452783, tolerance=1e-3)
    boot <- boot.mirt(model1a, R=5)
    cfs <- boot$t0
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.03027, -0.03172, 1.18018, 1.36523, 1.1824, 0.10923, 0.94459, 0.36701, 1.09781, -0.53252, 0.38686, 0.48484, 1.4409, -0.6983, 0.80128, -0.22106, 0.80456, 0.83075, 0.7102, 0.24493, 0.67924, 0.10731, 0.47428, 1.27867, 1.42737, 0.41311, 1.29134, 0.51946, 1.36572, -0.39269, 0.64571, 0.60175, 1.44829, -0.4938, 1.2331, -0.20888, 1.08338, 0.95549, 0.98608, 0.1816),
                 tolerance=1e-3)
    PLCI <- PLCI.mirt(mod1, parnum=c(1,2))
    expect_equal(c(PLCI$lower_2.5, PLCI$upper_97.5), c(0.7580446, 1.6843469, 1.2564076, 2.0529152),
                 tolerance=1e-3)
    DIFF <- DIF(model1a, which.par='d', items2test = 1:3)
    expect_is(DIFF, 'list')
    expect_equal(DIFF[[1L]][2,'logLik'], -12508.15, tolerance = 1e-3)    
    DIFF2 <- DIF(model2, which.par=c('a1', 'd'), items2test = 1:3, scheme='drop')
    expect_is(DIFF2, 'list')
    expect_equal(DIFF2[[1L]][,'logLik'], c(-12534.13, -12528.83), tolerance = 1e-3)    
    
    WALD <- DIF(model1a, which.par='d', items2test = 1:3, Wald=TRUE)
    expect_is(WALD, 'list')
    expect_equal(WALD[[1]]$W[1], 1.779187, tolerance = 1e-3)
    expect_equal(WALD[[1]]$p[1], .1822492, tolerance = 1e-3)
    WALD2 <- DIF(model1a, which.par=c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr') 
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
    expect_equal(ti[1:3, 1], c(0.06607895, 0.06918330, 0.07242870), tolerance = 1e-5)
    set.seed(1234)
    fs <- fscores(mod1, MI=20, verbose=FALSE)
    expect_is(fs, 'matrix')
    expect_equal(fs[1:3,'F1'], c(-1.853914, -1.509722, -1.514913), tolerance=1e-3)
})

