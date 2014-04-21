context('extras')

test_that('extras', {
    require(boot, quietly=TRUE, warn.conflicts=FALSE)
    data <- expand.table(LSAT7)
    data <- rbind(data, data)
    mod1 <- mirt(data, 1, verbose=FALSE)
    
    set.seed(12345)
    a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
    d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
    a2[1:2, ] <- a1[1:2, ]/3
    d1[c(1,3), ] <- d2[c(1,3), ]/4
    head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
    itemtype <- rep('dich', nrow(a1))
    N <- 1000
    dataset1 <- simdata(a1, d1, N, itemtype)
    dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
    dat <- rbind(dataset1, dataset2)
    group <- c(rep('D1', N), rep('D2', N))
    model1 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE)
    model2 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, 
                            invariance = c('slopes', 'intercepts', 'free_means', 'free_var'))
    
    boot <- boot.mirt(model1, R=5)
    cfs <- boot$t0
    cfs <- cfs[cfs != 0 & cfs != 1]
    expect_equal(cfs, c(1.1319, 0.1772, 1.3074, -0.6995, 0.9705, 0.1217, 1.1235, 0.8317, 1.2338, 0.3238, 0.4847, 0.4804, 1.1465, 1.0791, 0.8362, -0.3832, 0.8757, -1.0441, 0.8032, -1.0894, 0.9031, 1.1644, 1.5724, -0.1368, 1.3984, 0.6504, 1.0773, 0.4108, 0.8692, -0.0815, 0.5239, 0.6302, 0.5529, -0.5649, 1.3728, -0.2063, 1.2107, 1.0205, 1.258, 0.4231, 0.4662, 0.5413, 1.5702, 1.3373, 1.1292, -0.3902, 0.9303, -0.9375, 0.9328, -1.0265, 1.1988, 1.3824, 1.9066, -0.1398, 1.2711, 0.5115, 1.4832, 0.5471, 1.1101, -0.0948),
                 tolerance=1e-3)
    PLCI <- PLCI.mirt(mod1, parnum=c(1,2))
    expect_equal(c(PLCI$lower_2.5, PLCI$upper_97.5), c(0.7580446, 1.6843469, 1.2564076, 2.0529152),
                 tolerance=1e-3)
    DIFF <- DIF(model1, which.par='d', items2test = 1:3)
    expect_is(DIFF, 'list')
    expect_equal(DIFF[[1L]][2,'logLik'], -18035.15, tolerance = 1e-3)
    
    WALD <- DIF(model1, which.par='d', items2test = 1:3, Wald=TRUE)
    extr.2 <- extract.item(mod1, 2)
    Theta <- matrix(seq(-6,6, length.out=200))
    expected <- expected.item(extr.2, Theta)
    expect_equal(expected[1:3], c(0.003411474, 0.003640396, 0.003884619), tolerance=1e-4)
    expected <- expected.test(mod1, Theta=Theta)
    expect_equal(expected[1:3], c(0.1083644, 0.1133925, 0.1186483), tolerance=1e-4)
    data[1,1] <- NA
    data <- imputeMissing(mod1, Theta=fscores(mod1, full.scores=TRUE))
    expect_is(data, 'matrix')
    expect_true(!all(is.na(data)))
    pb <- probtrace(extr.2, Theta)
    expect_equal(pb[1:6, 1], c(0.9965885, 0.9963596, 0.9961154, 0.9958548, 0.9955769, 0.9952804),
                 tolerance = 1e-5)
    ti <- testinfo(mod1, Theta=Theta)
    expect_equal(ti[1:6, 1], c(0.06609361, 0.06919815, 0.07244372, 0.07583628, 0.07938197, 0.08308716),
                 tolerance = 1e-5)
})

