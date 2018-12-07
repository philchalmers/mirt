context('DIF')

test_that('DIF', {

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
    model1a <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE)
    model2 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE,
                            invariance = c('slopes', 'intercepts', 'free_means', 'free_var'))

    DIFF <- suppressMessages(DIF(model1a, which.par='d', items2test = c(1,3)))
    expect_is(DIFF, 'data.frame')
    expect_equal(DIFF[2,'AIC'], -3.414, tolerance = 1e-3)
    DIFF2 <- suppressMessages(DIF(model2, which.par=c('a1', 'd'), items2test = c(1,3), scheme='drop'))
    expect_is(DIFF2, 'data.frame')
    expect_equal(DIFF2[1L, 'X2'], 18.88, tolerance = 1e-3)

    WALD <- suppressMessages(DIF(model1a, which.par='d', items2test = 1:3, Wald=TRUE))
    expect_is(WALD, 'data.frame')
    expect_equal(WALD$W[1], 1.532718, tolerance = 1e-3)
    expect_equal(WALD$p[1], 0.2157049, tolerance = 1e-3)
    WALD2 <- suppressMessages(DIF(model1a, which.par=c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr'))
    expect_equal(as.numeric(WALD2$adj_pvals), c(0.02473307,0.01580279,0.1261503,0.4299811,0.4672727,0.4565218,0.05971613,0.5671398,0.4299811,0.4565218), tolerance = 1e-3)


    model1b <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, invariance = c(colnames(dat)[1:5],
                                                                                     'free_means', 'free_var'))

    out <- DIF(model1b, which.par = c('a1', 'd'), items2test = 6:10, Wald = TRUE)
    expect_is(out, 'data.frame')
    expect_equal(out$p, c(0.1090403,0.07219868,0.175297,0.1354964,0.6970995), tolerance = 1e-4)

    out <- DIF(model1b, which.par = c('a1', 'd'), items2test = 6:10)
    expect_is(out, 'data.frame')
    expect_equal(out$p, c(0.103812,0.04867894,0.1565903,0.1118176,0.6921238), tolerance = 1e-4)

    model <- "F = 1-10
              PRIOR = (6, a1, lnorm, .2, .3)"
    model2 <- multipleGroup(dat, model, group, SE = TRUE, verbose=FALSE, invariance = c(colnames(dat)[1:5],
                                                                                     'free_means', 'free_var'))
    out <- DIF(model2, which.par = c('a1', 'd'), items2test = 6:10, seq_stat = 'DIC')
    expect_is(out, 'data.frame')
    expect_equal(out$Bayes_Factor, c(0.08699142,0.03477748,0.1278902,0.08952973,0.6469454), tolerance = 1e-4)


})


