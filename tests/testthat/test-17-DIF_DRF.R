expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('DIF', {

    if(FALSE){
        rm(list=ls())
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
        save(dat, group, file = 'tests/tests/testdata/dif1.rds')
    }
    load('testdata/dif1.rds')
    model1a <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE)
    model2 <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE,
                            invariance = c('slopes', 'intercepts', 'free_means', 'free_var'))

    DIFF <- suppressMessages(DIF(model1a, which.par='d', items2test = c(1,3)))
    expect_class(DIFF, 'data.frame')
    expect_equal(DIFF[2,'AIC'], -3.414, tolerance = 1e-3)
    DIFF2 <- suppressMessages(DIF(model2, which.par=c('a1', 'd'), items2test = c(1,3), scheme='drop'))
    expect_class(DIFF2, 'data.frame')
    expect_equal(DIFF2[1L, 'X2'], 18.88, tolerance = 1e-3)

    WALD <- suppressMessages(DIF(model1a, which.par='d', items2test = 1:3, Wald=TRUE))
    expect_class(WALD, 'data.frame')
    expect_equal(WALD$W[1], 1.532718, tolerance = 1e-3)
    expect_equal(WALD$p[1], 0.2157049, tolerance = 1e-3)
    WALD2 <- suppressMessages(DIF(model1a, which.par=c('a1', 'd'), Wald=TRUE, p.adjust = 'fdr'))
    expect_equal(as.numeric(WALD2$adj_p), c(0.02473307,0.01580279,0.1261503,0.4299811,0.4672727,0.4565218,0.05971613,0.5671398,0.4299811,0.4565218), tolerance = 1e-3)


    model1b <- multipleGroup(dat, 1, group, SE = TRUE, verbose=FALSE, invariance = c(colnames(dat)[1:5],
                                                                                     'free_means', 'free_var'))

    out <- DIF(model1b, which.par = c('a1', 'd'), items2test = 6:10, Wald = TRUE)
    expect_class(out, 'data.frame')
    expect_equal(out$p, c(0.1090403,0.07219868,0.175297,0.1354964,0.6970995), tolerance = 1e-4)

    out <- DIF(model1b, which.par = c('a1', 'd'), items2test = 6:10)
    expect_class(out, 'data.frame')
    expect_equal(out$p, c(0.103812,0.04867894,0.1565903,0.1118176,0.6921238), tolerance = 1e-4)

    model <- "F = 1-10
              PRIOR = (6, a1, lnorm, .2, .3)"
    model2 <- multipleGroup(dat, model, group, SE = TRUE, verbose=FALSE, invariance = c(colnames(dat)[1:5],
                                                                                     'free_means', 'free_var'))
    out <- DIF(model2, which.par = c('a1', 'd'), items2test = 6:10, seq_stat = 'BIC')
    expect_class(out, 'data.frame')
    expect_equal(out$BIC, c(11.449,9.038,11.424,10.739,14.448), tolerance = 1e-3)

    drf <- DRF(model2)
    expect_equal(as.numeric(subset(drf, select=n_focal_items:dDRF)),
                            c(10.0000000,  -0.0403703,  0.1863619, 0.2107826),
                 tolerance = 1e-3)
    dif <- DRF(model2, DIF = TRUE)
    expect_equal(as.numeric(dif$sDIF[-c(1:5)]),
                 c(-0.04435560,  0.02738697, -0.02307566, -0.01429807,  0.01397205),
                 tolerance = 1e-4)
    expect_equal(as.numeric(dif$uDIF[-c(1:5)]),
                 c(0.04448119, 0.06111498, 0.04441609, 0.04817261, 0.01960173),
                 tolerance = 1e-4)

    set.seed(1234)
    drf <- DRF(model2, draws = 100)
    expect_equal(drf$X2, c(0.4536206, 16.6727036, NA), tolerance=1e-2)

})

test_that('MG_DIF', {

    set.seed(1234)
    nitem = 10

    n.g1 = 2000
    n.g2 = 1000
    n.g3 = 1000
    n.g4 = 1000

    groups = c(rep('g1', n.g1), rep('g2', n.g2), rep('g3', n.g3), rep('g4', n.g4))

    a <- rep(1, nitem)
    d <- matrix(rnorm(nitem,0,.7),ncol=1)
    item.type <- rep('dich', nitem)

    # Generate data
    dat1 <- simdata(a, d, n.g1, item.type)
    dat2 <- simdata(a, d, n.g2, item.type)
    dat3 <- simdata(a, d, n.g3, item.type)
    dat4 <- simdata(a, d, n.g4, item.type)

    dat <- rbind(dat1,dat2,dat3,dat4)

    anchor.set <- colnames(dat)[1:5]

    res.model <- multipleGroup(dat, 1,
                               group = groups,
                               invariance = c("free_var", "free_mean", anchor.set),
                               SE = TRUE, verbose=FALSE)


    # joint test
    out <- DIF(res.model, c('a1', 'd'), items2test=6)
    expect_equal(out$p, 0.2584821, tolerance = 1e-2)
    out2 <- DIF(res.model, c('a1', 'd'), Wald = TRUE, items2test=6:10)
    expect_equal(out2$p[1:2], c(0.3162816, 0.5118902), tolerance = 1e-2)

    # group subsets
    out3 <- DIF(res.model, c('d'), Wald = F, items2test = 6:7,
                groups2test = c('g1', 'g2'), verbose=FALSE)
    expect_equal(out3$p[1:2], c(0.5264586, 0.2571904), tolerance = 1e-2)
    out4 <- DIF(res.model, c('d'), Wald = T, items2test = 6:7, groups2test = c('g1', 'g2'))
    expect_equal(out4$p[1:2], c(0.5261373, 0.2604559), tolerance = 1e-2)

    out5 <- DIF(res.model, c('a1', 'd'), Wald = T, items2test = 6, groups2test = c('g1', 'g3'))
    expect_equal(out5$p, c(0.1533536), tolerance = 1e-2)

    out6 <- DIF(res.model, c('d'), Wald = F, items2test = 6, groups2test = c('g1', 'g2', 'g3'))
    expect_equal(out6$p, c(.369), tolerance = 1e-2)
    out7 <- DIF(res.model, c('d'), Wald = T, items2test = 6, groups2test = c('g1', 'g2', 'g3'))
    expect_equal(out7$p, c(.371), tolerance = 1e-2)



})


