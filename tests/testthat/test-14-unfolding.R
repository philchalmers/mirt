expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('GGUM', {

    mod <- suppressWarnings(mirt(Science, 1, c('graded', 'graded', 'graded', 'ggum'),
                                 verbose=FALSE, SE=TRUE))
    cfs <- as.vector(coef(mod, simplify=TRUE)$items)
    expect_equal(cfs, c(1.049862,1.218629,2.280344,0.7421556,4.869889,2.917697,5.214802,NA,2.645195,0.8993259,2.206224,NA,-1.469685,-2.262823,-1.957192,NA,NA,NA,NA,3.477198,NA,NA,NA,6.342889,NA,NA,NA,4.566828,NA,NA,NA,1.719266), tolerance = 1e-3)
    expect_equal(logLik(mod), -1609.924, tolerance = 1e-4)
    fs <- fscores(mod)
    expect_equal(as.vector(fs[1:6]), c(0.4267278,0.04564096,-0.9137894,-0.9137894,0.5612607,0.6853446), tolerance = 1e-2)
    fit <- itemfit(mod)
    expect_equal(fit$p.S_X2, c(0.6012163,0.3808717,0.5043588,0.4040644), tolerance = 1e-2)
    # expect_equal(extract.mirt(mod, 'condnum'), 2745.397, tolerance = 1e-4)

    mod2 <- mirt(Science, 2,c('graded', 'graded', 'graded', 'ggum'), TOL=.01,
                 verbose=FALSE)
    cfs <- as.vector(coef(mod2, simplify=TRUE)$items)
    expect_equal(cfs, c(-1.618017,-0.9619207,-1.670646,1.04929,-0.2230827,1.414658,1.191751,1e-04,5.601737,3.359375,4.890594,NA,3.118574,1.04238,2.061643,NA,-1.756586,-2.625442,-1.828546,NA,NA,NA,NA,-3.295782,NA,NA,NA,-0.8652751,NA,NA,NA,5.7103,NA,NA,NA,4.16543,NA,NA,NA,1.901035), tolerance = 1e-1)
    expect_equal(logLik(mod2), -1602.618, tolerance = 1e-4)
    fs <- fscores(mod2)
    expect_equal(as.vector(fs[1:6]), c(0.3068925, -1.0081743, -1.6638898, -1.6638898,  0.6702061, -0.1884406), tolerance = 1e-1)

    # mod3 <- mirt(Science, 1, 'ggum', optimizer='NR', pars=mod2values(mod), verbose=FALSE)
    # cfs <- as.vector(coef(mod3, simplify=TRUE)$items)
    # expect_equal(cfs, c(0.8238878,0.8183917,2.240166,0.6964504,3.476685,3.217257,2.800519,3.582939,6.824679,5.279991,4.8884,6.554798,6.4732,4.27415,3.774732,4.724218,1.778553,0.969088,1.961831,1.743103), tolerance = 1e-4)
    # expect_equal(logLik(mod3), -1611.484, tolerance = 1e-4)
    # fs <- fscores(mod3)
    # expect_equal(as.vector(fs[1:6]), c(0.3926059,0.04012361,-0.9082465,-0.9082465,0.6220411,0.6315483), tolerance = 1e-4)
    # fit <- itemfit(mod3)
    # expect_equal(fit$p.S_X2, c(0.3596684,0.2630703,0.3243456,0.3750705), tolerance = 1e-4)

    # unequal categories
    dat <- Science
    dat[,1] <- ifelse(Science[,1] == 1, 2, Science[,1]) - 1
    mod4 <- mirt(dat, 1, c('ggum', 'graded', 'graded', 'graded'), verbose=FALSE)
    expect_equal(logLik(mod4), -1595.793, tolerance = 1e-4)

})

test_that('unfolding', {
    dat <- expand.table(LSAT6)

    # ideal
    mod <- mirt(dat, 1, 'ideal', verbose=FALSE)
    expect_equal(-2466.925, logLik(mod), tolerance = .01)

    # hcm
    mod <- mirt(dat, 1, c(rep('2PL',4), 'hcm'), SE=TRUE, verbose=FALSE)
    expect_equal(-2466.898, logLik(mod), tolerance = .01)
    expect_equal(305.8971, extract.mirt(mod, 'condnum'), tolerance=.1)
    mod <- mirt(dat, 1, c(rep('2PL',4), 'ghcm'), verbose=FALSE)
    expect_equal(-2466.654, logLik(mod), tolerance = .01)

    mod2a <- mirt(dat, 'F1 = 1-3
                        F2 = 3-5', c(rep('2PL',4), 'hcm'), verbose=FALSE, SE=TRUE)
    expect_equal(-2469.595, logLik(mod2a), tolerance = .01)
    expect_equal(724.7384, extract.mirt(mod2a, 'condnum'), tolerance=.1)
    mod2b <- mirt(dat, 2, c(rep('2PL',4), 'hcm'), verbose=FALSE)
    expect_equal(-2465.09, logLik(mod2b), tolerance = .01)

    mod <- mirt(Science, 1, c('hcm', rep('graded',3)), verbose=FALSE)
    expect_equal(-1611.122, logLik(mod), tolerance = .01)
    mod <- mirt(Science, 1, c(rep('graded',3), 'ghcm'), verbose=FALSE)
    expect_equal(-1610.266, logLik(mod), tolerance = .01)


})

