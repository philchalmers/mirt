context('GGUM')

test_that('GGUM', {

    mod <- mirt(Science, 1, 'ggum', verbose=FALSE, SE=T, SE.type='Louis')
    cfs <- as.vector(coef(mod, simplify=TRUE)$items)
    expect_equal(cfs, c(0.8226008,0.8183312,2.246684,0.6955606,3.47901,3.216899,2.703884,3.581693,6.826873,5.279502,4.790232,6.555915,6.479501,4.273672,3.677056,4.723793,1.778566,0.968468,1.864631,1.739607), tolerance = 1e-4)
    expect_equal(logLik(mod), -1611.485, tolerance = 1e-4)
    fs <- fscores(mod)
    expect_equal(as.vector(fs[1:6]), c(0.3913489,0.03988078,-0.9089989,-0.9089989,0.6247694,0.6301075), tolerance = 1e-4)
    fit <- itemfit(mod)
    expect_equal(fit$p.S_X2, c(0.3577653,0.2629,0.3241444,0.3750244), tolerance = 1e-4)
    expect_equal(extract.mirt(mod, 'condnum'), 8328.319, tolerance = 1e-4)

    mod2 <- mirt(Science, 2, 'ggum', TOL = .1, verbose=FALSE)
    cfs <- as.vector(coef(mod2, simplify=TRUE)$items)
    expect_equal(cfs, c(0.6107421,1.629748,2.887708,0.4994923,2.205725,0.9528065,2.889943,1.896228,1.83997,0.5465146,0.526968,1.679464,-0.06901954,-0.5129543,-0.4140688,-0.1171546,2.254381,1.518071,1.506309,1.902736,1.784922,1.05284,1.049038,1.125199,0.2263336,-0.06741916,0.3692222,0.1366159), tolerance = 1e-4)
    expect_equal(logLik(mod2), -1609.922, tolerance = 1e-4)
    fs <- fscores(mod2)
    expect_equal(as.vector(fs[1:6]), c(0.3527256,0.2187546,-0.3718252,-0.3718252,0.5330427,0.2897776), tolerance = 1e-4)

    mod3 <- mirt(Science, 1, 'ggum', optimizer='NR', pars=mod2values(mod), verbose=FALSE)
    cfs <- as.vector(coef(mod3, simplify=TRUE)$items)
    expect_equal(cfs, c(0.8238878,0.8183917,2.240166,0.6964504,3.476685,3.217257,2.800519,3.582939,6.824679,5.279991,4.8884,6.554798,6.4732,4.27415,3.774732,4.724218,1.778553,0.969088,1.961831,1.743103), tolerance = 1e-4)
    expect_equal(logLik(mod3), -1611.485, tolerance = 1e-4)
    fs <- fscores(mod3)
    expect_equal(as.vector(fs[1:6]), c(0.3926059,0.04012361,-0.9082465,-0.9082465,0.6220411,0.6315483), tolerance = 1e-4)
    fit <- itemfit(mod3)
    expect_equal(fit$p.S_X2, c(0.3596684,0.2630703,0.3243456,0.3750705), tolerance = 1e-4)

    # unequal categories
    dat <- Science
    dat[,1] <- ifelse(Science[,1] == 1, 2, Science[,1]) - 1
    mod4 <- mirt(dat, 1, 'ggum', verbose=FALSE, TOL = .1)
    expect_equal(logLik(mod4), -1612.524, tolerance = 1e-4)

})

