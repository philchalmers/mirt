context('GGUM')

test_that('GGUM', {

    mod <- suppressWarnings(mirt(Science, 1, 'ggum', verbose=FALSE, TOL=1e-3, SE=TRUE))
    cfs <- as.vector(coef(mod, simplify=TRUE)$items)
    expect_equal(cfs, c(1.709598,0.9826335,2.98249,1.334567,0.1497495,0.7743619,0.5923796,0.1931122,2.908172,2.693359,2.540201,2.556675,2.349809,1.803254,1.610007,1.371799,-0.1264261,-0.8272219,0.2777304,-0.2542214), tolerance = 1e-3)
    expect_equal(logLik(mod), -1624.744, tolerance = 1e-4)
    fs <- fscores(mod)
    expect_equal(as.vector(fs[1:6]), c(0.3274676,0.1678313,-0.8463551,-0.8463551,0.9050643,0.2260045), tolerance = 1e-4)
    fit <- itemfit(mod)
    expect_equal(fit$p.S_X2, c(0.428228,0.3061431,0.2439318,0.2789786), tolerance = 1e-4)
    expect_equal(extract.mirt(mod, 'condnum'), 729.1741, tolerance = 1e-4)

    mod2 <- mirt(Science, 2, 'ggum', TOL = .1, verbose=FALSE)
    cfs <- as.vector(coef(mod2, simplify=TRUE)$items)
    expect_equal(cfs, c(0.3546753,1.463776,2.275139,0.7488499,1.964075,1.109024,2.353222,1.642901,3.491854,0.5043591,0.4489166,0.2560512,0.1463108,-0.357271,-0.05547154,0.1603908,2.57393,1.461018,1.514204,1.745221,2.074269,1.030846,1.049058,1.004315,0.271844,-0.07844877,0.3124729,0.02956934), tolerance = 1e-4)
    expect_equal(logLik(mod2), -1614.866, tolerance = 1e-4)
    fs <- fscores(mod2)
    expect_equal(as.vector(fs[1:6]), c(0.3408671,0.2480399,-0.2098457,-0.2098457,0.4436503,0.3223817), tolerance = 1e-4)

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
    mod4 <- mirt(dat, 1, 'ggum', verbose=FALSE, TOL = .1)
    expect_equal(logLik(mod4), -1612.458, tolerance = 1e-4)

})

