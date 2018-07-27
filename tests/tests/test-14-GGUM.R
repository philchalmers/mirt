context('GGUM')

test_that('GGUM', {

    mod <- mirt(Science, 1, 'ggum', verbose=FALSE, SE=T, SE.type='Louis')
    cfs <- as.vector(coef(mod, simplify=TRUE)$items)
    expect_equal(cfs, c(0.822457,0.8183556,2.245513,0.6954606,3.496507,3.216669,2.692917,3.592072,6.848083,5.279511,4.780072,6.566885,6.49711,4.273604,3.666343,4.734535,1.797227,0.9683621,1.853346,1.750892), tolerance = 1e-4)
    expect_equal(logLik(mod), -1611.485, tolerance = 1e-4)
    fs <- fscores(mod)
    expect_equal(as.vector(fs[1:6]), c(0.3912522,0.03987477,-0.9090671,-0.9090671,0.6249104,0.6301196), tolerance = 1e-4)
    fit <- itemfit(mod)
    expect_equal(fit$p.S_X2, c(0.3571234,0.2627762,0.3239219,0.3752965), tolerance = 1e-4)
    expect_equal(extract.mirt(mod, 'condnum'), 7811.406, tolerance = 1e-4)

    mod2 <- mirt(Science, 2, 'ggum', TOL = .1, verbose=FALSE)
    cfs <- as.vector(coef(mod2, simplify=TRUE)$items)
    expect_equal(cfs, c(0.597573,1.612019,2.93154,0.4885592,2.191113,0.9542118,2.929572,1.92599,1.893386,0.551031,0.528257,1.740618,-0.06953769,-0.4813906,-0.4118102,-0.1116948,2.264143,1.516003,1.504583,1.908433,1.796519,1.051688,1.047285,1.127507,0.2279314,-0.07262869,0.37068,0.1424873), tolerance = 1e-4)
    expect_equal(logLik(mod2), -1609.908, tolerance = 1e-4)
    fs <- fscores(mod2)
    expect_equal(as.vector(fs[1:6]), c(0.3565461,0.2227541,-0.3700079,-0.3700079,0.5399751,0.292352), tolerance = 1e-4)

    mod3 <- mirt(Science, 1, 'ggum', optimizer='NR', pars=mod2values(mod), verbose=FALSE)
    cfs <- as.vector(coef(mod3, simplify=TRUE)$items)
    expect_equal(cfs, c(0.8239157,0.8183987,2.239854,0.6964775,3.47662,3.217268,2.800517,3.582907,6.824561,5.27999,4.888454,6.554693,6.473052,4.274147,3.774769,4.724148,1.778522,0.9691021,1.961784,1.743114), tolerance = 1e-4)
    expect_equal(logLik(mod3), -1611.484, tolerance = 1e-4)
    fs <- fscores(mod3)
    expect_equal(as.vector(fs[1:6]), c(0.3926422,0.04013544,-0.9082043,-0.9082043,0.6219709,0.6316096), tolerance = 1e-4)
    fit <- itemfit(mod3)
    expect_equal(fit$p.S_X2, c(0.3597134,0.2630669,0.324352,0.3750729), tolerance = 1e-4)

    # unequal categories
    dat <- Science
    dat[,1] <- ifelse(Science[,1] == 1, 2, Science[,1]) - 1
    mod4 <- mirt(dat, 1, 'ggum', verbose=FALSE, TOL = .1)
    expect_equal(logLik(mod4), -1612.524, tolerance = 1e-4)

})

