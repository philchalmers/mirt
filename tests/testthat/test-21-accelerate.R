test_that('EM acceleration rejects invalid extrapolations', {
    # graded items with closely spaced intercepts: unchecked Ramsay/SQUAREM steps
    # used to push intercepts out of order, giving a NaN log-likelihood
    set.seed(2)
    J <- 30
    a <- matrix(rlnorm(J, .2, .3))
    d <- t(apply(matrix(rnorm(J*3), J), 1, sort, decreasing=TRUE))
    dat <- simdata(a, d, 200, itemtype='graded')

    mod_none <- mirt(dat, 1, verbose=FALSE, accelerate='none')
    mod_ramsay <- mirt(dat, 1, verbose=FALSE)
    mod_squarem <- mirt(dat, 1, verbose=FALSE, accelerate='squarem')
    expect_true(is.finite(logLik(mod_ramsay)))
    expect_true(is.finite(logLik(mod_squarem)))
    expect_equal(logLik(mod_ramsay), logLik(mod_none), tolerance = 1e-6)
    expect_equal(logLik(mod_squarem), logLik(mod_none), tolerance = 1e-6)
})
