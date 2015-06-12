context('wideData')

test_that('wide1dim', {
    set.seed(1234)
    n <- 1000
    N  <- 1000
    dat <- simdata(matrix(1, n), matrix(rnorm(n)), N, itemtype = 'dich')
    one2n <- 1:n
    one2N <- 1:N
    for(i in 1:100000)
        dat[sample(one2N, 1), sample(one2n, 1)] <- NA
    mod <- mirt(dat, 1, verbose = FALSE)
    expect_equal(mod@logLik, -492047, tolerance = .01)
    vals <- mod2values(mod)$value
    vals <- vals[vals != 0 & vals != 1]
    expect_equal(fivenum(vals), c(-3.4079643, 0.1270154, 0.8914533, 1.0007773, 3.5769222),
                 tolerance = 1e-4)
    EAP <- fscores(mod, verbose=FALSE)
    expect_equal(head(EAP[,'F1']), c(-1.5404, -2.6672, -1.6762, -2.0182, -2.6940, -1.3094),
                 tolerance=1e-4)
})


