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
    expect_equal(extract.mirt(mod, 'logLik'), -491985.9, tolerance = .01)
    vals <- mod2values(mod)$value
    vals <- vals[vals != 0 & vals != 1]
    expect_equal(fivenum(vals), c(-3.50602785, 0.03518035, 0.90189515, 1.01846863, 3.46084417),
                 tolerance = 1e-4)
    EAP <- fscores(mod, verbose=FALSE, full.scores=FALSE)
    expect_equal(head(EAP[,'F1']), c(-1.4040, -2.5166, -1.5617, -1.8788, -2.5433, -1.1932),
                 tolerance=1e-4)
})


