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
    expect_equal(mod@logLik, -492030.2, tolerance = .01)
    vals <- mod2values(mod)$value
    vals <- vals[vals != 0 & vals != 1]
    expect_equal(fivenum(vals), c(-3.2911582, 0.2112903, 0.8786638, 0.9858656, 3.6639469),
                 tolerance = 1e-4)
    EAP <- fscores(mod, verbose=FALSE)
    expect_equal(head(EAP[,'F1']), c(-1.633319, -2.799962, -1.797049, -2.146386, -2.826118, -1.409135),
                 tolerance=1e-4)
})


