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
    expect_equal(extract.mirt(mod, 'logLik'), -491760, tolerance = .01)
    vals <- mod2values(mod)$value
    vals <- vals[vals != 0 & vals != 1]
    expect_equal(fivenum(vals), c(-3.3647728,  0.1301622,  0.8968993,  1.0117605,  3.7725929),
                 tolerance = 1e-4)
    EAP <- fscores(mod, verbose=FALSE, full.scores=FALSE)
    expect_equal(head(EAP[,'F1']), c(-1.5322, -2.6723, -1.6393, -2.0958, -2.6657, -2.2742),
                 tolerance=1e-4)
})


