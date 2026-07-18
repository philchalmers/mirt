expect_class <- function(x, class) expect_true(inherits(x, class))

test_that('PIRT', {
    data(SAT12)
    data <- key2binary(SAT12,
                       key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
    specific <- c(2,3,2,3,3,2,1,2,1,1,1,3,1,3,1,2,1,1,3,3,1,1,3,1,3,3,1,3,2,3,1,2)
    mod <- bfactor(data, specific)

    pmod <- pirt(mod)
    expect_equal(logLik(pmod), -9492.956, tolerance = 1e-3)

    eap <- fscores(pmod)
    eap2 <- fscores(mod, project=1)
    expect_equal(as.numeric(cor(eap, eap2)), .9999, tolerance=1e-4)

    FA <- pirt(mod, estimator = 'FA')
    expect_equal(as.numeric(FA[1,]), c(0.7602271, -1.0383243), tolerance=1e-3)
    LKA <- pirt(mod, estimator = 'LKA')
    expect_equal(as.numeric(LKA[2,]), c(1.3430096, 0.4245588), tolerance=1e-3)

    pirt_short <- pirt(mod, shortform=c('Item.1', 'Item.5', 'Item.10'))
    expect_equal(as.numeric(coef(pirt_short, simplify=TRUE)$items[,1]),
                 c(0.7664865, 0.9071186, 0.9228290), tolerance=1e-3)
    pirt_short2 <- pirt(mod, shortform=c(1, 5, 10))
    expect_true(identical(logLik(pirt_short), logLik(pirt_short2)))

})