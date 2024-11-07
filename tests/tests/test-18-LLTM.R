context('LLTM')

test_that('LLTM', {
    set.seed(42)
    a <- matrix(rep(1,30))
    d <- rep(c(1,0, -1),each = 10)  # first easy, then medium, last difficult
    dat <- simdata(a, d, 1000, itemtype = '2PL')
    itemdesign <- data.frame(difficulty =
        factor(c(rep('easy', 10), rep('medium', 10), rep('hard', 10))))
    rownames(itemdesign) <- colnames(dat)

    lltm <- mirt(dat, itemtype = 'Rasch', SE=TRUE, verbose=FALSE,
        item.formula = ~ 0 + difficulty, itemdesign=itemdesign)
    expect_equal(as.numeric(coef(lltm, simplify=TRUE)$items[1:2,1]), c(0.9587247, 0.9587247), tol=1e-2)
    expect_equal(as.numeric(coef(lltm, printSE=TRUE)[[1]][2, 1:3]), c(0.03886556, 0.03898911, 0.03769205), tol=1e-2)
    expect_equal(extract.mirt(lltm, 'condnum'), 6.633144, tol=1e-2)

    # additional information for LLTM
    oo <- plot(lltm)
    expect_is(oo, 'trellis')
    ifit <- itemfit(lltm)
    expect_equal(ifit$S_X2[1:3], c(20.06072, 20.90161, 23.48163), tol=1e-2)
    eap <- fscores(lltm)
    expect_equal(eap[1:3], c(1.0141828, -0.2427042, -0.2427042))

    # using unconditional modeling for first four items
    itemdesign.sub <- itemdesign[5:nrow(itemdesign), , drop=FALSE]
    lltm.4 <- mirt(dat, itemtype = 'Rasch', SE=TRUE, verbose=FALSE,
        item.formula = ~ 0 + difficulty, itemdesign=itemdesign.sub)
    cfs <- coef(lltm.4, simplify=TRUE)$items
    expect_equal(as.vector(cfs[c(1,5), 1:3]), c(0, .9063842, numeric(4)), tol=1e-2)
    expect_equal(anova(lltm, lltm.4)$p[2], 0.04288353, tol=1e-2)

})

test_that('MLTM', {

})

