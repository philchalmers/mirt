context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- confmirt(fulldata, 1, verbose = FALSE)
    expect_is(onefact, 'ConfirmatoryClass')
    L <- matrix(0, 1, 22)
    L[1, c(1,5,9,13,17)] <- 1
    L2 <- matrix(0, 2, 22)
    L2[1, 1] <- L2[2, 5] <- 1
    L2[1, 9] <- L2[2, 13] <- -1
    W1 <- wald(L, onefact)
    W2 <- wald(L2, onefact)
    expect_is(W1, 'wald')
    expect_is(W2, 'wald')
    fitonefact <- fitIndices(onefact)
    expect_is(fitonefact, 'list')
    twofact <- confmirt(fulldata, 2, verbose = FALSE)
    expect_is(twofact, 'ExploratoryClass')
    fulldata[1,1] <- fulldata[2,2] <- NA
    onefactmissing <- confmirt(fulldata, 1, verbose = FALSE)
    expect_is(onefactmissing, 'ConfirmatoryClass')
        
    fs1 <- fscores(onefact, verbose = FALSE)
    expect_is(fs1, 'matrix')
    fs2 <- fscores(twofact, verbose = FALSE)
    expect_is(fs2, 'matrix') 
    fs3 <- fscores(onefactmissing, verbose = FALSE)
    expect_is(fs3, 'matrix')
})
 
