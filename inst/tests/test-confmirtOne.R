context('confmirtOne')

test_that('exploratory mods', {
    data(LSAT7)
    fulldata <- expand.table(LSAT7)
    onefact <- confmirt(fulldata, 1, verbose = FALSE)
    expect_is(onefact, 'ConfirmatoryClass')
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
 
