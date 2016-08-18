context('gpcm_mats')

test_that('gpcm_mats', {

    #----------
    # gpcm vs. gpcm_mats with optimizer="NR" to indirectly test 2nd derivatives
    dat<-Science
    ni<-ncol(dat)
    ncat<-4

    # gpcm using built-in gpcm item model and EM
    gpcm.em<-mirt(dat,1,itemtype=rep("gpcm",ni),optimizer="NR", verbose=FALSE)

    # gpcm_mats with custom scoring functions
    sf.matrix<-matrix(0:(ncat-1),nrow=ncat,ncol=1)
    sf.list<-list()
    for(i in 1:ni){ sf.list[[i]]<-sf.matrix}
    gpcm.sf.em<-mirt(dat,1,itemtype=rep("gpcm",ni),gpcm_mats=sf.list, optimizer="NR", verbose=FALSE)

    expect_equal(extract.mirt(gpcm.em, 'logLik'), extract.mirt(gpcm.sf.em, 'logLik'), tolerance=1e-4)

    #----------
    # gpcm_mats with 2 custom scoring functions
    dat<-Science
    ni<-ncol(dat)
    ncat<-4

    sf.matrix<-matrix(c(
        0,1,2,3,
        1,0,0,1),
        nrow=4,ncol=2)
    sf.list<-list()
    for(i in 1:ni){ sf.list[[i]]<-sf.matrix}
    model<-"
     F1 = 1-4
     F2 = 1-4
     COV = F1*F2
    "
    model<-mirt.model(model)
    gpcm.sf.em<-mirt(dat,model,itemtype=rep("gpcm",ni),gpcm_mats=sf.list, optimizer="NR", quadpts=21, verbose=FALSE)

    expect_equal(-1595.52, extract.mirt(gpcm.sf.em, 'logLik'), tolerance=.01)
    est<-coef(gpcm.sf.em)
    expect_equal(est[[1]][1:2], c(.721, .663),tolerance = .01)

})
