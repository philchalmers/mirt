context('gpcm_mats')

test_that('gpcm_mats', {

    #----------
    # gpcm vs. gpcm_mats with optimizer="NR" to indirectly test 2nd derivatives
    dat<-Science
    ni<-ncol(dat)
    ncat<-4

    # gpcm using built-in gpcm item model and EM
    gpcm.em<-mirt(dat,1,itemtype=rep("gpcm",ni),optimizer="NR")

    # gpcm_mats with custom scoring functions
    sf.matrix<-matrix(0:(ncat-1),nrow=ncat,ncol=1)
    sf.list<-list()
    for(i in 1:ni){ sf.list[[i]]<-sf.matrix}
    gpcm.sf.em<-mirt(dat,1,itemtype=rep("gpcm",ni),gpcm_mats=sf.list, optimizer="NR")

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
    gpcm.sf.em<-mirt(dat,model,itemtype=rep("gpcm",ni),gpcm_mats=sf.list, optimizer="NR", quadpts=49)

    expect_equal(-1594.08, extract.mirt(gpcm.sf.em, 'logLik'), tolerance=.01)
    est<-coef(gpcm.sf.em)
    expect_equal(est[[1]][1:2], c(.81, .68),tolerance = .01)

}
