# PDIF example


#################
# Case 1:
#  Differences in latent trait means, but no DIF.
#  When PIRT model formed, no DIF should be present

N <- 5000

# simulate data
a <- matrix(c(
    1,0.5,NA,
    1,0.5,NA,
    1,0.5,NA,
    1,0.5,NA,
    1,NA,0.5,
    1,NA,0.5,
    1,NA,0.5,
    1,NA,0.5),ncol=3,byrow=TRUE)

d <- rnorm(nrow(a), sd=1.5)
a2 <- a
d2 <- d


sigma <- diag(3)
dat1 <- simdata(a, d, N, itemtype='2PL',sigma=sigma)
dat2 <- simdata(a2, d2, N, itemtype='2PL',sigma=sigma, mu = c(0, 0, 0))

dat <- rbind(dat1, dat2)
group <- rep(c('G1', 'G2'), each=N)

#############

# configuration invariance type model setup for maximal E-table fit
specific <- "S1 = 1-4
             S2 = 5-8"
bmod <- bfactor(dat, specific, group = group)
coef(bmod, simplify=TRUE)

# still configural
pmod <- pirt(bmod)
coef(pmod, simplify=TRUE)

# equated
pmod_e <- pirt(bmod, invariance=c(colnames(dat), 'free_mean', 'free_var'))
coef(pmod_e, simplify=TRUE)
anova(pmod, pmod_e)  # LR test shouldn't be optimal, AIC still good though

# DIF test (logic is backwards here though...)
for(i in 1:nrow(a)){
    pmod_enest <- pirt(bmod, invariance=c(colnames(dat)[-i], 'free_mean', 'free_var'))
    coef(pmod_enest, simplify=TRUE)
    print(anova(pmod_enest, pmod_e)) # LR tests will be wrong, AIC still good
    out <- nonnest2::vuongtest(pmod_e, pmod_enest, nested=TRUE)
    print(out$p_omega)  # this should be better
}
