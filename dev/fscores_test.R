library(mirt)

a <- matrix(c(1.2, 1, 0,0,
              1.2, 1, 0,0,
              1,0,.8,0,
              1,0,.8,0,
              .8,0,0,1.2,
              .8,0,0,1.2), nrow=6, byrow=TRUE)
d <- c(-1, -.6, -.2, .2, .6, 1)

dat <- simdata(a, d, itemtype='2PL', N=10000)
itemstats(dat)

sv <- bfactor(dat, model=c(1,1,2,2,3,3), pars='values')
sv$value[grepl('^a', sv$name)] <- as.vector(t(a))
sv$value[grepl('^d', sv$name)] <- d

umod <- mirt(dat)
mod <- bfactor(dat, model=c(1,1,2,2,3,3), pars=sv, TOL=NaN)
coef(mod, simplify=TRUE)$items
pmod <- pirt(mod)

#################

# compare UIRT to marginalized (too optimistic due to ignored LD)
fscores(umod, method='EAPsum', full.scores=FALSE)

# compare to PIRT (less LD issues, but still too optimistic)
fscores(pmod, method='EAPsum', full.scores=FALSE)

# marginalizes over the specific factors first (better)
fscores(mod, method='EAPsum_2.0', full.scores=FALSE)

###############

# should give same theta_g and SE(theta_g) for all 0, all 1 combos, matching
# the LW 2.0 approach
fscores(mod, full.scores=FALSE, method='EAP_general')[c(1,64), ]
fscores(mod, method='EAPsum_2.0', full.scores=FALSE)

# dev version of projection at the trait estimate level.
# should give *approximately* the same theta_g as PIRT model
# for all 0, all 1 combos, but slightly higher SE(theta_g) due to
# the internal marginalization. This is approximate because the PIRT model
# is itself only an approximation, so the internal marginalization is
# technically more correct
#
# Issue is that the SEs are still too small as the LD persists ...
fscores(mod, full.scores=FALSE, project = 1)[c(1,64), ]
fscores(pmod, method='EAPsum', full.scores=FALSE)


#########################################################
# bigger example

library(mirt)

# much stronger specific factors (should affect the SEs more)
a <- matrix(c(1.2, 2, 0,0,
              1.2, 2, 0,0,
              1.2, 2, 0,0,
              1,0, 2.5,0,
              1,0, 2.5,0,
              1,0, 2.5,0,
              .8,0,0,3,
              .8,0,0,3,
              .8,0,0,3), nrow=9, byrow=TRUE)
d <- c(-1.5, -1, -.6, -.2, 0, .2, .6, 1, 1.5)

dat <- simdata(a, d, itemtype='2PL', N=10000)
itemstats(dat)

sv <- bfactor(dat, model=c(1,1,1,2,2,2,3,3,3), pars='values')
sv$value[grepl('^a', sv$name)] <- as.vector(t(a))
sv$value[grepl('^d', sv$name)] <- d

umod <- mirt(dat)
mod <- bfactor(dat, model=c(1,1,1,2,2,2,3,3,3), pars=sv, TOL=NaN)
coef(mod, simplify=TRUE)$items

(FApprox <- pirt(mod, estimator = 'FA'))
pmod <- pirt(mod)
coef(pmod, simplify=TRUE)

# set to FA estimates
sv <- mod2values(pmod)
sv$value[grepl('^a', sv$name)] <- FApprox[,1]
sv$value[grepl('^d', sv$name)] <- FApprox[,2]
pmod2 <- pirt(mod, pars=sv, TOL=NaN)
coef(pmod2, simplify=TRUE)

#################

# compare UIRT to marginalized (too optimistic due to ignored LD)
fscores(umod, method='EAPsum', full.scores=FALSE)

# compare to PIRT (less LD issues, but still too optimistic)
fscores(pmod, method='EAPsum', full.scores=FALSE)
fscores(pmod2, method='EAPsum', full.scores=FALSE)

# marginalizes over the specific factors first (better)
fscores(mod, method='EAPsum_2.0', full.scores=FALSE)

###############

# should give same theta_g and SE(theta_g) for all 0, all 1 combos, matching
# the LW 2.0 approach
fs <- fscores(mod, full.scores=FALSE, method='EAP_general')
fs[c(1, nrow(fs)),  c('G', 'SE_G')]
fscores(mod, method='EAPsum_2.0', full.scores=FALSE)

# dev version of projection at the trait estimate level.
# should give same theta_g as PIRT model
# for all 0, all 1 combos, but slightly higher SE(theta_g) due to
# the internal marginalization
fsp <- fscores(mod, full.scores=FALSE, project = 1, quadpts=31)
fsp[c(1, nrow(fsp)),  c('G', 'SE_G')]
fscores(pmod, method='EAPsum', full.scores=FALSE)
fscores(pmod2, method='EAPsum', full.scores=FALSE)
