# Description of LSAT6 data

Data from Thissen (1982); contains 5 dichotomously scored items obtained
from the Law School Admissions Test, section 6.

## References

Thissen, D. (1982). Marginal maximum likelihood estimation for the
one-parameter logistic model. *Psychometrika, 47*, 175-186.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

# \donttest{
dat <- expand.table(LSAT6)
head(dat)
#>   Item_1 Item_2 Item_3 Item_4 Item_5
#> 1      0      0      0      0      0
#> 2      0      0      0      0      0
#> 3      0      0      0      0      0
#> 4      0      0      0      0      1
#> 5      0      0      0      0      1
#> 6      0      0      0      0      1
itemstats(dat)
#> $overall
#>     N mean_total.score sd_total.score ave.r sd.r alpha SEM.alpha
#>  1000            3.819          1.035 0.077 0.03 0.295     0.869
#> 
#> $itemstats
#>           N K  mean    sd total.r total.r_if_rm alpha_if_rm
#> Item_1 1000 2 0.924 0.265   0.362         0.113       0.275
#> Item_2 1000 2 0.709 0.454   0.567         0.153       0.238
#> Item_3 1000 2 0.553 0.497   0.618         0.173       0.217
#> Item_4 1000 2 0.763 0.425   0.534         0.144       0.246
#> Item_5 1000 2 0.870 0.336   0.435         0.122       0.266
#> 
#> $proportions
#>            0     1
#> Item_1 0.076 0.924
#> Item_2 0.291 0.709
#> Item_3 0.447 0.553
#> Item_4 0.237 0.763
#> Item_5 0.130 0.870
#> 

model <- 'F = 1-5
         CONSTRAIN = (1-5, a1)'
(mod <- mirt(dat, model))
#> 
#> Call:
#> mirt(data = dat, model = model)
#> 
#> Full-information item factor analysis with 1 factor(s).
#> Converged within 1e-04 tolerance after 12 EM iterations.
#> mirt version: 1.46.4 
#> M-step optimizer: BFGS 
#> EM acceleration: Ramsay 
#> Number of rectangular quadrature: 61
#> Latent density type: Gaussian 
#> 
#> Log-likelihood = -2466.938
#> Estimated parameters: 6 
#> AIC = 4945.875
#> BIC = 4975.322; SABIC = 4956.265
#> G2 (25) = 21.8, p = 0.6474
#> RMSEA = 0, CFI = NaN, TLI = NaN
M2(mod)
#>          M2 df     p RMSEA RMSEA_5 RMSEA_95 SRMSR   TLI CFI
#> stats 5.293  9 0.808     0       0    0.023 0.022 1.073   1
itemfit(mod)
#> Called from: EAPsum(x, S_X2 = TRUE, gp = gp, CUSTOM.IND = x@Internals$CUSTOM.IND, 
#>     den_fun = mirt_dmvnorm, quadpts = quadpts, theta_lim = theta_lim, 
#>     discrete = discrete, QMC = QMC, mixture = mixture, pis = pis, 
#>     which.items = which.items, use_dentype_estimate = use_dentype_estimate)
#> debug: if (version2) {
#>     if (length(CUSTOM.IND)) 
#>         stop("Custom items not yet supported for EAPsum_2.0", 
#>             call. = FALSE)
#>     for (i in seq_len(nspec)) {
#>         pick <- blist$specific == i
#>         if (i == 1) 
#>             pick <- blist$specific == i | is.na(blist$specific)
#>         tmpitemloc <- c(1, cumsum(K[pick]) + 1)
#>         itemtrace <- computeItemtrace(pars = pars[c(which(pick), 
#>             length(pars))], Theta = Theta, itemloc = tmpitemloc, 
#>             CUSTOM.IND = CUSTOM.IND, pis = pis)
#>         item_weights_long <- rep(item_weights[pick], K[pick])
#>         itemtrace <- t(itemtrace)^item_weights_long
#>         tmp <- calcL1(itemtrace = itemtrace, K = K[pick], itemloc = tmpitemloc)
#>         L1 <- t(tmp$L1)
#>         stage2K[i] <- length(tmp$Sum.Scores)
#>         subL1 <- matrix(0, ncol(L1), length(theta))
#>         for (j in 1:length(theta)) subL1[, j] <- colSums(L1[Theta[, 
#>             1] == theta[j], ] * sprior)
#>         L1_lst[[i]] <- subL1
#>     }
#>     itemtrace <- do.call(rbind, L1_lst)
#>     K <- stage2K
#>     itemloc <- c(1, cumsum(K) + 1)
#>     tmp <- calcL1(itemtrace = itemtrace, K = K, itemloc = itemloc)
#>     L1 <- tmp$L1
#>     Sum.Scores <- tmp$Sum.Scores
#>     Theta <- ThetaShort <- matrix(theta)
#>     prior <- den_fun(Theta, mean = gp$gmeans[1], sigma = gp$gcov[1, 
#>         1], ...)
#>     prior <- prior/sum(prior)
#>     nfact <- 1
#> } else {
#>     itemtrace <- computeItemtrace(pars = pars, Theta = Theta, 
#>         itemloc = itemloc, CUSTOM.IND = CUSTOM.IND, pis = pis)
#>     item_weights_long <- rep(item_weights, K)
#>     itemtrace <- t(itemtrace)^item_weights_long
#>     tmp <- calcL1(itemtrace = itemtrace, K = K, itemloc = itemloc)
#>     L1 <- tmp$L1
#>     Sum.Scores <- tmp$Sum.Scores
#> }
#> debug: itemtrace <- computeItemtrace(pars = pars, Theta = Theta, itemloc = itemloc, 
#>     CUSTOM.IND = CUSTOM.IND, pis = pis)
#> debug: item_weights_long <- rep(item_weights, K)
#> debug: itemtrace <- t(itemtrace)^item_weights_long
#> debug: tmp <- calcL1(itemtrace = itemtrace, K = K, itemloc = itemloc)
#> debug: L1 <- tmp$L1
#> debug: Sum.Scores <- tmp$Sum.Scores
#> debug: if (S_X2) {
#>     L1total <- L1 %*% prior
#>     Elist <- vector("list", J)
#>     for (i in which.items) {
#>         KK <- K[-i]
#>         T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>             , drop = FALSE]
#>         itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 
#>             1L)), , drop = FALSE]
#>         if (i != J) {
#>             itemloc2 <- itemloc[-i]
#>             itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#>         }
#>         else itemloc2 <- itemloc[-(J + 1)]
#>         tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#>         E <- matrix(NA, nrow(L1total), nrow(T))
#>         for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>             (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, 
#>             ]
#>         Elist[[i]] <- E[-c(1L, nrow(E)), ]
#>     }
#>     return(Elist)
#> }
#> debug: L1total <- L1 %*% prior
#> debug: Elist <- vector("list", J)
#> debug: for (i in which.items) {
#>     KK <- K[-i]
#>     T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#>     itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 
#>         1L)), , drop = FALSE]
#>     if (i != J) {
#>         itemloc2 <- itemloc[-i]
#>         itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#>     }
#>     else itemloc2 <- itemloc[-(J + 1)]
#>     tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#>     E <- matrix(NA, nrow(L1total), nrow(T))
#>     for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>         (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#>     Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> }
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-(J + 1)]
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: return(Elist)
#>     item  S_X2 df.S_X2 RMSEA.S_X2 p.S_X2
#> 1 Item_1 0.436       2          0  0.804
#> 2 Item_2 1.576       2          0  0.455
#> 3 Item_3 0.871       1          0  0.351
#> 4 Item_4 0.190       2          0  0.909
#> 5 Item_5 0.190       2          0  0.909
coef(mod, simplify=TRUE)
#> $items
#>           a1     d g u
#> Item_1 0.755 2.730 0 1
#> Item_2 0.755 0.999 0 1
#> Item_3 0.755 0.240 0 1
#> Item_4 0.755 1.307 0 1
#> Item_5 0.755 2.100 0 1
#> 
#> $means
#> F 
#> 0 
#> 
#> $cov
#>   F
#> F 1
#> 

# equivalentely, but with a different parameterization
mod2 <- mirt(dat, 1, itemtype = 'Rasch')
anova(mod, mod2) #equal
#>           AIC    SABIC       HQ      BIC    logLik X2 df   p
#> mod  4945.875 4956.265 4957.067 4975.322 -2466.938          
#> mod2 4945.875 4956.266 4957.067 4975.322 -2466.938  0  0 NaN
M2(mod2)
#>          M2 df     p RMSEA RMSEA_5 RMSEA_95 SRMSR   TLI CFI
#> stats 5.293  9 0.808     0       0    0.023 0.022 1.073   1
itemfit(mod2)
#> Called from: EAPsum(x, S_X2 = TRUE, gp = gp, CUSTOM.IND = x@Internals$CUSTOM.IND, 
#>     den_fun = mirt_dmvnorm, quadpts = quadpts, theta_lim = theta_lim, 
#>     discrete = discrete, QMC = QMC, mixture = mixture, pis = pis, 
#>     which.items = which.items, use_dentype_estimate = use_dentype_estimate)
#> debug: if (version2) {
#>     if (length(CUSTOM.IND)) 
#>         stop("Custom items not yet supported for EAPsum_2.0", 
#>             call. = FALSE)
#>     for (i in seq_len(nspec)) {
#>         pick <- blist$specific == i
#>         if (i == 1) 
#>             pick <- blist$specific == i | is.na(blist$specific)
#>         tmpitemloc <- c(1, cumsum(K[pick]) + 1)
#>         itemtrace <- computeItemtrace(pars = pars[c(which(pick), 
#>             length(pars))], Theta = Theta, itemloc = tmpitemloc, 
#>             CUSTOM.IND = CUSTOM.IND, pis = pis)
#>         item_weights_long <- rep(item_weights[pick], K[pick])
#>         itemtrace <- t(itemtrace)^item_weights_long
#>         tmp <- calcL1(itemtrace = itemtrace, K = K[pick], itemloc = tmpitemloc)
#>         L1 <- t(tmp$L1)
#>         stage2K[i] <- length(tmp$Sum.Scores)
#>         subL1 <- matrix(0, ncol(L1), length(theta))
#>         for (j in 1:length(theta)) subL1[, j] <- colSums(L1[Theta[, 
#>             1] == theta[j], ] * sprior)
#>         L1_lst[[i]] <- subL1
#>     }
#>     itemtrace <- do.call(rbind, L1_lst)
#>     K <- stage2K
#>     itemloc <- c(1, cumsum(K) + 1)
#>     tmp <- calcL1(itemtrace = itemtrace, K = K, itemloc = itemloc)
#>     L1 <- tmp$L1
#>     Sum.Scores <- tmp$Sum.Scores
#>     Theta <- ThetaShort <- matrix(theta)
#>     prior <- den_fun(Theta, mean = gp$gmeans[1], sigma = gp$gcov[1, 
#>         1], ...)
#>     prior <- prior/sum(prior)
#>     nfact <- 1
#> } else {
#>     itemtrace <- computeItemtrace(pars = pars, Theta = Theta, 
#>         itemloc = itemloc, CUSTOM.IND = CUSTOM.IND, pis = pis)
#>     item_weights_long <- rep(item_weights, K)
#>     itemtrace <- t(itemtrace)^item_weights_long
#>     tmp <- calcL1(itemtrace = itemtrace, K = K, itemloc = itemloc)
#>     L1 <- tmp$L1
#>     Sum.Scores <- tmp$Sum.Scores
#> }
#> debug: itemtrace <- computeItemtrace(pars = pars, Theta = Theta, itemloc = itemloc, 
#>     CUSTOM.IND = CUSTOM.IND, pis = pis)
#> debug: item_weights_long <- rep(item_weights, K)
#> debug: itemtrace <- t(itemtrace)^item_weights_long
#> debug: tmp <- calcL1(itemtrace = itemtrace, K = K, itemloc = itemloc)
#> debug: L1 <- tmp$L1
#> debug: Sum.Scores <- tmp$Sum.Scores
#> debug: if (S_X2) {
#>     L1total <- L1 %*% prior
#>     Elist <- vector("list", J)
#>     for (i in which.items) {
#>         KK <- K[-i]
#>         T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>             , drop = FALSE]
#>         itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 
#>             1L)), , drop = FALSE]
#>         if (i != J) {
#>             itemloc2 <- itemloc[-i]
#>             itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#>         }
#>         else itemloc2 <- itemloc[-(J + 1)]
#>         tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#>         E <- matrix(NA, nrow(L1total), nrow(T))
#>         for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>             (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, 
#>             ]
#>         Elist[[i]] <- E[-c(1L, nrow(E)), ]
#>     }
#>     return(Elist)
#> }
#> debug: L1total <- L1 %*% prior
#> debug: Elist <- vector("list", J)
#> debug: for (i in which.items) {
#>     KK <- K[-i]
#>     T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#>     itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 
#>         1L)), , drop = FALSE]
#>     if (i != J) {
#>         itemloc2 <- itemloc[-i]
#>         itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#>     }
#>     else itemloc2 <- itemloc[-(J + 1)]
#>     tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#>     E <- matrix(NA, nrow(L1total), nrow(T))
#>     for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>         (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#>     Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> }
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-i]
#> debug: itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: KK <- K[-i]
#> debug: T <- itemtrace[c(itemloc[i]:(itemloc[i + 1L] - 1L)), , drop = FALSE]
#> debug: itemtrace2 <- itemtrace[-c(itemloc[i]:(itemloc[i + 1L] - 1L)), 
#>     , drop = FALSE]
#> debug: if (i != J) {
#>     itemloc2 <- itemloc[-i]
#>     itemloc2[i:J] <- itemloc2[i:J] - nrow(T)
#> } else itemloc2 <- itemloc[-(J + 1)]
#> debug: itemloc2 <- itemloc[-(J + 1)]
#> debug: tmp <- calcL1(itemtrace = itemtrace2, K = KK, itemloc = itemloc2)
#> debug: E <- matrix(NA, nrow(L1total), nrow(T))
#> debug: for (j in 1L:(nrow(T))) E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% 
#>     (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: E[1L:nrow(tmp$L1) + j - 1L, j] <- tmp$L1 %*% (T[j, ] * prior)/L1total[1L:nrow(tmp$L1) + 
#>     j - 1L, ]
#> debug: Elist[[i]] <- E[-c(1L, nrow(E)), ]
#> debug: return(Elist)
#>     item  S_X2 df.S_X2 RMSEA.S_X2 p.S_X2
#> 1 Item_1 0.436       2          0  0.804
#> 2 Item_2 1.576       2          0  0.455
#> 3 Item_3 0.872       1          0  0.351
#> 4 Item_4 0.190       2          0  0.909
#> 5 Item_5 0.190       2          0  0.909
coef(mod2, simplify=TRUE)
#> $items
#>        a1     d g u
#> Item_1  1 2.731 0 1
#> Item_2  1 0.999 0 1
#> Item_3  1 0.240 0 1
#> Item_4  1 1.307 0 1
#> Item_5  1 2.100 0 1
#> 
#> $means
#> F1 
#>  0 
#> 
#> $cov
#>       F1
#> F1 0.572
#> 
sqrt(coef(mod2)$GroupPars[2]) #latent SD equal to the slope in mod
#> [1] 0.7561877

# }
```
