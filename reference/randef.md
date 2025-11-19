# Compute posterior estimates of random effect

Stochastically compute random effects for `MixedClass` objects with
Metropolis-Hastings samplers and averaging over the draws to obtain
expected a posteriori predictions. Returns a list of the estimated
effects.

## Usage

``` r
randef(x, ndraws = 1000, thin = 10, return.draws = FALSE)
```

## Arguments

- x:

  an estimated model object from the
  [`mixedmirt`](https://philchalmers.github.io/mirt/reference/mixedmirt.md)
  function

- ndraws:

  total number of draws to perform. Default is 1000

- thin:

  amount of thinning to apply. Default is to use every 10th draw

- return.draws:

  logical; return a list containing the thinned draws of the posterior?

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29.

Chalmers, R. P. (2015). Extended Mixed-Effects Item Response Models with
the MH-RM Algorithm. *Journal of Educational Measurement, 52*, 200-222.
[doi:10.1111/jedm.12072](https://doi.org/10.1111/jedm.12072)
[doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{
# make an arbitrary groups
covdat <- data.frame(group = rep(paste0('group', 1:49), each=nrow(Science)/49))

# partial credit model
mod <- mixedmirt(Science, covdat, model=1, random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1915, Max-Change = 0.1573, Max-Change = 0.1007, Max-Change = 0.0881, Max-Change = 0.0682, Max-Change = 0.0656, Max-Change = 0.0356, Max-Change = 0.0255, Max-Change = 0.0751, Max-Change = 0.0268, Max-Change = 0.0304, Max-Change = 0.0229, Max-Change = 0.0200, Max-Change = 0.0445, Max-Change = 0.0308, Max-Change = 0.0444, Max-Change = 0.0135, Max-Change = 0.0410, Max-Change = 0.0253, Max-Change = 0.0290, Max-Change = 0.0157, Max-Change = 0.0171, Max-Change = 0.0140, Max-Change = 0.0117, Max-Change = 0.0350, Max-Change = 0.0598, Max-Change = 0.0187, Max-Change = 0.0132, Max-Change = 0.0261, Max-Change = 0.0297, Max-Change = 0.0074, Max-Change = 0.0052, Max-Change = 0.0281, Max-Change = 0.0186, Max-Change = 0.0302, Max-Change = 0.0417, Max-Change = 0.0215, Max-Change = 0.0501, Max-Change = 0.0178, Max-Change = 0.0388, Max-Change = 0.0248, Max-Change = 0.0528, Max-Change = 0.0227, Max-Change = 0.0224, Max-Change = 0.0226, Max-Change = 0.0182, Max-Change = 0.0037, Max-Change = 0.0091, Max-Change = 0.0254, Max-Change = 0.0213, Max-Change = 0.0231, Max-Change = 0.0188, Max-Change = 0.0261, Max-Change = 0.0165, Max-Change = 0.0136, Max-Change = 0.0445, Max-Change = 0.0255, Max-Change = 0.0269, Max-Change = 0.0198, Max-Change = 0.0410, Max-Change = 0.0160, Max-Change = 0.0145, Max-Change = 0.0096, Max-Change = 0.0160, Max-Change = 0.0495, Max-Change = 0.0382, Max-Change = 0.0194, Max-Change = 0.0273, Max-Change = 0.0555, Max-Change = 0.0118, Max-Change = 0.0238, Max-Change = 0.0052, Max-Change = 0.0371, Max-Change = 0.0250, Max-Change = 0.0414, Max-Change = 0.0157, Max-Change = 0.0294, Max-Change = 0.0270, Max-Change = 0.0228, Max-Change = 0.0059, Max-Change = 0.0323, Max-Change = 0.0277, Max-Change = 0.0092, Max-Change = 0.0094, Max-Change = 0.0214, Max-Change = 0.0300, Max-Change = 0.0302, Max-Change = 0.0401, Max-Change = 0.0091, Max-Change = 0.0820, Max-Change = 0.0138, Max-Change = 0.0163, Max-Change = 0.0619, Max-Change = 0.0539, Max-Change = 0.0059, Max-Change = 0.0250, Max-Change = 0.0126, Max-Change = 0.0472, Max-Change = 0.0278, Max-Change = 0.0999, Max-Change = 0.0248, Max-Change = 0.0933, Max-Change = 0.2000, Max-Change = 0.0784, Max-Change = 0.0183, Max-Change = 0.0312, Max-Change = 0.0170, Max-Change = 0.0667, Max-Change = 0.0366, Max-Change = 0.1640, Max-Change = 0.0398, Max-Change = 0.0143, Max-Change = 0.0223, Max-Change = 0.0374, Max-Change = 0.0143, Max-Change = 0.0153, Max-Change = 0.0348, Max-Change = 0.0118, Max-Change = 0.0261, Max-Change = 0.1477, Max-Change = 0.0443, Max-Change = 0.0091, Max-Change = 0.0144, Max-Change = 0.0168, Max-Change = 0.0444, Max-Change = 0.0190, Max-Change = 0.0157, Max-Change = 0.0231, Max-Change = 0.0213, Max-Change = 0.0231, Max-Change = 0.0278, Max-Change = 0.0082, Max-Change = 0.0440, Max-Change = 0.0320, Max-Change = 0.0238, Max-Change = 0.0408, Max-Change = 0.0082, Max-Change = 0.0251, Max-Change = 0.0223, Max-Change = 0.2000, Max-Change = 0.0111, Max-Change = 0.0245, Max-Change = 0.0052, Max-Change = 0.0256, Max-Change = 0.0124, Max-Change = 0.0316, Max-Change = 0.0363, Max-Change = 0.0302, Max-Change = 0.0088, Max-Change = 0.0356, Max-Change = 0.0298, Max-Change = 0.0329, Max-Change = 0.0517, Max-Change = 0.0219, Max-Change = 0.0085, Max-Change = 0.0108, Max-Change = 0.0283, Max-Change = 0.0234, Max-Change = 0.0235, Max-Change = 0.0381, Max-Change = 0.0311, Max-Change = 0.0693, Max-Change = 0.0397, Max-Change = 0.0596, Max-Change = 0.0471, Max-Change = 0.0343, Max-Change = 0.0371, Max-Change = 0.0371, Max-Change = 0.0160, Max-Change = 0.0327, Max-Change = 0.0152, Max-Change = 0.0242, Max-Change = 0.0370, Max-Change = 0.0250, Max-Change = 0.0365, Max-Change = 0.0145, Max-Change = 0.0242, Max-Change = 0.0182, Max-Change = 0.0203, Max-Change = 0.0343, Max-Change = 0.0251, Max-Change = 0.0231, Max-Change = 0.0061, Max-Change = 0.0145, Max-Change = 0.0328, Max-Change = 0.0212, Max-Change = 0.0181, Max-Change = 0.0260, Max-Change = 0.0101, Max-Change = 0.0139, Max-Change = 0.0107, Max-Change = 0.0053, Max-Change = 0.0552, Max-Change = 0.0186, Max-Change = 0.0305, Max-Change = 0.0174, Max-Change = 0.0232, Max-Change = 0.0176, Max-Change = 0.0055, Max-Change = 0.0257, Max-Change = 0.0202, Max-Change = 0.0191, Max-Change = 0.0085, Max-Change = 0.0236, Max-Change = 0.0117, Max-Change = 0.0123, Max-Change = 0.0412, Max-Change = 0.0324, Max-Change = 0.0294, Max-Change = 0.0184, Max-Change = 0.0159, Max-Change = 0.0174, Max-Change = 0.0132, Max-Change = 0.0290, Max-Change = 0.0222, Max-Change = 0.0331, Max-Change = 0.0446, Max-Change = 0.0547, Max-Change = 0.0381, Max-Change = 0.0320, Max-Change = 0.0527, Max-Change = 0.0478, Max-Change = 0.0351, Max-Change = 0.0380, Max-Change = 0.0211, Max-Change = 0.0213, Max-Change = 0.0547, Max-Change = 0.0155, Max-Change = 0.0092, Max-Change = 0.0178, Max-Change = 0.0393, Max-Change = 0.0074, Max-Change = 0.0178, Max-Change = 0.0199, Max-Change = 0.0740, Max-Change = 0.0120, Max-Change = 0.0230, Max-Change = 0.0285, Max-Change = 0.0158, Max-Change = 0.0081, Max-Change = 0.0257, Max-Change = 0.0076, Max-Change = 0.0203, Max-Change = 0.0119, Max-Change = 0.0317, Max-Change = 0.0153, Max-Change = 0.0327, Max-Change = 0.0362, Max-Change = 0.0118, Max-Change = 0.0187, Max-Change = 0.0571, Max-Change = 0.0288, Max-Change = 0.0133, Max-Change = 0.0365, Max-Change = 0.0157, Max-Change = 0.0135, Max-Change = 0.0211, Max-Change = 0.0368, Max-Change = 0.0104, Max-Change = 0.0265, Max-Change = 0.0138, Max-Change = 0.0203, Max-Change = 0.0366, Max-Change = 0.0227, Max-Change = 0.0356, Max-Change = 0.0329, Max-Change = 0.0566, Max-Change = 0.0124, Max-Change = 0.0127, Max-Change = 0.0163, Max-Change = 0.0133, Max-Change = 0.0307, Max-Change = 0.0247, Max-Change = 0.0069, Max-Change = 0.0221, Max-Change = 0.0158, Max-Change = 0.0242, Max-Change = 0.0052, Max-Change = 0.0228, Max-Change = 0.0170, Max-Change = 0.0107, Max-Change = 0.0293, Max-Change = 0.0277, Max-Change = 0.0145, Max-Change = 0.0471, Max-Change = 0.0053, Max-Change = 0.0236, Max-Change = 0.0284, Max-Change = 0.0147, Max-Change = 0.0189, Max-Change = 0.0147, Max-Change = 0.0253, Max-Change = 0.0429, Max-Change = 0.0081, Max-Change = 0.0276, Max-Change = 0.0364, Max-Change = 0.0184, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0217, gam = 0.1057, Max-Change = 0.0110, gam = 0.0780, Max-Change = 0.0088, gam = 0.0629, Max-Change = 0.0031, gam = 0.0532, Max-Change = 0.0044, gam = 0.0464, Max-Change = 0.0037, gam = 0.0413, Max-Change = 0.0045, gam = 0.0374, Max-Change = 0.0044, gam = 0.0342, Max-Change = 0.0029, gam = 0.0316, Max-Change = 0.0062, gam = 0.0294, Max-Change = 0.0045, gam = 0.0276, Max-Change = 0.0019, gam = 0.0260, Max-Change = 0.0010, gam = 0.0246, Max-Change = 0.0030, gam = 0.0233, Max-Change = 0.0011, gam = 0.0222, Max-Change = 0.0026, gam = 0.0212, Max-Change = 0.0026, gam = 0.0203, Max-Change = 0.0017, gam = 0.0195, Max-Change = 0.0023, gam = 0.0188, Max-Change = 0.0007, gam = 0.0181, Max-Change = 0.0036, gam = 0.0175, Max-Change = 0.0007, gam = 0.0169, Max-Change = 0.0025, gam = 0.0164, Max-Change = 0.0033, gam = 0.0159, Max-Change = 0.0020, gam = 0.0154, Max-Change = 0.0011, gam = 0.0150, Max-Change = 0.0030, gam = 0.0146, Max-Change = 0.0009, gam = 0.0142, Max-Change = 0.0022, gam = 0.0139, Max-Change = 0.0029, gam = 0.0135, Max-Change = 0.0014, gam = 0.0132, Max-Change = 0.0010, gam = 0.0129, Max-Change = 0.0010, gam = 0.0126, Max-Change = 0.0017, gam = 0.0124, Max-Change = 0.0010, gam = 0.0121, Max-Change = 0.0014, gam = 0.0119, Max-Change = 0.0025, gam = 0.0116, Max-Change = 0.0012, gam = 0.0114, Max-Change = 0.0010, gam = 0.0112, Max-Change = 0.0003, gam = 0.0110, Max-Change = 0.0009
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod)
#> 
#> Call:
#> mixedmirt(data = Science, covdata = covdat, model = 1, random = ~1 | 
#>     group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>      F1
#> F1 0.95
#> 
#> $group
#>           COV_group
#> COV_group    0.0185
#> 

effects <- randef(mod, ndraws = 2000, thin = 20)
head(effects$Theta)
#>               F1
#> [1,]  0.41948968
#> [2,]  0.06777752
#> [3,] -0.67178790
#> [4,] -0.60653538
#> [5,]  0.02488654
#> [6,]  0.82010861
head(effects$group)
#>                group
#> group1  1.798943e-02
#> group2 -7.928247e-03
#> group3  8.101692e-05
#> group4  2.228795e-02
#> group5  2.733886e-02
#> group6 -1.161409e-02

# lr.random input
mod2 <- mixedmirt(Science, covdat, model=1, lr.random = ~ 1|group)
#> , Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.2000, Max-Change = 0.1915, Max-Change = 0.1573, Max-Change = 0.1007, Max-Change = 0.0881, Max-Change = 0.0682, Max-Change = 0.0656, Max-Change = 0.0356, Max-Change = 0.0255, Max-Change = 0.0751, Max-Change = 0.0268, Max-Change = 0.0304, Max-Change = 0.0229, Max-Change = 0.0200, Max-Change = 0.0445, Max-Change = 0.0308, Max-Change = 0.0444, Max-Change = 0.0135, Max-Change = 0.0410, Max-Change = 0.0253, Max-Change = 0.0290, Max-Change = 0.0157, Max-Change = 0.0171, Max-Change = 0.0140, Max-Change = 0.0117, Max-Change = 0.0350, Max-Change = 0.0598, Max-Change = 0.0187, Max-Change = 0.0132, Max-Change = 0.0261, Max-Change = 0.0297, Max-Change = 0.0074, Max-Change = 0.0052, Max-Change = 0.0281, Max-Change = 0.0186, Max-Change = 0.0302, Max-Change = 0.0417, Max-Change = 0.0215, Max-Change = 0.0501, Max-Change = 0.0178, Max-Change = 0.0388, Max-Change = 0.0248, Max-Change = 0.0528, Max-Change = 0.0227, Max-Change = 0.0224, Max-Change = 0.0226, Max-Change = 0.0182, Max-Change = 0.0037, Max-Change = 0.0091, Max-Change = 0.0254, Max-Change = 0.0213, Max-Change = 0.0231, Max-Change = 0.0188, Max-Change = 0.0261, Max-Change = 0.0165, Max-Change = 0.0136, Max-Change = 0.0445, Max-Change = 0.0255, Max-Change = 0.0269, Max-Change = 0.0198, Max-Change = 0.0410, Max-Change = 0.0160, Max-Change = 0.0145, Max-Change = 0.0096, Max-Change = 0.0160, Max-Change = 0.0495, Max-Change = 0.0382, Max-Change = 0.0194, Max-Change = 0.0273, Max-Change = 0.0555, Max-Change = 0.0118, Max-Change = 0.0238, Max-Change = 0.0052, Max-Change = 0.0371, Max-Change = 0.0250, Max-Change = 0.0414, Max-Change = 0.0157, Max-Change = 0.0294, Max-Change = 0.0270, Max-Change = 0.0228, Max-Change = 0.0059, Max-Change = 0.0323, Max-Change = 0.0277, Max-Change = 0.0092, Max-Change = 0.0094, Max-Change = 0.0214, Max-Change = 0.0300, Max-Change = 0.0302, Max-Change = 0.0401, Max-Change = 0.0091, Max-Change = 0.0820, Max-Change = 0.0138, Max-Change = 0.0163, Max-Change = 0.0619, Max-Change = 0.0539, Max-Change = 0.0059, Max-Change = 0.0250, Max-Change = 0.1899, Max-Change = 0.1933, Max-Change = 0.1280, Max-Change = 0.1086, Max-Change = 0.0546, Max-Change = 0.0750, Max-Change = 0.1542, Max-Change = 0.1966, Max-Change = 0.2000, Max-Change = 0.0711, Max-Change = 0.0324, Max-Change = 0.0147, Max-Change = 0.0462, Max-Change = 0.0279, Max-Change = 0.0309, Max-Change = 0.0194, Max-Change = 0.0733, Max-Change = 0.1374, Max-Change = 0.0180, Max-Change = 0.0347, Max-Change = 0.1197, Max-Change = 0.0132, Max-Change = 0.0149, Max-Change = 0.0380, Max-Change = 0.0629, Max-Change = 0.0111, Max-Change = 0.0315, Max-Change = 0.0404, Max-Change = 0.0238, Max-Change = 0.0206, Max-Change = 0.0161, Max-Change = 0.0327, Max-Change = 0.0365, Max-Change = 0.0163, Max-Change = 0.0540, Max-Change = 0.0262, Max-Change = 0.0420, Max-Change = 0.0280, Max-Change = 0.0219, Max-Change = 0.0658, Max-Change = 0.0245, Max-Change = 0.0234, Max-Change = 0.0276, Max-Change = 0.0224, Max-Change = 0.0245, Max-Change = 0.0237, Max-Change = 0.0184, Max-Change = 0.0427, Max-Change = 0.0290, Max-Change = 0.0313, Max-Change = 0.0370, Max-Change = 0.0373, Max-Change = 0.0347, Max-Change = 0.0522, Max-Change = 0.0189, Max-Change = 0.0250, Max-Change = 0.0195, Max-Change = 0.0284, Max-Change = 0.0366, Max-Change = 0.0176, Max-Change = 0.0289, Max-Change = 0.0262, Max-Change = 0.0490, Max-Change = 0.0336, Max-Change = 0.0446, Max-Change = 0.0192, Max-Change = 0.0195, Max-Change = 0.0348, Max-Change = 0.0537, Max-Change = 0.0329, Max-Change = 0.0317, Max-Change = 0.0398, Max-Change = 0.0264, Max-Change = 0.0277, Max-Change = 0.0317, Max-Change = 0.0294, Max-Change = 0.0340, Max-Change = 0.0406, Max-Change = 0.0338, Max-Change = 0.0178, Max-Change = 0.0254, Max-Change = 0.0089, Max-Change = 0.0878, Max-Change = 0.0468, Max-Change = 0.0329, Max-Change = 0.0246, Max-Change = 0.0071, Max-Change = 0.0205, Max-Change = 0.0308, Max-Change = 0.0178, Max-Change = 0.0359, Max-Change = 0.0225, Max-Change = 0.0100, Max-Change = 0.0261, Max-Change = 0.0228, Max-Change = 0.0223, Max-Change = 0.0326, Max-Change = 0.0343, Max-Change = 0.0368, Max-Change = 0.0338, Max-Change = 0.0122, Max-Change = 0.0081, Max-Change = 0.0109, Max-Change = 0.0453, Max-Change = 0.0249, Max-Change = 0.0423, Max-Change = 0.0133, Max-Change = 0.0187, Max-Change = 0.0289, Max-Change = 0.0234, Max-Change = 0.0390, Max-Change = 0.0520, Max-Change = 0.0172, Max-Change = 0.0222, Max-Change = 0.0316, Max-Change = 0.0104, Max-Change = 0.0097, Max-Change = 0.0311, Max-Change = 0.0140, Max-Change = 0.0150, Max-Change = 0.0247, Max-Change = 0.0106, Max-Change = 0.0261, Max-Change = 0.0283, Max-Change = 0.0419, Max-Change = 0.0289, Max-Change = 0.0121, Max-Change = 0.0154, Max-Change = 0.0470, Max-Change = 0.0236, Max-Change = 0.0899, Max-Change = 0.0252, Max-Change = 0.0080, Max-Change = 0.0071, Max-Change = 0.0292, Max-Change = 0.0153, Max-Change = 0.0314, Max-Change = 0.0378, Max-Change = 0.0465, Max-Change = 0.0288, Max-Change = 0.0138, Max-Change = 0.0226, Max-Change = 0.0142, Max-Change = 0.0498, Max-Change = 0.0157, Max-Change = 0.0099, Max-Change = 0.0156, Max-Change = 0.0520, Max-Change = 0.0221, Max-Change = 0.0362, Max-Change = 0.0481, gam = 0.0000, Max-Change = 0.0000, gam = 0.1778, Max-Change = 0.0163, gam = 0.1057, Max-Change = 0.0270, gam = 0.0780, Max-Change = 0.0130, gam = 0.0629, Max-Change = 0.0047, gam = 0.0532, Max-Change = 0.0046, gam = 0.0464, Max-Change = 0.0076, gam = 0.0413, Max-Change = 0.0048, gam = 0.0374, Max-Change = 0.0059, gam = 0.0342, Max-Change = 0.0076, gam = 0.0316, Max-Change = 0.0091, gam = 0.0294, Max-Change = 0.0013, gam = 0.0276, Max-Change = 0.0029, gam = 0.0260, Max-Change = 0.0009, gam = 0.0246, Max-Change = 0.0028, gam = 0.0233, Max-Change = 0.0025, gam = 0.0222, Max-Change = 0.0031, gam = 0.0212, Max-Change = 0.0024, gam = 0.0203, Max-Change = 0.0025, gam = 0.0195, Max-Change = 0.0061, gam = 0.0188, Max-Change = 0.0031, gam = 0.0181, Max-Change = 0.0017, gam = 0.0175, Max-Change = 0.0028, gam = 0.0169, Max-Change = 0.0032, gam = 0.0164, Max-Change = 0.0024, gam = 0.0159, Max-Change = 0.0015, gam = 0.0154, Max-Change = 0.0014, gam = 0.0150, Max-Change = 0.0011, gam = 0.0146, Max-Change = 0.0020, gam = 0.0142, Max-Change = 0.0023, gam = 0.0139, Max-Change = 0.0007, gam = 0.0135, Max-Change = 0.0029, gam = 0.0132, Max-Change = 0.0012, gam = 0.0129, Max-Change = 0.0026, gam = 0.0126, Max-Change = 0.0017, gam = 0.0124, Max-Change = 0.0004, gam = 0.0121, Max-Change = 0.0012, gam = 0.0119, Max-Change = 0.0013, gam = 0.0116, Max-Change = 0.0010, gam = 0.0114, Max-Change = 0.0014, gam = 0.0112, Max-Change = 0.0007, gam = 0.0110, Max-Change = 0.0017, gam = 0.0108, Max-Change = 0.0011, gam = 0.0106, Max-Change = 0.0008, gam = 0.0104, Max-Change = 0.0011, gam = 0.0102, Max-Change = 0.0021, gam = 0.0101, Max-Change = 0.0004, gam = 0.0099, Max-Change = 0.0005, gam = 0.0098, Max-Change = 0.0005
#> 
#> Calculating information matrix...
#> 
#> Calculating log-likelihood...
summary(mod2)
#> 
#> Call:
#> mixedmirt(data = Science, covdata = covdat, model = 1, lr.random = ~1 | 
#>     group)
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $Theta
#>       F1
#> F1 0.969
#> 
#> 
#> --------------
#> RANDOM EFFECT COVARIANCE(S):
#> Correlations on upper diagonal
#> 
#> $group
#>           COV_group
#> COV_group     0.284
#> 

effects <- randef(mod2, ndraws = 2000)
head(effects$Theta)
#>              F1
#> [1,]  0.5542568
#> [2,]  0.1301837
#> [3,] -0.5320751
#> [4,] -0.7467891
#> [5,]  0.1070476
#> [6,]  0.8924793
head(effects$group)
#>              group
#> group1  0.05385179
#> group2 -0.04401677
#> group3 -0.39288394
#> group4 -0.19423063
#> group5  0.14682569
#> group6 -0.36128757

# }
```
