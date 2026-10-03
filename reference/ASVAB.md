# Description of ASVAB data

Table of counts extracted from Mislvey (1985). Data the 16 possible
response patterns observed for four items from the arithmetic reasoning
test of the Armed Services Vocational Aptitude Battery (ASVAB), Form 8A,
from samples of white males and females and black males and females.

## References

Mislevy, R. J. (1985). Estimation of latent group effects. *Journal of
the American Statistical Association, 80*, 993-997.

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r

data(ASVAB)
datWM <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, White_Male)))
datWF <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, White_Female)))
datBM <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, Black_Male)))
datBF <- expand.table(subset(ASVAB, select=c(Item.1:Item.4, Black_Female)))

dat <- rbind(datWM, datWF, datBM, datBF)
sex <- rep(c("Male", "Female", "Male", "Female"),
  times=c(nrow(datWM), nrow(datWF), nrow(datBM), nrow(datBF))) |> factor()
color <- rep(c("White", "Black"),
  times=c(nrow(datWM) + nrow(datWF), nrow(datBM) + nrow(datBF))) |> factor()
group <- sex:color

itemstats(dat, group=group)
#> Error in eval(substitute(expr), data, enclos = parent.frame()): object 'sd_total' not found
```
