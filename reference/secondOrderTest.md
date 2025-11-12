# Second-order test of convergence

Test whether terminated estimation criteria for a given model passes the
second order test by checking the positive definiteness of the resulting
Hessian matrix. This function, which accepts the symmetric
Hessian/information matrix as the input, returns `TRUE` if the matrix is
positive definite and `FALSE` otherwise.

## Usage

``` r
secondOrderTest(mat, ..., method = "eigen")
```

## Arguments

- mat:

  symmetric matrix to test for positive definiteness (typically the
  Hessian at the highest point of model estimator, such as MLE or MAP)

- ...:

  arguments passed to either
  [`eigen`](https://rdrr.io/r/base/eigen.html),
  [`chol`](https://rdrr.io/r/base/chol.html), or `'det'` for the
  positiveness of the eigen values, positiveness of leading minors via
  the Cholesky decomposition, or evaluation of whether the determinant
  is greater than 0

- method:

  method to use to test positive definiteness. Default is `'eigen'`

## Value

a matrix with all possible combinations

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
# \donttest{

# PD matrix
mod <- mirt(Science, 1, SE=TRUE)
info <- solve(vcov(mod))   ## observed information
secondOrderTest(info)
#> [1] TRUE
secondOrderTest(info, method = 'chol')
#> [1] TRUE
secondOrderTest(info, method = 'det')
#> [1] TRUE

# non-PD matrix
mat <- matrix(c(1,0,0,0,1,1,0,1,1), ncol=3)
mat
#>      [,1] [,2] [,3]
#> [1,]    1    0    0
#> [2,]    0    1    1
#> [3,]    0    1    1
secondOrderTest(mat)
#> [1] FALSE
secondOrderTest(mat, method = 'chol')
#> [1] FALSE
secondOrderTest(mat, method = 'det')
#> [1] FALSE

# }
```
