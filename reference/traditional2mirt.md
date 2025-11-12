# Convert traditional IRT metric into slope-intercept form used in mirt

This is a helper function for users who have previously available
traditional/classical IRT parameters and want to know the equivalent
slope-intercept translation used in `mirt`. Note that this function
assumes that the supplied models are unidimensional by definition (i.e.,
will have only one slope/discrimination) and in the logistic metric
(i.e., logistic-ogive scaling coefficient D=1). If there is no supported
slope-intercept transformation available then the original vector of
parameters will be returned by default.

## Usage

``` r
traditional2mirt(x, cls, ncat)
```

## Arguments

- x:

  a vector of parameters to transform

- cls:

  the class or itemtype of the supplied model

- ncat:

  the number of categories implied by the IRT model

## Value

a named vector of slope-intercept parameters (if supported)

## Details

Supported class transformations for the `cls` input are:

- Rasch, 2PL, 3PL, 3PLu, 4PL:

  Form must be: (discrimination, difficulty, lower-bound, upper-bound)

- graded:

  Form must be: (discrimination, difficulty 1, difficulty 2, ...,
  difficulty k-1)

- gpcm:

  Form must be: (discrimination, difficulty 1, difficulty 2, ...,
  difficulty k-1)

- nominal:

  Form must be: (discrimination 1, discrimination 2, ..., discrimination
  k, difficulty 1, difficulty 2, ..., difficulty k)

## Examples

``` r
# classical 3PL model
vec <- c(a=1.5, b=-1, g=.1, u=1)
slopeint <- traditional2mirt(vec, '3PL', ncat=2)
slopeint
#>  a1   d   g   u 
#> 1.5 1.5 0.1 1.0 

# classical graded model (four category)
vec <- c(a=1.5, b1=-1, b2=0, b3=1.5)
slopeint <- traditional2mirt(vec, 'graded', ncat=4)
slopeint
#>    a1    d1    d2    d3 
#>  1.50  1.50  0.00 -2.25 

# classical generalize partial credit model (four category)
vec <- c(a=1.5, b1=-1, b2=0, b3=1.5)
slopeint <- traditional2mirt(vec, 'gpcm', ncat=4)
slopeint
#>    a1   ak0   ak1   ak2   ak3    d0    d1    d2    d3 
#>  1.50  0.00  1.00  2.00  3.00  0.00  1.50  1.50 -0.75 

# classical nominal model (4 category)
vec <- c(a1=.5, a2 = -1, a3=1, a4=-.5, d1=1, d2=-1, d3=-.5, d4=.5)
slopeint <- traditional2mirt(vec, 'nominal', ncat=4)
slopeint
#>         a1        ak0        ak1        ak2        ak3         d0         d1 
#> -0.3333333  0.0000000  4.5000000 -1.5000000  3.0000000  0.0000000 -2.0000000 
#>         d2         d3 
#> -1.5000000 -0.5000000 

```
