# Print generic for customized matrix console output

Provides a nicer output for most printed `matrix` objects defined by
functions in `mirt`.

## Usage

``` r
# S3 method for class 'mirt_matrix'
print(x, digits = 3, ...)

# S3 method for class 'matrix'
print(x, ...)
```

## Arguments

- x:

  object of class `'mirt_matrix'`

- digits:

  number of digits to round

- ...:

  additional arguments passed to `print(...)`
