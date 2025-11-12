# Define a parallel cluster object to be used in internal functions

This function defines a object that is placed in a relevant internal
environment defined in mirt. Internal functions such as `calcLogLik`,
`fscores`, etc, will utilize this object automatically to capitalize on
parallel processing architecture. The object defined is a call from
[`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).
Note that if you are defining other parallel objects (for simulation
designs, for example) it is not recommended to define a mirtCluster.

## Usage

``` r
mirtCluster(spec, omp_threads, remove = FALSE, ...)
```

## Arguments

- spec:

  input that is passed to
  [`parallel::makeCluster()`](https://rdrr.io/r/parallel/makeCluster.html).
  If no input is given the maximum number of available local cores minus
  1 will be used. Setting this to NULL will skip a new definition
  (allows `omp_threads` to be used independently)

- omp_threads:

  number of OpenMP threads to use (currently applies to E-step
  computations only). Not used when argument input is missing

- remove:

  logical; remove previously defined `mirtCluster()`?

- ...:

  additional arguments to pass to
  [`parallel::makeCluster`](https://rdrr.io/r/parallel/makeCluster.html)

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>

## Examples

``` r
if (FALSE) { # \dontrun{
if(interactive()){
  # use all available cores
  mirtCluster()
  mirtCluster(remove = TRUE)

  # make 4 cores available for parallel computing
  mirtCluster(4)
  mirtCluster(remove = TRUE)

  # create 3 core architecture in R, and 4 thread architecture with OpenMP
  mirtCluster(spec = 3, omp_threads = 4)

  # leave previous multicore objects, but change omp_threads
  mirtCluster(spec = NULL, omp_threads = 2)
}

} # }
```
