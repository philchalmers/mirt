# Class "MixedClass"

Defines the object returned from
[`mixedmirt`](https://philchalmers.github.io/mirt/reference/mixedmirt.md).

## Slots

- `Call`::

  function call

- `Data`::

  list of data, sometimes in different forms

- `Options`::

  list of estimation options

- `Fit`::

  a list of fit information

- `Model`::

  a list of model-based information

- `ParObjects`::

  a list of the S4 objects used during estimation

- `OptimInfo`::

  a list of arguments from the optimization process

- `Internals`::

  a list of internal arguments for secondary computations (inspecting
  this object is generally not required)

- `vcov`::

  a matrix represented the asymptotic covariance matrix of the parameter
  estimates

- `time`::

  a data.frame indicating the breakdown of computation times in seconds

## Methods

- coef:

  `signature(object = "MixedClass")`

- print:

  `signature(x = "MixedClass")`

- residuals:

  `signature(object = "MixedClass")`

- show:

  `signature(object = "MixedClass")`

- summary:

  `signature(object = "MixedClass")`

- logLik:

  `signature(object = "MixedClass")`

- anova:

  `signature(object = "MixedClass")`

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## Author

Phil Chalmers <rphilip.chalmers@gmail.com>
