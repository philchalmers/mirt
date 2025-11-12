# Extract an item object from `mirt` objects

Extract the internal `mirt` objects from any estimated model.

## Usage

``` r
extract.item(x, item, group = NULL, drop.zeros = FALSE)
```

## Arguments

- x:

  mirt model of class `'SingleGroupClass'`, `'MultipleGroupClass'`, or
  `'MixtureClass'`

- item:

  a number or character signifying which item to extract

- group:

  which group the item should be extracted from (applies to
  `'MultipleGroupClass'` and `'MixtureClass'` only). Can be a numeric
  value or the name of the group to be extracted

- drop.zeros:

  logical; drop slope values that are numerically close to zero to
  reduce dimensionality? Useful in objects returned from
  [`bfactor`](https://philchalmers.github.io/mirt/reference/bfactor.md)
  or other confirmatory models that contain several zero slopes

## References

Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`extract.group`](https://philchalmers.github.io/mirt/reference/extract.group.md),
[`extract.mirt`](https://philchalmers.github.io/mirt/reference/extract.mirt.md)

## Examples

``` r
# \donttest{
mod <- mirt(Science, 1)
extr.1 <- extract.item(mod, 1)
# }
```
