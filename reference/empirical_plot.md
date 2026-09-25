# Function to generate empirical unidimensional item and test plots

Given a dataset containing item responses this function will construct
empirical graphics using the observed responses, potentially conditioned
on the (reduced) total score. When individual item plots are requested
then the total score will be formed without the item of interest (i.e.,
the total score without that item).

## Usage

``` r
empirical_plot(
  data,
  which.items = NULL,
  type = "prop",
  smooth = FALSE,
  sort = FALSE,
  formula = resp ~ s(TS, k = 5),
  org.data = NULL,
  discrim.cut = 0.2,
  main = NULL,
  par.strip.text = list(cex = 0.7),
  par.settings = list(strip.background = list(col = "#9ECAE1"), strip.border = list(col =
    "black")),
  auto.key = list(space = "right", points = FALSE, lines = TRUE),
  ...
)
```

## Arguments

- data:

  a `data.frame` or `matrix` of item responses (see
  [`mirt`](https://philchalmers.github.io/mirt/reference/mirt.md) for
  typical input)

- which.items:

  a numeric vector indicating which items to plot in a faceted image
  plot. If NULL then empirical test plots will be constructed instead

- type:

  character vector specifying type of plot to draw, some of which change
  as a function of the `which.item` input.

  'prop' (default)

  :   cumulative total-score proportions. If `which.items` specified
      then item-level conditional proportions are instead plotted
      against the reduced total scores

  'hist'

  :   histogram of total scores

  'freq'

  :   item response frequencies (supports `which.item`)

  'discrim'

  :   reduced item-total correlations to visualize discrimination
      effects

  'difficulty'

  :   mean/proportion of each item

  'discrim_diff'

  :   reduced item-total correlation against item difficulty

  'boxplot'

  :   conditional boxplots of reduced total scores (supports
      `which.items`)

  'bubble'

  :   bivariate frequency bubble plots (requires that `which.items` is
      exactly of length two)

- smooth:

  logical; include a GAM smoother instead of the raw proportions?
  Default is FALSE

- sort:

  logical; when applicable, sort the items first (e.g., in
  discrimination plot)?

- formula:

  formula used for the GAM smoother

- org.data:

  identical to `data`, but contains the "unscored" response options
  (e.g., the original coding in a multiple-choice test). Used in various
  item-level plots for diagnostic purposes, such as in distractor
  analyses. This will also automatically add size and linetype changes
  to highlight the detected scored categories

- discrim.cut:

  horizontal cut-off line to use when `type = 'discrim'`. Default is .2
  (to omit, use `NA`)

- main:

  the main title for the plot. If NULL an internal default will be used

- par.strip.text:

  plotting argument passed to
  [`lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html)

- par.settings:

  plotting argument passed to
  [`lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html)

- auto.key:

  plotting argument passed to
  [`lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html)

- ...:

  additional arguments to be passed to
  [`lattice`](https://rdrr.io/pkg/lattice/man/Lattice.html) and
  [`coef()`](https://rdrr.io/r/stats/coef.html)

## Details

Note that some of these plot types should only be used for
unidimensional tests with monotonically increasing item response
functions. If monotonicity is not true for all items, however, then
these plots may serve as a visual diagnostic tool so long as the
majority of items are indeed monotonic.

## References

Chalmers, R. P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. *Journal of Statistical Software, 48*(6),
1-29. [doi:10.18637/jss.v048.i06](https://doi.org/10.18637/jss.v048.i06)

## See also

[`itemstats`](https://philchalmers.github.io/mirt/reference/itemstats.md),
[`itemplot`](https://philchalmers.github.io/mirt/reference/itemplot.md),
[`itemGAM`](https://philchalmers.github.io/mirt/reference/itemGAM.md)

## Examples

``` r

# \donttest{

SAT12[SAT12 == 8] <- NA
data <- key2binary(SAT12,
   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))

# test plot
empirical_plot(data)

empirical_plot(data, type = 'hist')

empirical_plot(data, type = 'hist', breaks=20)

empirical_plot(data, type = 'discrim')

empirical_plot(data, type = 'discrim', sort=TRUE)

empirical_plot(data, type = 'difficulty')

empirical_plot(data, type = 'difficulty', sort=TRUE)

empirical_plot(data, type = 'discrim_diff')

empirical_plot(data, type = 'freq')


# items 1, 2 and 5
empirical_plot(data, c(1, 2, 5), type = 'freq')

empirical_plot(data, c(1, 2, 5))

empirical_plot(data, c(1, 2, 5), smooth = TRUE)

empirical_plot(data, c(1, 2, 5), type = 'boxplot')

empirical_plot(data, c(1, 2), type = 'bubble')

empirical_plot(data, c(1, 5), type = 'bubble')


# replace weird looking items with unscored versions for diagnostics
empirical_plot(data, 32)

data2 <- data
data2[,32] <- SAT12[,32]
empirical_plot(data2, 32)

empirical_plot(data2, 32, smooth = TRUE)


# alternatively, distractor analyses using original dataset
empirical_plot(data, which.items=32, org.data=SAT12)

empirical_plot(data, which.items=32, org.data=SAT12, smooth=TRUE)

empirical_plot(data, which.items=1:12, org.data=SAT12, smooth=TRUE)

empirical_plot(data, which.items=13:32, org.data=SAT12, smooth=TRUE)



#################
# polytomous response data
empirical_plot(Science)

empirical_plot(Science, type = 'hist', breaks=20)

empirical_plot(Science, type = 'freq')

empirical_plot(Science, type = 'difficulty')

empirical_plot(Science, type = 'discrim_diff')

empirical_plot(Science, type = 'boxplot')


# item-level
empirical_plot(Science, c(1, 2), type = 'bubble')

empirical_plot(Science, c(1, 3), type = 'bubble')

empirical_plot(Science, which.items = 1:4, type = 'prop')

empirical_plot(Science, which.items = 1:4, type = 'prop', smooth=TRUE)


# last plot very similar to model-based approach (though conditioned
#   on reduced total scores rather than scaled latent trait)
mod <- mirt(Science)
plot(mod, type='trace')


# when missing values present
Science[1:3, 1] <- NA
Science[6:8, 2] <- NA
empirical_plot(Science, type = 'freq')



# }
```
