# Correlation matrix plot

Heatmap version of correlation matrix to visualize bivariate
correlations.

## Usage

``` r
plot_cors(
  r,
  numbers = TRUE,
  colors = TRUE,
  n = 51,
  main = NULL,
  zlim = c(-1, 1),
  show.legend = TRUE,
  labels = NULL,
  select = NULL,
  pval = NULL,
  cuts = c(0.001, 0.01),
  stars = FALSE,
  pdigits = 3,
  scale = TRUE,
  cex = 1,
  thresh = 0.2,
  min.cex = 0.6,
  drop.zero = TRUE,
  fit = TRUE,
  fill = 0.9,
  upper = TRUE,
  diag = TRUE,
  symmetric = TRUE,
  xaxis = 1,
  xlas = 90,
  gr = NULL,
  digits = 2,
  alpha = 0.75,
  border = "white",
  ...
)
```

## Arguments

- r:

  square matrix (correlations, or any matrix with values on a common
  scale), or a data.frame/non-square matrix of raw data, in which case
  `cor(r, use = "pairwise.complete.obs")` is computed

- numbers:

  print the (rounded) values in the cells

- colors:

  if FALSE, a grey scale is used

- n:

  number of colour bins

- main:

  title

- zlim:

  range of the colour scale (values are clamped to it)

- show.legend:

  draw the colour key

- labels:

  variable names (default: dimnames of r)

- select:

  indices/names of the variables to show

- pval:

  optional matrix of p-values (same dim as r). If supplied, the upper
  triangle shows the p-values as plain text with no colour fill, and the
  lower triangle/diagonal show the coloured correlations

- cuts:

  p-value cutpoints for stars, e.g. c(.001, .01) -\> \*\*\* / \*\* / \*
  (a final single \* is used for p \< .05)

- stars:

  optionally draw significance stars after the lower-triangle
  correlations (requires pval; FALSE by default, since the p-values
  themselves are shown in the upper triangle)

- pdigits:

  decimals for printed p-values (smaller values print as \<.001)

- scale:

  if TRUE, text is drawn at size cex for every cell (diagonal included)
  whose \|r\| \>= thresh, and shrinks linearly with \|r\| below thresh

- cex:

  A common shrink factor is set from the bulk of the labels (not the
  single widest one), then any label that still overflows its own cell
  is capped individually

- thresh:

  `|r|` below which the text starts to shrink

- min.cex:

  smallest relative text size (reached as \|r\| -\> 0)

- drop.zero:

  drop the leading zero of the printed values (.45, -.45)

- fit:

  if `TRUE`, text is shrunk (never enlarged) to fit the drawn cells; the
  size is measured at draw time, so it adapts to the device.

- fill:

  proportion of a cell width/height the typical label may fill

- upper:

  draw the upper triangle?

- diag:

  draw the diagonal?

- symmetric:

  if `TRUE` and `r` is square, row i is plotted in the same position as
  column i

- xaxis:

  1 = x labels at the bottom, 3 = at the top

- xlas:

  rotation (degrees) of x labels

- gr:

  colour-ramp function, e.g.,
  `colorRampPalette(c("red","white","blue"))`

- digits:

  digits shown in cells

- alpha:

  transparency of the cell colours

- border:

  colour of cell borders (NA/"transparent" for none)

- ...:

  passed on to
  [`lattice::levelplot()`](https://rdrr.io/pkg/lattice/man/levelplot.html)
