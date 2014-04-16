\name{DIF}
\alias{DIF}
\title{Differential item functioning tests}
\usage{
DIF(MGmodel, which.par, scheme = "add", items2test = 1:ncol(MGmodel@data),
  seq_stat = "SABIC", Wald = FALSE, p.adjust = "none", max_run = Inf,
  plotdif = FALSE, type = "trace", verbose = TRUE, ...)
}
\arguments{
  \item{MGmodel}{an object returned from
  \code{\link{multipleGroup}} to be used as the reference
  model}

  \item{which.par}{a character vector containing the
  parameter names which will be inspected for DIF}

  \item{Wald}{logical; perform Wald tests for DIF instead
  of likelihood ratio test?}

  \item{items2test}{a numeric vector, or character vector
  containing the item names, indicating which items will be
  tested for DIF. In models where anchor items are known,
  omit them from this vector. For example, if items 1 and 2
  are anchors in a 10 item test, then \code{items2test =
  3:10} would work for testing the remaining items
  (important to remember when using sequential schemes)}

  \item{scheme}{type of DIF analysis to perform, either by
  adding or dropping constraints across groups.  These can
  be: \describe{ \item{'add'}{parameters in
  \code{which.par} will be constrained each item one at a
  time for items that are specified in \code{items2test}.
  This is beneficial when examining DIF from a model with
  parameters freely estimated across groups, and when
  inspecting differences via the Wald test}
  \item{'drop'}{parameters in \code{which.par} will be
  freely estimated for items that are specified in
  \code{items2test}. This is useful when supplying an
  overly restrictive model and attempting to detect DIF
  with a slightly less restrictive model}
  \item{'add_sequential'}{sequentially loop over the items
  being tested, and at the end of the loop treat DIF tests
  that satisfy the \code{seq_stat} criteria as invariant.
  The loop is then re-run on the remaining invariant items
  to determine if they are now displaying DIF in the less
  constrained model, and when no new invariant item is
  found the algorithm stops and returns the items that
  displayed DIF} \item{'drop_sequential'}{sequentially loop
  over the items being tested, and at the end of the loop
  treat items that violate the \code{seq_stat} criteria as
  demonstrating DIF. The loop is then re-run, leaving the
  items that previously demonstrated DIF as variable across
  groups, and the remaining test items that previously
  showed invariance are re-tested. The algorithm stops when
  no more items showing DIF are found and returns the items
  that displayed DIF} }}

  \item{seq_stat}{select a statistic to test for in the
  sequential schemes. Potential values are (in descending
  order of power) \code{'AIC'}, \code{'AICc'},
  \code{'SABIC'}, and \code{'BIC'}.  If a numeric value is
  input that ranges between 0 and 1, the 'p' value will be
  tested (e.g., \code{seq_stat = .05} will test for the
  difference of p < .05 in the add scheme, or p > .05 in
  the drop scheme), along with the specified
  \code{p.adjust} input}

  \item{max_run}{a number indicating the maximum number of
  cycles to perform in sequential searches.  The default is
  to perform search until no further DIF is found}

  \item{plotdif}{logical; create itemplots for items that
  are displaying DIF according to the \code{seq_stat}
  criteria? Only available for 'add' type schemes}

  \item{type}{the \code{type} of plot argument passed to
  \code{plot()}. Default is 'trace', though another good
  option is 'infotrace'. For ease of viewing, the
  \code{facet_item} argument to mirt's \code{plot()}
  function is set to \code{TRUE}}

  \item{p.adjust}{string to be passed to the
  \code{\link{p.adjust}} function to adjust p-values.
  Adjustments are located in the \code{adj_pvals} element
  in the returned list}

  \item{verbose}{logical print extra information to the
  console?}

  \item{...}{additional arguments to be passed to
  \code{\link{multipleGroup}} and \code{plot}}
}
\description{
This function runs the Wald and likelihood-ratio approaches
for testing differential item functioning (DIF). This is
primarily a convenience wrapper to the
\code{\link{multipleGroup}} function for performing
standard DIF procedures. Models can be estimated in
parallel automatically by defining a parallel object with
\code{\link{mirtCluster}}, which will help to decrease the
runtime.
}
\examples{
\dontrun{

#simulate data where group 2 has a smaller slopes and more extreme intercepts
set.seed(12345)
a1 <- a2 <- matrix(abs(rnorm(15,1,.3)), ncol=1)
d1 <- d2 <- matrix(rnorm(15,0,.7),ncol=1)
a2[1:2, ] <- a1[1:2, ]/3
d1[c(1,3), ] <- d2[c(1,3), ]/4
head(data.frame(a.group1 = a1, a.group2 = a2, d.group1 = d1, d.group2 = d2))
itemtype <- rep('dich', nrow(a1))
N <- 1000
dataset1 <- simdata(a1, d1, N, itemtype)
dataset2 <- simdata(a2, d2, N, itemtype, mu = .1, sigma = matrix(1.5))
dat <- rbind(dataset1, dataset2)
group <- c(rep('D1', N), rep('D2', N))

#### no anchors, all items tested for DIF by adding item constrains one item at a time.
# define a parallel cluster (optional) to help speed up internal functions
mirtCluster()

#  Information matrix with S-EM
model <- multipleGroup(dat, 1, group, SE = TRUE)

#test whether adding slopes and intercepts constraints results in DIF. Plot items showing DIF
resulta1d <- DIF(model, c('a1', 'd'), plotdif = TRUE)
resulta1d

#same as above, but using Wald tests with Benjamini & Hochberg adjustment
resulta1dWald <- DIF(model, c('a1', 'd'), Wald = TRUE, p.adjust = 'fdr')
resulta1dWald
round(resulta1dWald$adj_pvals, 4)

#test whether adding only slope constraints results in DIF for all items
resulta1 <- DIF(model, 'a1')
resulta1

#following up on resulta1d, to determine whether it's a1 or d parameter causing DIF
(a1s <- DIF(model, 'a1', items2test = 1:3))
(ds <- DIF(model, 'd', items2test = 1:3))

#### using items 4 to 15 as anchors
itemnames <- colnames(dat)
model_anchor <- multipleGroup(dat, model = 1, group = group,
  invariance = c(itemnames[4:15], 'free_means', 'free_var'))
anchor <- DIF(model_anchor, c('a1', 'd'), items2test = 1:3)
anchor

### drop down approach (freely estimating parameters accross groups) when
### specifying a highly constrained model with estimated latent parameters
model_constrained <- multipleGroup(dat, 1, group,
  invariance = c(colnames(dat), 'free_means', 'free_var'))
dropdown <- DIF(model_constrained, 'd', scheme = 'drop')
dropdown

### sequential searches using SABIC as the selection criteria
# starting from completely different models
model <- multipleGroup(dat, 1, group)
stepup <- DIF(model, c('a1', 'd'), scheme = 'add_sequential')
stepup

#step down procedure for highly constrained model
model <- multipleGroup(dat, 1, group, invariance = itemnames)
stepdown <- DIF(model, c('a1', 'd'), scheme = 'drop_sequential')
stepdown

}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\seealso{
\code{\link{multipleGroup}}
}
\keyword{DIF}

