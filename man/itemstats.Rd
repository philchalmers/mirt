% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/itemstats.R
\name{itemstats}
\alias{itemstats}
\title{Generic item summary statistics}
\usage{
itemstats(
  data,
  group = NULL,
  use_ts = TRUE,
  itemfreq = "proportions",
  ts.tables = FALSE
)
}
\arguments{
\item{data}{An object of class \code{data.frame} or \code{matrix}
with the response patterns}

\item{group}{optional grouping variable to condition on when computing
summary information}

\item{use_ts}{logical; include information that is conditional on a
meaningful total score?}

\item{itemfreq}{character vector indicting whether to
include item response \code{"proportions"} or \code{"counts"}
for each item. If set to \code{'none'} then this will be omitted}

\item{ts.tables}{logical; include mean/sd summary information
pertaining to the unweighted total score?}
}
\value{
Returns a list containing the summary statistics
}
\description{
Function to compute generic item summary statistics that do not require
prior fitting of IRT models. Contains information about coefficient alpha
(and alpha if an item is deleted), mean/SD and frequency of total scores,
reduced item-total correlations, average/sd of the correlation between items,
response frequencies, and conditional mean/sd information given the
unweighted sum scores. Summary information involving the total scores
only included for responses with no missing data to ensure the metric is
meaningful, however standardized statistics (e.g., correlations) utilize
all possible response information.
}
\examples{

# dichotomous data example
LSAT7full <- expand.table(LSAT7)
head(LSAT7full)
itemstats(LSAT7full)
itemstats(LSAT7full, itemfreq='counts')

# behaviour with missing data
LSAT7full[1:5,1] <- NA
itemstats(LSAT7full)

# data with no meaningful total score
head(SAT12)
itemstats(SAT12, use_ts=FALSE)

# extra total scores tables
dat <- key2binary(SAT12,
                   key = c(1,4,5,2,3,1,2,1,3,1,2,4,2,1,
                           5,3,4,4,1,4,3,3,4,1,3,5,1,3,1,5,4,5))
itemstats(dat, ts.tables=TRUE)

# grouping information
group <- gl(2, 300, labels=c('G1', 'G2'))
itemstats(dat, group=group)


#####
# polytomous data example
itemstats(Science)

# polytomous data with missing
newScience <- Science
newScience[1:5,1] <- NA
itemstats(newScience)

# unequal categories
newScience[,1] <- ifelse(Science[,1] == 1, NA, Science[,1])
itemstats(newScience)

merged <- data.frame(LSAT7full[1:392,], Science)
itemstats(merged)

}
\references{
Chalmers, R., P. (2012). mirt: A Multidimensional Item Response Theory
Package for the R Environment. \emph{Journal of Statistical Software, 48}(6), 1-29.
\doi{10.18637/jss.v048.i06}
}
\seealso{
\code{\link{empirical_plot}}
}
\author{
Phil Chalmers \email{rphilip.chalmers@gmail.com}
}
\keyword{data}
