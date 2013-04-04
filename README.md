#mirt

Multidimensional item response theory in R. 

##Description

Analysis of dichotomous and polytomous response data using latent
trait models under the Item Response Theory paradigm. Includes univariate
and multivariate one-, two-, three-, and four-parameter logistic models,
graded response models, rating scale graded response models, (generalized)
partial credit models, rating scale models, nominal models, multiple choice
models, and multivariate partially-compensatory models. Many of these models 
can be used in an exploratory or confirmatory manner with optional user defined 
constraints. Exploratory models can be estimated via quadrature or
stochastic methods, a generalized confirmatory bi-factor analysis is
included, and confirmatory models can be fit with a Metropolis-Hastings
Robbins-Monro algorithm which can include polynomial or product constructed
latent traits. Additionally, multiple group analysis may be performed for
unidimensional or multidimensional item response models for detecting
differential item functioning.

##Installing from source

It's recommended to use the development version of this package since it is more likely to be up to date 
than the version on CRAN. To install this package from source: 

1) Obtain recent gcc and g++ compilers. Windows users can install the
   [Rtools](http://cran.r-project.org/bin/windows/Rtools/) suite while Mac users will have to
   download the necessary tools from the [Xcode](https://itunes.apple.com/ca/app/xcode/id497799835?mt=12) suite and its
   related command line tools (found within Xcode's Preference Pane under Downloads/Components); most Linux
   distributions should already have up to date compilers (or if not they can be updated easily). 

2) Install the `devtools` (if necessary). In R, paste the following into the console:

```r
install.packages('devtools')
```

3) Load the `devtools` package and install from the github source code. 
 
```r
library('devtools')
install_github('mirt', 'philchalmers', quick = TRUE)
```

##Installing from a binary (Windows only)

For those having difficulty installing the package on Windows, binary installation (.zip) files 
for 32- or 64-bit Windows may be installed with:

```r
download.file('http://dl.dropbox.com/u/10780530/mirt/mirt.zip', 'mirt.zip')
install.packages('mirt.zip', repos=NULL)
```

Note that this binary file is updated periodically and is not guarenteed to be in sync with 
the source code. 

#Extra

Bug reports are always welcome and the preferred way to address these bugs is through
the github 'issues'. Feel free to submit issues or feature requests on the site, and I'll 
address them ASAP. Cheers!
