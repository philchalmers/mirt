# Changes in mirt 1.44

- Standardized factor loadings for the multidimensional nominal response model
  now report consistent values regardless of the category ordering

# Changes in mirt 1.43

- `M2()` family no longer requires row-wise removal of missing data to behave
  correctly. As such, the `na.rm` argument has been removed as it is no longer
  required (requested by Ulrich Schroeders)

- Added support for latent regression ACOV/SE estimation with Oakes method
  in `mirt()`

- Related to both points below, general MLTM (Embretson, 1984) added when 
  itemtype is specified as `PC1PL` and an `itemdesign` set is used, where 
  formula must include the name of the factor in the formula expressions. See
  examples in the `mirt` documentation (requested by Susan Embretson)

- Added `PC1PL` itemtype to more easily specify conjunctive models with
  slopes fixed to 1 and estimation of the latent variance term, mimicking the 
  `Rasch` itemtype family

- `mirt()` and `multipleGroup()` gain `itemdesign` and `item.formula` arguments 
  to fit fixed item design characteristics (e.g. LLTMs; Fischer, 1983) 
  to all or a subset of items. Arguments are similar to those in `mixedmirt()`, 
  though currently not as flexible

- Partially-compensatory family of `itemtypes` now behave more consistently
  when loading structures specified where trace lines products are only 
  computed for dimensions with non-zero slopes

- `RCI()` gains a `shiny` logical to create an interactive scoring interface

# Changes in mirt 1.42

- The `model` argument in `bfactor()` can now be specified using the `mirt.model()`
  syntax to include more cognitively friendly tracking of item names and
  respective locations (requested by Afshin Khosravi)

- Add `reverse.score()` function for reverse scoring specific items within
  a `matrix` or `data.frame` 

- Fixed issue related to missing data patterns that resulted in bias 
  when estimating the hyper-parameters in single and multi-group models 
  (reported by Paul Jewsbury)

- `mirt.model()` syntax gains a negation operator for omitting specific 
  observed/latent groups from specifications. For example, the following will
  omit "Group3" identifies from between groups equality constraint definitions
  `CONSTRAINB[-Group3] = ...`

- `RMSD_DIF()` now supports datasets that follow vertical scaling structures
  (i.e., when groups answer some items but not others). 
  Requested by Alexandre Jaloto

- `M2()` functions now compute null model and SRMR fall all models whenever
  possible, including the latent class variance (reported by Hynek Cigler)
  
- VCOV memory leak bugfix for mixture models (see Github issue #247)

- Standardized residuals for point estimates now returned in `personfit()` 
  when passing `return.resids=TRUE` (requested by Raymond Hernandez)

# Changes in mirt 1.41

- Fix for `DIF()` when sparse data included with mixed item formats (reported by
  Heather Leigh Kayton)

- When computing category-level information curves include the negative Hessian
  in computations (reported by Milica Kabic)

- Allow missing data patterns in `personfit()`, as well as a new option 
  to return all raw item by person residuals (requested by George Karabatsos)

- Fix Zero-inflated model example in `multipleGroup()`, which required the 
  discontinuous trait location to be populated explicitly with a
  `customTheta` syntax (reported by Brooke Magnus)

- Empirical reliability estimates in `fscores()` and `empirical_rxx()` 
  include option to use the true score variance as an estimate of 
  the observed score variance (suggested by Hynek Cigler)

# Changes in mirt 1.40

- `technical` list gains a `nconstrain` argument for specifying equality 
  constraints with negative relationships (e.g., `a12 = -a21`). Requested by 
  Berend Terluin

- Added unipolar log‑logistic model (Lucke, 2015)  itemtype, specified 
  as `itemtype = 'ULL'`. Note that this automatically changes a number of 
  internal defaults, such as using a log-normal(0,1) density for the latent
  traits, and where the `theta_lim` is specified to be positive

- Added complementary log‑log model (Shim, Bonifay, and Wiedermann, 2022) 
  itemtype, specified as `itemtype = 'CLL'`

- Added `itemtype = '5PL'` model for unidimensional dichotomous data to included 
  asymmetric response functions. Example in `help(mirt)` also demonstrates 
  asymmetric 2PL model as the 5PL itself is very unstable and requires 
  strong priors

- Methods using Quasi-Monte Carlo integration post-convergence were 
  not respecting correlated latent variable structures 
  (reported by George Kephart when using `M2()`)

- Bugfix for `fscores()` when supplying mixture models 
  that was introduced by changing previous classification default for 
  latent class models (reported by Karel Veldkamp)
  
- `residuals()` gains a `p.adjust` argument for FWE control

- `DRF()` gains a `DIF.cat` argument to compute statistics on a per-category
  basis when studying polytomous items
  
- `expected.test()` gains a `probs.only` logical to return probability 
  functions for each category (only used when `individual = TRUE`)
  
- Small bug fixes in C++ code that resulted in memory leaks

# Changes in mirt 1.39

- For models fit using `mdirt()` the `fscores()` EAP and EAPsum methods now 
  always returns classification probabilities as the default (reported 
  by Matthew Madison)

- `SIBTEST()` gains a `DIF` logical to perform DIF tests across `suspect_set`
  
- `DIF()` and `SIBTEST()` gain a `pairwise` logical input to perform pairwise 
  post-hoc comparisons for multi-group applications
  
- `DRF()` gains `groups2test` argument and friends for multi-group models

# Changes in mirt 1.38.1

- infit and outfit statistics can now be computed in `itemfit()` when missing 
  data are present (requested by Hanif on the mirt-package forum: 
  https://groups.google.com/g/mirt-package/c/_mA3YbMmbzM/m/CydOl-F4BQAJ?utm_medium=email&utm_source=footer)

- `coef(..., IRTpars=TRUE)` is now applied to multidimensional IRT models, 
  provided that the item contains simple structure (suggested by Sverre Ofstad)

- Fixed `match()` bug in `SIBTEST()` when total score is missing
  (reported by Ziying Li)

- `fscores(..., method ='EAPsum')` now supports returning the ACOV matrices, 
  matching the behaviour of the other estimators

- Store previously defined `customItems` and `customGroup` lists for use in
  secondary functions (e.g., `DIF()`, `boot.mirt()`, etc). 
  Reported by Nataly Beribisky

- Combining priors with equality constraints no longer uses multiple prior 
  definitions in the likelihood computations. Hence, constrained parameters
  are now treated as though they are a single parameter with only one prior 
  distribution (reported by Matthias von Davier in the context of 
  multiple-group models with between group item priors)

- Added a `groups2test` argument to `DIF()` to isolate individual grouping
  variable specification when using more than 2 groups

- Implicit argument 'invariance' stored in multiple-group objects
  now automatically used in `boot.mirt()` (previously had to be 
  manually passed)
  
- Bugfix when using `items2test` in DIF when input is a character vector
  (reported by @jbuncher)
  
- Bug fixes for multiple-group DIF testing with `DIF()` when using more 
  than two groups (reported by Ruben Neda and Davin Díaz García)

# Changes in mirt 1.37.1

- `boot.mirt()` gains a `boot.fun` argument to accept user-defined functions for 
  extracting the associated statistics to bootstrap 

- When `verbose = TRUE` in `residuals()` a set of summary statistics is reported
  for easier flagging 

- `itemfit()` arguments changed to accommodate outputting tables more
  consistently. Now a single `return.tables` argument is used to specify
  which tables to return

- `anova()` removes support for the `verbose` flag, and instead labels 
  the rows of the resulting output to identify the models

- `X2` and `G2` classes of item-fit statistics now better deal with large 
  missing value vectors on a per-item basis for better consistency

- `technical` list gains a `storeEMhistory` flag to store the EM history
  (requested by @netique)

- `DRF()` gains best-fitting prior support (currently limited to Gaussian distributions)

- Correct index subset caused by tmp row removals in MG objects (fixes #227)

# Changes in mirt 1.36

- Progress bar added automatically (controlled via the `verbose` argument) 
  when using several of the package's secondary functions (e.g., `fscores()`, 
  `DIF()`, `'DRF()`, `mdirt()`, etc)

- Added `itemstats()` function to give basic item information statistics

- Item-EFA models now automatically flips negative signs in rotate solutions
  (e.g., via `summary()`) according to the sign of the largest observed 
  loading (allows easier interpretation of the resulting correlation matrix)

- `response.pattern` deals with completely missing vectors now (issue #220)

- `residuals()` gains a `approx.z` logical to transform LD values into 
  approximate z-ratios

- `mirt()`, `mixedmirt()`, and `multipleGroup()` now have `model = 1` to 
  fit a unidimensional IRT model by default

# Changes in mirt 1.35.1

- Added `covdata` argument to `fscores()` to allow latent regression covariate information as well. 
  Example added to `fscores()` documentation to demonstrate this addition

- Added `RCI()` function to compute reliable change index via IRT modelling

- Added delta method SE in `coef(., IRTpars = TRUE)` for the nominal and nested-logit models

- `itemfit()` gains a `S_X2.plot` argument to visualize the expected-observed probability 
  differences based on the S-X2 conditional sum-score strategy

- Added `type = 'EAPsum'` to `plot()` generic to view an expected vs observed sum-scores plot

- `itemfit()` gains a `p.adjust` argument to allow for p-value adjustments in the 
  output for all methods

- `anova()` generic now supports a `...` input to compare many nested models, compared in sequence

- Added `type = 'threshold'` to `itemplot()` to plot cumulative probability information (requested 
  by Azman Sami)
  
- Fixed Bug `Error in if (any(SEtmp < 0))` that appeared due new R 4.0+ behaviour (reported by 
  Ziying Li and Caroline Böhm)
  
- Fix bug in `itemfit()` when plotting multiple-group objects

- Bugfix in `fscores()` report on which row failed to converge when datsets contain response 
  patterns that were completely missing

# Changes in mirt 1.34

- Previous `technical = list(removeEmptyRows = TRUE)` input now deprecated. Response patterns
  that are now completely missing are supplied NA placeholders within estimation and 
  post-estimation supporting functions (e.g., `fscores()`, `personfit()`, `fixed()`, etc)

- Added `converged` element in `DIF()` output to evaluate whether the nested model iteration converged

- Added support for plausible-value draws in `fscores()` when using `response.pattern` argument

- Fix `SE.type = 'Fisher'` computation in multi-group models (reported by Felix Zimmer)

- Switch `par` and `f` inputs in `numerical_deriv()`

# Changes in mirt 1.33.1

- Added `gen.difficulty()` to compute the generalized difficulty statistics 
   described by Ali, Chang, and Anderson (2015) for polytomous response models (suggested by Alexander Freund)

- Added `RMSD_DIF()` to compute marginal effect size measure recently used in PISA anlayses when investigating 
  'badness-of-fit' DIF effects when using constrained multiple-group models

- `extract.group()` now explicitly requires the group name to be passed rather than the group number (this 
  is a far more natural route)

- `plot(..., type =)` now supports `'trace'`, `'infotrace'`, `'itemscore'`, and `''` for two-dimensional 
  models to create faceted graphics

- Added `read.mirt()` function back to package now that `plink` is again available on CRAN

- Syntax input from the `car` package's `lht()` function adopted within `mirt`'s `wald()` function
  for easier specifications (see examples)

- Better cope with syntax definitions of models in `DIF()`, particularly with the `CONSTRAINB` form
  (reported by Hao Wu)

- Corrected outer-product summation for `SE.type = 'Fisher'` computation (reported by Felix Zimmer)

- Added `fixedCalib()` function to perform the five fixed-calibration methods describe by Kim (2006)

- Empirical histogram `dentype` convergence tolerance no longer modified (default now the same as the Gaussian
  `dentype` criteria)

- Fix for GGUMs using model syntax input (was ignoring the slope loading specifications; reported by Ben Listyg)

- fixed `traditional2mirt()` math for gpcm when 5 or more category items are supplied (reported by Aiden Loe)

# Changes in mirt 1.32.1

- OpenMP support added to E-step portion of the package, where number of threads can be 
  specified via the `mirtCluster()` function argument `omp_threads`. Special thanks to 
  Matthias von Davier for providing the `omp reduction` code in the `Estep.cpp` file

- Behaviour of `mirt(..., large)` has now been modified, where `large = TRUE` now skips computing
  the unique response patterns for datasets that likely contain little to no repeated response patterns
  (suggested by Matthias von Davier). The previous two-step behaviour is now achieved by passing 
  `large = 'return'`, storing this list object, and passing it back to the `large` input argument 

- Positive/negative sign remove from chi-square components in `residuals(type = 'LD')` 
  (requested by Cengiz Zopluoglu to help avoid confusion). Sign is still however present in the
  standardized correlation estimates

- `itemtype = 'rsm'` reported the incorrect information functions due to use of - instead of +
  from `traditional2mirt()` (reported by Nasser Hasan)
  
- column names of the `fscores()` results now correspond to the model syntax definition names instead
  of the previous F# convention

- fix `method = 'classify'` option in `fscores()` when more than two mixtures are fitted
  (reported by Lisa Limeri)

- fix bug in `'drop_sequential'` scheme in `DIF()` introduced in the previous 
  version of mirt due to some internal organization changes (reported by Balal Izanloo)

- allow infit/outfit statistics to be computed for non-Rasch models (suggested by 
  Alexander Freund for use with GGUMs)

- added `p.adjust` argument to `DRF()` (requested by Keri J. S. Brady)

- support for computation of the ACOV matrix when the variance of the specific 
  factors are freely estimated in `bfactor()` 

- fix for `invariance = 'free_var'` argument in `multipleGroup()` for multidimensional 
  models with correlated traits, which previously fixed the correlation parameters 
  inadvertently (reported by Ruoyi Zhu)

- use proper `mins` internal when using `extract.group()` to keep the original minimum
  response scoring pattern (reported by Adam Ťápal)

- bugfix for single-group models for `draw_parameters()` (reported by Keri Brady 
  and @ddueber)

- numeric model specification in `bfactor()` bug patched when intervals were not 
  1 unit apart due to NA placeholders (reported by Luis Manuel Lozano)

- latent trait/class names now are forced to be different than the data column names 
  (bug reported by Nathan Carter)
  
- fixed `X2*_df` and `PV_Q1*` when missing data pattern resulted in dropped categories
  (reported by Mac Pank)

# Changes in mirt 1.31

- added `likert2int()` to convert Likert-type character/factor responses to integer data

- `estfun()` gains a `centering` argument to center the scores (contributed by Rudolf Debelak)

- `impute` argument in `itemfit()` and `M2()` have been deprecated in favour of removing data
  row-wise via `na.rm=TRUE`

- Acceptance ratio when using MH samplers now returned prior to 'Stage 2' during estimation so that 
  these ratios are better behaved. As well, an heuristic improved method for increasing/decreasing
  the acceptance ratios is now implemented

- Added `return_seq_model` to  `DIF()` to return the final MG model on the last iteration of the
  sequential search schemes

- Bugfix in `DIF()` when sequential scheme was selected but no items contained DIF on the first
  iteration (reported by Scott Withrow)

- `SIBTEST()` gains a `plot` argument to create various plots depicting the (weighted) differences 
  between the focal subtest versus the matched subtest information

- `residuals()` gains a `'JSI'` type to compute the JSI statistics proposed 
  by Edwards et al. (2018)

- `residuals()` gains an `'expfull'` type to compute an expected value table for 
  all possible response patterns (not just those observed in the data)
  
- Fix for `key` variable for nested-logit models when data are collapsed to have equal intervals 
  (reported by Emil Kirkegaard)
  
- Added delta method for IRT parameter transformations when using multiple-group models 
  (reported by Alex Miller)

# Changes in mirt 1.30

- `empirical.poly.collapse` argument added to `itemfit()` to plot expected score functions for polytomous
   items (suggested by Keri Brady)

- SRMSR now reported in `M2()` for GGUMs (suggested by Bo on the mirt-package forum)

- `weights` argument added to `estfun.AllModelClass` to allow for the inclusion
  of `survey.weights` to calculate the scores

- `DIF()` now simplifies the output by default rather than returning lists from `anova()`. Wald tests are always simplified

- Where applicable, RMSEA statistics are reported in `itemfit()` for tests that return suitable 
  X2 and df components
  
- Fix negative TLI and CFI values when using the C2 statistic from the `M2()` function (reported 
  by Jake Kraska and Charlie Iaconangelo)

- Fix delta method SEs for `'gpcm'` itemtype (reported by Lennart Schneider)

# Changes in mirt 1.29

- When lower/upper bounded parameters are included the default optimizer is now 'nlminb' rather 
  than 'L-BFGS-B'. This is mainly due to the instability in the 'L-BFGS-B' algorithm which 
  is prone to converging instantly for unknown reasons

- `mdirt()` gains a `item.Q` list to specify Q-matrices at the item-category level for each item
  
- `createItem()` functions gain an optional argument to the function definitions to allow for
  list-specified data from functions such as `mirt()` via a silent `mirt(..., customItemsData)`
  argument

- lattice `auto.key` default now reports lines rather than points. This is now more consistent
  when, for example, color theme is changed to black and white in the trellis window

- Added Differential Response Function (DRF) statistics from upcoming publication (Chalmers, accepted)
  in a new function entitled `DRF()`. These are related to compensatory and non-compensatory measures 
  of response bias for DIF, DBF, and DTF available from the SIBTEST framework but 
  for IRT model fitted within the multiple-group estimation framework

- `structure` argument added to `mdirt()` function to allow log-linear models for 
  simplifying the profile probability model computations
  
- export internally used `traditional2mirt()` function to transform a small selection of 
  classical IRT parameterizations into the slope-intercept form
  
- fix `survey.weights` input for multiple group models (reported by Leigh Allison)

- fix `itemtype = "rsm"` block restriction when items contain unequal category lengths 
  (reported by Aiden Loe)

- `SIBTEST()` computation of beta coefficient changed to match Shealy and Stout's (1993) 
  form of `p_k * (Y_R - Y_F)` (was previously `p_k * (Y_F - Y_R)`; reported by Craig Wells). 
  As well, `Jmin` default is increased to 5 to avoid conservative Type I error behavior in 
  longer tests
  
- Fix negative chi-square differences in `DIF()` function due to non-converged sub-models 
  (reported by Daniel McKelvey)

# Changes in mirt 1.28

- `M2()` function gains a `type` input to distinguish between the univariate-bivariate 
  collapsed M2* statistic and the bivariate only collapsed C2 statistic (Cai and Monro, 2014).
  C2 can be useful for polytomous items when there are too few degrees of freedom to compute
  the fully collapsed M2* 

- `multipleGroup()` gains the `dentype` argument to allow for mixture IRT models to be 
  fitted (e.g., `dentype = 'mixture-3'` fits a three-class mixture model). This also
  allow modifications such as the zero-inflated IRT model to be fitted

- `technical` gains a `zeroExtreme` logical flag to assign survey weights of 0 to extreme 
  response patterns (FALSE by default). This may be required when
  Woods' extrapolation-interpolation method is used with empirical histograms 
  to avoid ill defined extrapolated densities

- `fscores()`, `itemfit()`, `M2()`, and `residuals()` gain a `use_dentype_estimate` argument to 
  compute EAP-based scores whenever the latent trait density was estimated (e.g., via empirical histograms)

- Empirical histograms can now be now scaled to [0,1] using Woods' extrapolation-interpolation method
  via the input `dentype = 'empiricalhist_Woods'`. Degrees of freedom updated to reflect this change, 
  and 121 quadrature points are used instead of the previous 199 for better stability
  
- Semi-parametric Davidian curve estimation of the shape of the latent trait distribution in
  unidimensional IRT models was contributed by Oguzhan Ogreden, as well the associated 
  components used within this framework (such as the interpolation-extrapolation method 
  described by Woods, 2006). This estimation method is available through the new `dentype` input. 
  mirt also now links to the `dcurver` package to obtain the associated computation functions 
  in the EM algorithm
  
- `M2()`, `itemfit()`, `SIBTEST()`, and `fscores()` gain an `na.rm` logical to remove rows of missing data

- `fscores()` gains a `append_response.pattern` logical to indicate whether response patterns via the
  `response.pattern` input should be appended to the factor score results

- new `dentype` argument added to estimation-based functions to specify the density structure of the 
  latent traits (default is `'Gaussian'`). This update breaks the previous `empiricalhist` logical option 

- `anova()` will accept a single fitted model object and return information related to AIC, BIC, 
  log-likelihood, etc
  
- Hannan–Quinn (HQ) Criterion added to `anova()`

# Changes in mirt 1.27

- Added multidimensional version of sequential response model (e.g., Tutz, 1990). Includes 
  `itemtype = 'sequential'` for the multidimensional 2PL variant, and `itemtype = 'Tutz'`
  for the Rasch variant

- Printing IRT parameters via `coef(mod, IRTpars = TRUE)` now computes the delta method 
  for the `g` and `u` terms as well. Interpreting these is generally not recommended 
  due to their bounded parameter nature (CIs can be outside the range [0,1]), 
  but are included for posterity

- `createItem()` gains a `bytecompile` flag to indicate whether the internal functions should 
  be byte-compiled before using (default is TRUE)

- Special `GROUP` location holder in `mirt.model()` to index the group-level
  hyper-parameter terms

- `key2binary()` gains a `score_missing` flag to indicate whether missing values should be scored 
  as 0 or left as NA

- `createItem()` gains support for `derivType = 'symbolic'` and 
  `derivType.hss = 'symbolic'` to symbolically compute the gradient/Hessian 
  functions (template code-base contributed by Chen-Wei Liu)

- `createItem()` gains a `derivType.hss` argument to distinguish gradient from 
  Hessian numerical computations

- `mdirt()` gains support for `createItem()` inputs

- More plotting points added to default `plot()` and `itemplot()` generics to create 
  smoother traceline functions

## Bug fixes
  
- fix `simdata()` bug for new `ggum` itemtype

- fix new grouping syntax specification in `mirt.model()` when combining START and FIXED 
  (reported by Garron Gianopulos)
  
- fix `IRTpars = TRUE` input when itemtype was `Rasch` (reported by Benjamin Shear)

# Changes in mirt 1.26.3

- `mod2values()` and passing `pars = 'values'` now return `data.frame` objects without any
  factor variables (previously the defaults to `data.frame()` were used, which created factors
  for all categorical variables by default)

- Add `monopoly` itemtype to fit unidimensional monotonic polynomial item response model
  for polytomous data (see Falk and Cai, 2016)
  
- Add `ggum` itemtype to fit unidimensional/multidimensional graded unfolding model 
  (e.g., Roberts & Laughlin, 1996). Special thanks to David King for providing the 
  necessary C++ derivative functions and starting values

- Square brackets can now be included in the `mirt.model()` syntax to indicate
  group-specific constraints, priors, starting/fixed values, and so on. These
  are all of the general form `"CONSTRAIN [group1, group2] = ..."` or 
  `"FIXED [group1] = ..."`

- Added delta method for several classical IRT parameterization 
  (via `coef(model, IRTpars = TRUE)`) when a suitable information matrix 
  was previously estimated

- `numDeriv` dependency removed because `numerical_deriv()` now supports a local 
  Richardson extrapolation type. For best accuracy, this is now used as the default
  throughout the package

- `createItem()` and `lagrange()` now use Richardson extrapolation as default 
  instead of the less accurate forward/central difference method

- `estfun()` function added to extract gradient information directly from fitted objects 
  (contributed by Lennart Schneider)
  
- `simdata()` gains an `equal.K` argument to redraw data until $K$ categories are 
  populated for a given item
  
- Fix initialization of `fscores()` when using 'MH' plausible value imputations (reported by 
  Charlie Iaconangelo)

- Various other small bug fixes and performance improvements, fixes for Solaris compatibility,
  and run a small number of examples during R CMD check

# Changes in mirt 1.25

- `mdirt()` now supports latent regression covariate predictors. Associated function
  (e.g., `fscores()`) also include the latent regression information for discrete 
  models by default

- `SIBTEST()` replaced with the asymptotic sampling distribution version of 
   CSIBTEST described by Chalmers (accepted)

- `calcNull` set to `FALSE` by default

- Sandwich ACOV estimate now uses the Oakes estimate in the computations rather than the intensive 
  Louis form (which require low-level coding of the item-level Hessian terms). Added a 
  new `SE.type = 'sandwich.Louis'` for the original sandwich VCOV estimate in the previous version of
  mirt

- fix latent regression models with QMCEM and MCEM algorithms (reported by Seongho Bae)

- `fscores()` gains a `max_theta` argument to apply upper/lower bounds to iterative searching 
  algorithms (issue reported by Sebastian Born), and a `start` input to set the starting values
  as well (primarily useful in mirtCAT to reduce iterations)

- `alabama` package optimizer no longer used. Replaced with generic interface from `nloptr`
  package to support numerous optimizers with greater control instead. Associated inputs
  (e.g., `alabama_args`) replaced as well
  
- Export missing S4 methods for external R packages to import

# Changes in mirt 1.24

- MDIFF and MDISC no longer in normal ogive metric (1.702 scaling value removed)

- added `QMC` option to `residuals()` for `LD` and `LDG2` methods. Also globally set the number of QMC points to 
  5000 throughout the package for consistency

- `info_if_converged` and `logLik_if_converged` added to `technical` list to indicate whether the information matrix
  and stochastic log-likelihood should be computed only when the model converges. Default is now `TRUE` for both

- added `'MCEM'` method for Monte Carlo EM. An associated `MCEM_draws` function added to the `technical` list
  as well to control the number of draws throughout the EM cycles

- support for information matrix computations for QMCEM method added (e.g., Oakes, crossprod, Louis)

- globally improve numerical efficiency of QMC methods, including the QMCEM estimator

- include missing data values in `itemfit()` for parametric bootstrap methods to replicate missing
  data pattern

- ensure that nest-logit models have at least 3 categories (reported by Seongho Bae)

- convergence set to FALSE if any `g > u` is found in the 4PL model

- in verbose console output the log-posterior is printed when priors are included in the EM (previously was only the 
  marginal likelihood)
  
- various bug fixes to SIBTEST, particularly for very small sample sizes

# Changes in mirt 1.23

- `anova()` LRT comparison gains a `bounded` logical to indicate whether a bounded parameter is being compared, 
  as well as a `mix` argument to indicate the mixture of chi-squared distributions 

- MH-RM estimation `optimizer` argument can now be modified to `BFGS`, `L-BFGS-B`, and `NR` 
  instead of the default `NR1`

- a distinction between the `NR` optimizer in the EM and MH-RM applications is included, where the MH-RM now
  defaults to `NR1` to indicate a single Newton-Raphson update that uses an RM filtered Hessian term

- `method = 'SEM'` added to perform the stochastic EM algorithm (first two stages of the MH-RM algorithm setup).  
  Alternatively, setting `technical = list(NCYCLES = NA)` when using the MH-RM algorithm now returns 
  the stochastic EM results

- added `multidim_matrix` option to `iteminfo()` to expose computation of information matrices

- bounded parameter spaces handled better when using the NR optimizer

- various bug fixes and performance improvements

# Changes in mirt 1.22

- `SE.type = 'Oakes'` set as the new default when computing standard errors via the ACOV matrix when
  using the EM algorithm

- new `SE.type = 'Oakes'` to compute Oakes' 1999 form of the observed information matrix using 
  central difference approximation. Applicable for all IRT models (including customized IRT types)

- added support for `gpcmIRT` and `rsm` itemtypes for traditional generalized partial credit model and 
  Rasch rating scale model (which may be modified for a generalized rating scale model by freeing the 
  slope parameters)

- `SE.type = 'Fisher'` now supports the inclusion of latent distribution hyper-parameters. 
  Officially, all SE-types now provide proper hyper-parameter influence in the information matrices

- wrapped various output objects as `mirt_df`, `mirt_matrix`, and `mirt_list` class to 
  avoid the need for passing a `digits` argument for rounding output in the console. 
  Now, returned objects are never rounded, which makes writing Monte Carlo
  simulation code safer in that rounded results will not appear in the results

- added Stone's (2000) fit statistics and forthcoming PV-Q1 fit statistics to `itemfit()`

## Bug fixes

- patched underflow bug in `fscores()` when EAP estimates were used in extremely long (1000+ item) 
  tests. Error now reported when this happens. Using MAP estimates in these extreme situations
  is essentially equivalent and now recommended

# Changes in mirt 1.21

- add information about the number of freely estimated parameters to `print()` generic

- in `plot()`, `auto.key` is only disabled when `facet_items = FALSE` for dichotomous items. Also, adjusted 
  ordering of `plot(mod, type = 'itemscore')` to reflect actual item ordering in the data

- Stretched the theoretical bounds of the y-axis for score-based functions in `plot()` and `itemplot()`
  (e.g., 3PL models will now always stretch to S(theta) = 0)

- `plot(mod, type = 'score')` not supports the `which.items` input to make expected score plots for 
  bundles of items

- penalized term added to EM algorithm estimation subroutines to help keep the covariance matrix
  of the latent trait parameters positive definite in the M-step (helps convergence 
  properties of the optimizers, especially 'L-BFGS-B'). To turn this penalized 
  term off use `technical = list(keep_vcov_PD = FALSE)`

- added `type = 'itemscore'` to `plot()` generic to plot faceted version of the item
  scoring functions. Particularly useful when investigating DIF with `multipleGroup()`
  
- better support for `splines` itemtype in multiple-group models
  
## Bug fixes

- fix problem with 'EAPsum' in `fscores()` when `response.pattern` input 
  supplied (reported by Eva de Schipper)
  
- `plot(mod, type = 'rxx')` now uses the latent variance in the computations (reported by 
  Amin Mousavi)
  
- fix syntax input when customized IRT models are supplied

# Changes in mirt 1.20.1

- `df` adjustment for the `S_X2` item-fit statistic for models where the latent trait
  hyper-parameters have been estimated

- `itemfit()` and `personfit()` properly detect dichotomous Rasch models which 
  have been defined with the constrained slopes approach

- argument `'fit_stats'` now used in `itemfit()` to replace longer list of logicals 
  (e.g., `itemfit(mod, S_X2 = FALSE, X2 = TRUE, infit = FALSE, ...)`). Now fit stats
  are explicitly requested through a character vector input. Default still uses the 
  S_X2 statistic

- when using `'lnorm'` prior lower bound automatically set to 0, and with `'beta'` prior
  the lower and upper bounds are set to [0,1]

- `mdirt()` now uses `optimizer = 'nlminb'` by default

- revert using default 'penalized version of the BFGS algorithm' instead of L-BFGS-B 
  when box-constraints are used (introduced in version 1.19)

- Neale & Miller 1997 approximation added to `PLCI()` (default still computes exact PL CIs)

- `type = 'score'` supported for multiple group models in `itemplot()`

- added `poly2dich` function to quickly change polytomous response data to a comparable matrix of 
  dichotomous response data

# Changes in mirt 1.19

- a penalized version of the BFGS algorithm is now used instead of the L-BFGS-B when upper and lower
  bounds are included (provides more robust estimates)

- the variances of the orthogonal factors in `bfactor()` can now be freely
  estimated. This allows modeling of designs such as the testlet response model (example included in
  the documentation)

- new `spline` itemtype to model B-spline response functions for dichotomous models. Useful for 
  diagnostic purposes after detecting item-misfit. Additional arguments can be passed to the 
  `spline_args` list input to control the behaviour of the splines for each item. Currently limited
  to unidimensional models only

- `fscores()` gains a `plausible.type` argument to select between normal approximation PVs or 
  Metropolis-Hastings samples (suggested by Yang Liu)

- `mdirt()` has been modified to support DINA, DINO, located latent class, 
   and other diagnostic classification models. Additionally, the `customTheta` input required to build 
   customized latent class patterns has been changed from the previously cumbersome  
   `mdirt(..., technical = list(customTheta = Theta))` to simply `mdirt(..., customTheta = Theta)`

- `simdata()` gains a `prob.list` input to supply a list of matrices with probability values to be sampled
  from (useful when specialized response functions outside the package are required)

- `simdata()` supports 'lca' itemtypes for latent class model generation

- improved M2 accuracy when latent trait variances are estimated

- corrected behaviour of `M2()` when linear constraints are applied (M2 test was previously too conservative
  when constraints were used). This affects single as well as multiple-group models (reported by Rudolf Debelak)
  
- add plausible values for latent class and related models estimated from `mdirt()`

## BUG FIXES

- `multipleGroup()` throws proper error when vertical scaling is not identified correctly due to NAs

- S-X2 itemfit statistic fix when very rare expected categories appear (reported by Seongho Bae)

# Changes in mirt 1.18

- `mdirt()` function now includes explicit parameters for the latent class intercepts (in log-form).
  This implies that correct standard errors can be computed using various methods (e.g., SEM, Richardson,
  etc)

- new `customGroup()` function to define hyper-parameter objects for the latent trait distributions
  (generally assumed to be Gaussian with a mean and covariance structure)

- new `boot.LR()` function to perform a parametric bootstrap likelihood-ratio test between 
  nested models. Useful when testing nested models which contain bounded parameters (e.g.,
  testing a 3PL versus a 2PL model)

- adjust the `lagrange()` function to use the full information matrix (was previously only a 
  quasi-lagrange approximation)

- greatly improved speed in `simdata()`, consequently changes the default seed

## BUG FIXES

- fix crash error in `mirtmirt()` for multidimensional models with lr.random effects (reported 
  by Diah Wihardini)
  
- `expbeta` prior starting values fix by setting to the mean of the prior rather than the mode 
  (reported by Insu Paek)

# Changes in mirt 1.17.1

- `itemfit()` function reworked so that all statistics have their own input flag (e.g., `Zh = TRUE`,
  `infit = TRUE`, etc). Additionally, only S-X2 is computed by default and X2/G2 (and the associated
  graphics and tables) are computed using 10 fixed bins

- added `empirical.table` argument to return tables of expected/observed values for `X2` and `G2`

- `group.bins` and `group.fun` argument added to `itemfit()` to control the size of the 
  bins and the central tendancy function for `X2` and `G2` computations

- `'expbeta'` option added to implement a beta prior specifically for the `g` and `u` parameters which
  internally have been transformed to logits (performes the back transformation before computing the 
  values)

- check whether multiple-group models contain enough data to estimate parameters uniquely 
  when no constraints are applied 
  
- set the starting values the same when using `parprior` list or `mirt.model()` syntax
  (reported by Insu Paek)
  
- `empirical_ES()` function added for effect size estimates in DIF/DBF/DTF analyses (contributed by
  Adam Meade)

## BUG FIXES

- standardized loadings not correct when factor correlations included in confirmatory models 
  (reported by Seongho Bae)

- `MDISC` and `MDIFF` values were missing the 1.702 multiplicitive constant (reported by 
  Yi-Ling Cheng)

- fix information trace-lines in multiple-group plots (reported by Conal Monaghan)

- suppress standard errors in exploratory models when `rotate != 'none'` (suggested by Hao Wu)

- sequential schemes in `DIF()` generated the wrong results (reported by Adam Meade)

- `M2()` was not properly accounting for latent variance terms (reported by Ismail Cuhadar)

# Changes in mirt 1.16

- enable `lr.random` input to `mixedmirt()` for multilevel-IRT models which are not from the Rasch family

- add common `vcov()` and `logLik()` methods

- latent regression EM models now have standard error computation supporte with the 'complete', 
  'forward', 'central', and 'Richardson' methods

- new `areainfo()` function to compute the area under information curves within specified ranges 
  (suggested by Conal Monaghan)

- `method = 'BL'` supported for multiple-group models. As well, `SE.type = 'numerical'` included to return
  the observed-data ACOV matrix from the call to `optim()` (can only be used when the `BL` method is selected)

- new `SE.type = 'FMHRM'` to compute information matrix based on a fixed number of MHRM draws, and an
  associated `technical = list(MHRM_SE_draws)` argument has been added to control the number of draws

- added `lagrange` (i.e., score) test function for testing whether parameters should be freed in single
  and multiple group models estimated with the EM algorithm

- `numerical_deriv` function made available for simple numerical derivatives, which may be useful when
  defining fast custom itemtype derivative terms
  
- `SE.type` used to compute the ACOV matrix gained three numerical estimates for the 
  forward difference ('forward'), central difference ('central'), and Richardson extropolation ('Richardson')
  
## BUG FIXES

- SE methods based on the Louis (1982) computations no longer contain NA placeholders for the latent trait 
  hyper-parameters

# Changes in mirt 1.15

## MINOR CHANGES

- added SIBTEST and crossed-SIBTEST procedures with the new function `SIBTEST()`

- added `empirical_plot` function for building empirical plots (with potential smoothing) 
  when conditioning on the total score

- more low-level elements included in `extract.mirt()` function

- added `grsmIRT` itemtype for classical graded rating scale form (contributed by KwonHyun Kim)

- added missing analytic Hessian terms when `gpcm_mats` are used (contributed by Carl Falk)

## BUG FIXES

- fixed row-removal bug when using `technical = list(removeEmptyRows = TRUE)` (reported by Aaron Kaat)

# Changes in mirt 1.14

## MAJOR CHANGES

- the structure of the output objects now contains considerably fewer S4 slots, and instead are 
  organized into more structured list elements such as `Data`, `Model`, `Fit`, and so on. Additionally,
  the information matrix has slot has been removed in favour of providing the asymptotic covariance 
  matrix (a.k.a., the inverse of the information matrix)

## MINOR CHANGES

- added `extract.mirt()` function to allow more convenient extracting of internal elements

- `crossprod` SE.type now incorporates latent variable information (replaces NA placeholders)

- changed the default `full.scores = FALSE` argument to `TRUE` in `fscores()`

- added `profile` argument to `plot()` for `mdirt()` objects so that profile plots can be generated

- `converge_info` option added to `fscores()` to return convergence information

- add `removeEmptyRows` option to `technical` list

## BUG FIXES

- return a vector of `NA`s when WLE estimation has a Fisher information matrix with a determinant of 0 
  (reported by Christopher Gess)

- fix df in multiple-group models with crossed between/within constrains (reported by Leah Feuerstahler)

- compute residuals when responses are sparse, and return `NaN` when residual could not be computed 
  (reported by Aaron Kaat)

# Changes in mirt 1.13

- adjust plausible values format for multiple group objects

- `simdata()` gains a `model` input to impute data from pre-organized models (useful in conjunction
  with mirtCAT or to generate datasets from already converged models). Also gains a `mins` argument
  to specify what the lowest category should be for each item if `model` is not supplied (default is 0)

- number of `SEMCYCLES` increased from 50 to 100 in the MH-RM algorithm, and RM gain rate changed from
  `c(.15, .65)` to `c(.1, .75)` 

- further improve item fit statistics when using imputations

- facet plots now try to keep the items in their respective order

- panel theme for lattice plots changed from default to a lighter blue colour, and legend now
  automatically placed on the right hand side rather than the top
  
## BUG FIXES

- fix for Q3 computations (noticed by Katherine Castellano)

# Changes in mirt 1.10

- when using prior distributions, starting values now automatically set equal to the mode of 
  the prior distribution, and appropriate lower and upper parameter bounds are supplied

- added `NEXPLORE` term to `mirt.model()` to specify exploratory models via the syntax

- add `itemGAM()` function to provide a non-linear smoother for better understanding mis-functioning 
  items (and without loosing established precision by reverting to purely non-parametric IRT methods)

- category scores are now automatically recoded to have spaces of 1, and a message is printed if/when 
  this occurs

- added `MDISC()` and `MDIFF()` functions

- the inclusion of prior parameter distributions will now report the log-posterior rather than
  the log-likelihood. Functions such as `anova()` will also report Bayesian criteria rather than
  the previous likelihood-based model comparison statistics

- `impute` argument in `itemfit()` and `M2()` now use plausible values instead of point estimates

- `START` syntax element in `mirt.model()` now supports multiple parameters, and `FIXED` 
  argument added to declare parameters as 'fixed' at their staring values

- added `LBOUND` and `UBOUND` syntax support in `mirt.model()`

- report proper lower and upper bounds in starting values data frame and from `mod2values()`

- `invariance` argument to `bfactor()` now automatically indexes the second-tier factors to make 
  multiple-group testing with `bfactor()` easier

- remove `rotate` and `Target` arguments from model objects, and pass these only through axillary
  functions such as `summary()`, `fscores()`, etc

- `model` based arguments now can be strings, which are passed to `mirt.model()`. This is now
  the preferred method for defining models syntactically, though the previous methods 
  will still work

- integration range (`theta_lim`) globally set to `c(-6, 6)`, and number of default quadrature nodes
  have systematically increased in parameter estimation functions. This will slightly change some 
  numerical results, but provides more consistence throughout the package

- add `theta_lim` arguments to various functions 

- better control of QMC grid, and more effective usage for higher dimensions

- internal code organization now makes it easier to add user defined `itemtype`s (which can be natively 
  added into the package, if requested)

## BUG FIXES

- fix conservative imputation standard errors in `itemfit()` and `M2()` (reported by Irshad Mujawar)

- fixed plausible value draws for multidimensional latent regression models (reported by Tongyun Li)

- don't allow crossprod, Louis, or sandwich information matrices when using custom item types (reported by
  Charlie Rutgers)

# Changes in mirt 1.9

- when using `coef(mod, printSE=TRUE)` the `g` and `u` parameters are relabeled to `logit(g)` and
  `logit(u)` to represent the internal labels

- added various facet plots for three dimensional models to `plot()` generic

- support `optimizer = 'nlminb'`, and pass optimizer control arguments to a `contol` list

- added `fixef()` function to extract expected values implied by the fixed effect parameters
  in latent regression models

- added `gpcm_mats` argument to estimation functions for specifying a customize scoring pattern
  for multidimensional generalized partial credit models

- added `custom_theta` input to `fscores()` for including customized integration grids

- add a `suppress` argument to `residuals()` and `M2()` to suppress local dependence
  values less than this specific value
  
- print a message in `DIF()` and `DTF()` when hyper-parameters are not freely
  estimated in focal groups

- constraits for hetorogenous item names added to `mirt.model()` syntax

- WLE support for multidimensional models added

- added `'SEcontour'` argument to `plot()` generic

- use NA's in `fscores()` when response patterns contain all NA responses (suggested by Tomasz Zoltak)

## BUG FIX

- S-X2 in `itemfit()` now returns appropriate values for multiple-group models

- multidimensional plausible value imputation fix (reported by KK Sasa)

- `plot(..., type = 'infotrace')` for multiple group objects fixed (reported by Danilo Pereira)

# Changes in mirt 1.8

- `fscores()` nows accepts `method = "plausible"` to draw a single plausible value set

- `plot()` default type is now `score`, and will accept rotation arguments for exploratory models 
  (default rotation is `'none'`)

- `imputeMissing()` supports a list of plausible values to generate multiple complete datasets

- new `custom_den` input to `fscores()` to use custom prior density functions for Bayesian estimates

- more optimized version of the 'WLE' estimator in `fscores()`

- empirical reliability added when `method = 'EAPsum'` in `fscores()`

- new `START` argument in `mirt.model()` for specifying simple starting values one parameter at 
  a time

## BUG FIX

- fix carryover print-out error in `summary()` when confirmatory models were estimated

- bound contraints not were not included for group hyper-parameters (reported by KK Sasa)

# Changes in mirt 1.7

## MAJOR CHANGES

- improved estimation efficiency when using MH-RM algorithm. As a result, the default seed 
  was changed, therefore results from previous versions will be slightly different
  
- objects of class 'ExploratoryClass' and 'ConfirmatoryClass' have been merged into a single class
  'SingleGroupClass' with an `exploratory` logical slot
  
- the `technical = list(SEtol)` criteria for approximating the information matrix 
  was lowered to 1e-4 in `mixedmirt()` to provide better standard error estiamtes

## NEW FEATURES

- `boot.mirt` now uses the optimizer used to estimate the model (default previously was EM)

- `mixedmirt` now supports interaction effects in random intercepts, including cross-level
  interactions

- added `averageMI()` function to compute multiple imputation averages for the plausible
  values methodology using Rubin's 1987 method

- plausible value imputation now available in `fscores()` using the new `plausible.draws` 
  numeric input

- add `return.models` argument to `DIF()` to return estimated models with free/constrained 
  parameters

- latent regression models added to `mixedmirt()` for non-Rasch models using the 
  new `lr.formula` input

- `mirt.model()` syntax can now define within individual item equality constraints
  by using more than 1 parameter specification name in the syntax

- latent regression models added to `mirt()` function by using the new `covdata` and `formula` 
  inputs

- added confidence envelope plots to `PLCI.mirt`, and throw warnings when intervals could not be
  located

- `coef()` now accepts a `simplify` logical, indicating whether the items should be collapsed
  to a matrix and returned as a list of length 2 (suggested by Michael Friendly)
  
## BUG FIXES

- bias correction in variance estimates `mixedmirt` when random effects are 
  included (reported by KK Sasa)

- fix missing data imputation bug in `itemfit()` (reported by KK Sasa)

- M2 statistic for bifactor/two-tier models was overly conservative

- better checks for numerical underflow issues

- use triangle 0's for identifying exploratory IFA models.
  As such, standard errors/condition numbers for exploratory models can be estimated again

# Changes in mirt 1.6.1

## MAJOR CHANGES

- `sirt` package added to suggests list. Special thanks to Alexander Robitzsch (author of `sirt`) 
  for developing useful wrapper functions for mirt such as `mirt.wrapper.coef()`, `tam2mirt()`, and   
  `lavaan2mirt()`. As well, many examples in `sirt` demonstrate the possibility of estimating 
  specialized IRT models with `mirt`, such as the: Ramsay quotient, latent class,
  mixed Rasch, located latent class, probabilistic Guttman, nonparametric, discrete graded 
  membership, and multidimensional IRT discrete traits, DINA, and Rasch copula models. 
  
- exploratory IRT models are no longer rotated by default in `coef()`, and now requires an 
  explicit `rotate` argument
  
- computation of `S_X2` statistic in `itemfit` now much more stable for polytomous item types

- support for the `plink` package now unofficially dropped because it was removed from CRAN

- data inputs are now required to have category spacing codings exactly equal to 1 
  (e.g., [0, 1, 2, ...]; patterns such as [0, 2, 3] which are implicitly missing spaces 
  are now invalid)
  
## NEW FEATURES

- `mdirt` function added to model discrete latent variables such as latent class analysis for
  dichotomous and polytomous items. Can be used to model several other discrete IRT models as well,
  such as the located latent class model, multidimensional IRT with discrete traits, DINA models,
  etc. See the examples and documentation for details
  
- axillary support for `DiscreteClass` objects added to `itemfit()`, `M2()`, `fscores()`,
  `wald()`, and `boot.mirt()`
  
- the S-X2 statistic available in `itemfit()` has been generalized to include multidimensional
  models
  
- the method `'QMCEM'` has been added for quasi-Monte Carlo integration in `mirt()` 
  and `multipleGroup()` for estimating higher dimensional models with greater accuracy 
  (suggested by Alexander Robitzsch). Several axillary function such as `fscores()`, 
  `itemfit()`, and `M2()` also now contain an `QMC` argument (or will accept one through the ... 
  argument) to use the same integration scheme for better accuracy in higher dimensional models
  
- nonlinear parameter constraints for EM estimation can be specified by using the `Rsolnp` 
  and `alabama` packages by passing `optimizer = 'solnp'` and `optimizer = 'alabama'`, 
  as well as the relevant package arguments through the `solnp_ags` and `alabama_ags` list inputs
  
- `itemnames` argument added to `mirt.model()` to allow model specifications using raw 
  item names rather than location indicators
  
- `accelerate` argument changed from logical to character vector, now allowing three 
  potential options: 'Ramsay' (default), 'squarem', and 'none' for modifying the 
  EM acceleration approach

## BUG FIXES

- fixed bug in `bfactor()` starting values when NAs were specified in the `model` argument

- adjust overly optimistic termination criteria in EM algorithm

# Changes in mirt 1.5

## MAJOR CHANGES

- for efficiency, the Hessian is no longer computed in `fscores()` unless it is required in the 
  returned object
  
- estimation with `method = 'MHRM'` now requires and explicitly `SE=TRUE` call to compute the
  information matrix. The matrix is now computed using the ML estimates rather than 
  approximated sequentially after each iteration (very unstable), and therefore a 
  separate stage is performed. This provides much better accuracy in the computations
  
## NEW FEATURES

- new `extract.group()` function to extract a single group object from an objects 
  previously returned by `multipleGroup()` 

- return the SRMSR statistic in `M2()` along with the residual matrix (suggested by Dave Flora)

- accept `Etable` default input in `customPriorFun` (suggested by Alexander Robitzsch)

- vignette files for the package examples are now hosted on Github and can be accessed 
  by following the link mentioned in the vignette location in the index or `?mirt` help file

- E-step is now computed in parallel (if available) following a `mirtCluster()` definition

- run no M-step optimizations by passing `TOL = NaN`. Useful to have the model converge 
  instantly with all parameters exactly equal to the starting values

- confidence envelope plots in `itemplot()` generate shaded regions instead of dotted lines,
  and confidence interval plots added to `plot()` generic through the `MI` input
  
- passes to `fscores()` slightly more optimized for upcoming mirtCAT package release

- `method = 'EAPsum'` argument to `fscores()` support for multidimensional models
  
## BUG FIXES

- fix forcing all SEs MHRM information matrix computations to be positive

- `imputeMissing()` crash fix for multiple-group models

- fix divide-by-0 bug in the E-step when number of items is large

- fix crash in EM estimation with `SE.type = 'MHRM'` 
  
# Changes in mirt 1.4

## MAJOR CHANGES

- calculating the information matrix for exploratory item factor analysis models has been disabled
  since the rotational indeterminacy of the model results in improper parameter variation

- changed default `theta_lim` to `c(-6,6)` and number of quadrature defaults 
  increased as well

- `@Data` slot added for organizing data based arguments. Removed several data slots from
  estimated objects as a consequence

- removed 'Freq' column when passing a `response.pattern` argument to `fscores()`

- increase number of Mstep iterations proportionally in quasi-Newton algorithms 
  as the estimation approaches the ML location
  
- 'rsm' itemtype removed for now until optimized version is implemented
  
## NEW FEATURES

- link to `mirt` vignettes on Github have been registered with the `knitr` package and 
  are now available through the package index

- `optimizer` argument added to estimation function to switch the default optimizer. Multiple 
  optimizers are now available, including the BFGS (EM default), L-BFGS-B, 
  Newton-Raphson, Nelder-Mead, and SANN

- new `survey.weights` argument can be passed to parameter estimation functions (i.e., `mirt()`) to 
  apply so-called stratification/survey-weights during estimation

- `returnList` argument added to `simdata()` to return a list containing the S4 item objects, 
  Theta matrix, and simulated data

- support custom item type `fscores()` computations when `response.pattern` is passed instead
  of the original data

- `impute` option for `itemfit()` and `M2()` to estimate statistics via plausible imputation when 
  missing data are present

- multidimensional ideal-point models added for dichotomous items

- M2* statistic added for polytomous item types

- Bock and Lieberman (`'BL'`) method argument added (not recommend for serious use)

## BUG FIXES

- large bias correction in information matrix and standard errors for models that 
  contain equality constraints (standard errors were too high)

- drop dimensions fix for nested logit models

# Changes in mirt 1.3

## MAJOR CHANGES

- default `SE.type` changed to `crossprod` since it is better at detecting when models are not 
  identified compared to `SEM`, and is generally much cheaper to compute for larger models

- M-step optimizer now automatically selected to be 'BFGS' if there are no bounded parameters, and
  'L-BFGS-B' otherwise. Some models will have notably different parameter estimates because of 
  this, but should have nearly identical model log-likelihoods

- better shiny UI which adapts to the itemtype specifically, and allows for classical parameter
  inputs (special thanks to Jonathan Lehrfeld for providing code that inspired both these changes)

- scores.only option now set to `TRUE` in `fscores()`

- `type = 'score'` for plot generics no longer adjusts the categories for expected test scores

- M-step optimizer in EM now deters out-of-order graded response model intercepts (was a problem 
  if the startvalues were too far from the ML estimate in graded models)

## NEW FEATURES

- `return.acov` logical added to `fscores()` to return a list of matrices containing the 
  ACOV theta values used to compute the SEs (suggested by Shiyang Su)

- `printCI` logical option to `summary()` to print confidence intervals for standardized loadings

- new `expected.test()` function, which is an extension of `expected.item()` but for the whole test

- `mirt.model()` syntax supports multiple * combinations in `COV = ` for more easily specifying 
  covariation blocks between factors. Also allows variances to be freed by specifying the same
  factor name, e.g., `F*F`

- `full.scores.SE` logical option for `fscores()` to return standard errors for each respondent

- multiple imputation (MI) option in `fscores()`, useful for obtaining less biased factor score 
  estimates when model parameter variability is large (usually due to smaller sample size)

- group-level (i.e., means/covariances) equality constrains are now available for 
  the EM algorithm

- `theta_lim` input to `plot()`, `itemplot()`, and `fscores()` for modifying range of 
  latent values evaluated

## BUG FIXES

- `personfit()` crash for multipleGroup objects since itemtype slot was not filled (reported
  by Michael Hunter)

- fix crash in two-tier models when correlations are estimated (reported by David Wu)

- R 3.1.0 appears to evaluate List objects differently at the c level causing strange behaviour, 
  therefore slower R versions of some internal function (such as mirt:::reloadPars()) will be used
  until a patch is formed
  
- behaviour of `mvtnorm::dmvnorm` changed as of version 0.9-9999, causing widely different 
  convergence results. Similar versions of older mvtnorm functions are now implemented instead

# Changes in mirt 1.2.1

## MAJOR CHANGES

- `fitIndices()` replaced with `M2()` function, and currently limited to only dichotomous items
  of class 'dich'

- `bfactor()` default SE.type set to 'crossprod' rather than 'SEM'

- generalized partial credit models now display fixed scoring coefs

- `TOL` convergence criteria moved outside of the `technical` input to its own argument

- `restype` argument to `residuals()` changed to `type` to be more consistent with the package

- removed `fitted()` since `residuals(model, type = 'exp')` gives essentially the same output

- mixedmirt has `SE` set to `TRUE` by default to help construct a more accurate information matrix

- if not specified, S-EM `TOL` dropped to `1e-6` in the EM, and `SEtol = .001` for each
  parameter to better approximate the information matrix

## NEW FEATURES

- two new `SE.type` inputs: 'Louis' and 'sandwich' for computing Louis' 1982 computation of
  the observed information matrix, and for the sandwich estimate of the covariance matrix

- `as.data.frame` logical option for `coef()` to convert list to a row-stacked data.frame

- `type = 'scorecontour'` added to `plot()` for a contour plot with the expected total scores

- `type = 'infotrace'` added to `itemplot()` to plot trace lines and information on the same plot,
  and `type = 'tracecontour'` for a contour plot using trace lines (suggested by Armi Lantano)

- `mirt.model()` support for multi-line inputs

- new `type = 'LDG2'` input for `residuals()` to compute local dependence stat based on G2
  instead of X2, and `type = 'Q3'` added as well

- S-EM computation of the information matrix support for latent parameters, which previously
  was only effective when estimation item-level parameters. A technical option has also been
  added to force the information matrix to be symmetric (default is set to `TRUE` for better
  numerical stability)

- new `empirical.CI` argument in `itemfit()` used when plotting confidence intervals for
  dichotomous items (suggested by Okan Bulut)

- `printSE` argument can now be passed to `coef()` for printing the standard errors instead of
  confidence intervals. As a consequence, `rawug` is automatically set to `TRUE` (suggested
  by Olivia Bertelli)

- second-order test and condition number added to estimated objects when an information
  matrix is computed

- `tables` argument can be passed to `residuals()` to return all observed and expected tables
  used in computing the LD statistics

## BUG FIXES

- using `scores.only = TRUE` for multipleGroup objects returns the correct person ordering
  (reported by Mateusz Zoltak)

- `read.mirt()` crash fix for multiple group analyses objects (reported by Felix Hansen)

- updated math for `SE.type = 'crossprod'`

# Changes in mirt 1.1

## NEW FEATURES

- `facet_items` argument added to plot() to control whether separate plots should be constructed
  for each item or to merge them onto a single plot

- three dimensional models supported in `itemplot()` for types `trace`, `score`, `info`, and `SE`

- new DIF() function to quicky calculate common differential item functioning routines, similar
  to how IRTLRDIF worked. Supports likelihood ratio testings as well as the Wald approach, and includes
  forward and backword sequential DIF searching methods

- added a `shiny = TRUE` option to `itemplot()` to run the interactive shiny applet.
  Useful for instructive purposes, as well as understanding how the internal parameters of mirt behave

- `type = 'trace'` and `type = 'infotrace'` support added to `plot` generic for multiple group objects

- `fscores(..., method = 'EAPsum')` returns observed and expected values, along with general fit
  statistics that are printed to the console and returned as a 'fit' attribute

- removed multinomial constant in log-likelihood since it has no influence on nested model comparisons

- `SE.type = 'crossprod'` and `Fisher` added for computing the parameter information matrix based on the
  variance of the Fisher scoring vector and complete Fisher information matrix, respectively

- `customPriorFun` input to technical list now available for utilizing user defined prior distribution
  functions in the EM algorithm

- empirical histogram estimation now available in `mirt()` and `multipleGroup()` for unidimensional
  models. Additional plot `type = 'empiricalhist'` also added to the `plot()` generic

- re-implement `read.mirt()` with better consistency checking between the `plink` package

## BUG FIXES

- starting values for `multipleGroup()` now returns proper estimated parameter information from
  the `invariance` input argument

- remove `as.integer()` in MultipleGroup df slot

- pass proper item type when using custom pattern calles in `fscores()`

- return proper object in personfit when gpcm models used

# Changes in mirt 1.0

## NEW FEATURES

- `GenRandomPars` logical argument now supported in the `technical = list()` input. This will generate
  random starting values for freely estimated parameters, and can be helpful to determine if obtained
  solutions are local minimums

- seperate `free_var` and `free_cov` invariance options available in multipleGroup

- new `CONSTRAIN` and `CONSTRAINB` arguments in `mirt.model()` syntax for specifying equality
  constraints explicitly for parameters accross items and groups. Also the `PRIOR = ...`
  specification was brought back and uses a similar format as the new CONSTRAIN options

- `plot(..., type = 'trace')` now supports polytomous and dichotomous tracelines, and `type = 'infotrace'`
  has a better y-axis range

- removed the '1PL' itemtype since the name was too ambiguous. Still possible to obtain however by
  applying slope constraints to the 2PL/graded response models

- `plot()` contains a which.items argument to specify which items to plot in aggregate type, such as
  `'infotrace'` and `'trace'`

- `fitIndicies()` will return `CFI.M2` and `TLI.M2` if the argument `calcNull = TRUE` is passed. CFI stats also
  normed to fall between 0 and 1

- data.frame returned from `mod2values()` and `pars = 'values'` now contains a column indicating
  the internal item class

- use `ginv()` from MASS package to improve accuracy in `fitIndices()` calculation of M2

## BUG FIXES

- fix error thrown in `PLCI.mirt` when parameter value is equal to the bound

- fix the global df values, and restrict G2 statistic when tables are too sparse

# Changes in mirt 0.9.0

## NEW FEATURES

- `PLCI.mirt()` function added for computing profiled likelihood standard errors. Currently only applicable
  to unidimensional models

- prior distributions returned in the `pars = 'values'` data.frame along with the input parameters,
  and can be edited and returned as well

- full.scores option for `residuals()` to compute residuals for each row in the original data

- `bfactor()` can include an additional model argument for modeling two-tier structures introduced
  by Cai (2010), and now supports a `'group'` input for multiple group analyses

- added a general Ramsey (1975) acceleration to EM estimation by default. Can be disable with
  `accelerate = FALSE` (and is done so automatically when estimating SEM standard errors)

- renamed response.vector to response.pattern in `fscores()`, and now supports matrix input for
  computing factor scores on larger data sets (suggested by Felix Hansen)

- total.info logical added to `iteminfo()` to return either total item information or information
  from each category

- `mirt.model()` supports the so-called Q-matrix input format, along with a matrix input for the
  covariance terms

- MH-RM algorithm now accessible by passing `mirt(..., method = 'MHRM')`, and `confmirt()` function
  removed completely. `confmirt.model()` also renamed to `mirt.model()`

- support for polynomial and interaction terms in EM estimation

- lognormal priors may now be passed to parprior

- iterative computations in `fscores()` can now be run in parallel automatically following a
  `mirtCluster()` definition

- `mirtCluster()` function added to make utilizing parallel cores more convenient. Globally removed
  the cl argument for multi-core objects

- updated documentation for data sets by adding relevant examples, and added Bock1997 data set
  for replicating table 3 in van der Linden, W. J. & Hambleton, R. K. (1997)
  Handbook of modern item response theory

- general speed improvements in all functions

## BUG FIXES

- WLE estimation fixed and now estimates extreme response patterns

- exploratory starting values no longer crash in datasets with a huge number of NAs, which caused
  standard deviations to be zero

- math fix for beta priors

# Changes in mirt 0.8.0

## NEW FEATURES

- support for random effect predictors now available in `mixedmirt()`, along with a `randef()` function
  for computing MAP predictions for the random parameters

- EAPsum support in `fscores()` for mixed item types

- for consistency with current IRT software (rather than TESTFACT and POLYFACT), the scaling
  constant has been set to D = 1 and fixed at this value

- nominal.highlow option added to specify which categories are the highest and lowest in nominal models.
  Often provide better numerical stability when utilized. Default is still to use the highest and lowest
  categories

- increase number of draws in the Monte Carlo calculation of the log-likelihood from 3000 to 5000

- when itemtype all equal 'Rasch' or 'rsm' models the latent variance parameter(s) are automatically
  freed and estimated

- `mixedmirt()` more supportive of user defined R formulas, and now includes an internal 'items'
  argument to create the item design matrix used to estimate the intercepts. More closely mirrors
  the results from lme4 for Rasch models as well (special thanks to Kevin Joldersma for
  testing and debugging)

- `drop.zeros` option added to extract.item and itemplot to reduce dimensionality of factor structures
  that contain slopes equal to zero

- EM tolerance (TOL argument) default dropped to .0001 (originally .001)

- `type = 'score'` and `type = 'infoSE'` added to `plot()` generic for expected total score and joint test
  standard error/information

- custom latent mean and covariance matrix can be passed to `fscores()` for EAP, MAP, and EAPsum methods.
  Also applies to `personfit()` and `itemfit()` diagnostics

- scores.only option to `fscores()` for returning just the estimated factor scores

- bfactor can include NA values in the model to omit the estimation of specific factors for the
  corresponding item

## BUG FIXES

- limiting values in z.outfit and z.infit statistics for small sample sizes (fix suggested by Mike Linacre)

- missing data gradient bug fix in MH-RM for dichotomous item models

- global df fix for multidimensional confirmatory models

- SEM information matrix computed with more accuracy (M-step was not identical to original EM),
  and fixed when equality constrains are imposed

# Changes in mirt 0.7.0

## NEW FEATURES

- new `'#PLNRM'` models to fit Suh & Bolt (2010) nested logistic models

- `'large'` option added to estimation functions. Useful when the datasets being analysed are very
  large and organizing the data becomes a computationally burdensome task that should be avoided when
  fitting new models. Also, overall faster handling of datasets

- `plot()`, `fitted()`, and `residuals()` generic support added for MultipleGroup objects

- CFI and X2 model statistics added, and output now includes fit stats w.r.t. both G2 and X2

- z stats added for itemfit/personfit infit and outfit statistics

- supplemented EM ('SEM') added for calculating information matrix from EM history. By default
  the TOL value is dropped to help make the EM iterations longer and more stable. Supports
  parallel computing

- added return empirical reliability (`returnER`) option to `fscores()`

- `plot()` supports individual item information trace lines on the same graph (dichotomous items only) with
  the option `type = 'infotrace'`

- `createItem()` function available for defining item types that can be passed to estimation functions.
  This can be used to model items not available in the package (or anywhere for that matter) with the
  EM or MHRM. Derivatives are computed numerically by default using the numDeriv package for defining
  item types on the fly

- Mstep in EM moved to quasi-Newton instead of my home grown MV Newton-Raphson approach. Gives
  more stability during estimation when the Hessian is ill-conditioned, and will provide an easier
  front-end for defining user rolled IRT models

## BUG FIXES

- small bias fix in Hessian and gradients in `mirt()` implementation causing the likelihood to not always be
  increasing near maximum

- fix input to `itemplot()` when object is a list of model objects

- fixed implementation of infit and outfit Rasch statistics

- order of nominal category intercepts were sometimes backwards. Fixed now

- S_X2 collapsed cells too much and caused negative df

- `response.vector` input now supports NA inputs (reported by Neil Rubens)

# Changes in mirt 0.6.0

## NEW FEATURES

- S-X2 statistic computed automatically for unidimensional models via itemfit()

- EAP for sum-scores added to fscores() with method = 'EAPsum'. Works with full.scores option
  as well

- improve speed of estimation in multipleGroup() when latent means/variances are estimated

- multipleGroup(invariance = '') can include item names to specify which items are to
  be considered invariant across groups. Useful for anchoring and DIF testing

- type = 'trace' option added to plot() to display all item trace lines on a single
  graph (dichotomous items only)

- default estimation method in multipleGroup() switched to 'EM'

- boot.mirt() function added for computing bootstrapped standard errors with via the boot
  package (which supports parallel computing as well), as well as a new option SE.type = ''
  for choosing between Bock and Lieberman or MHRM type information matrix computations

- indexing items in itemplot, itemfit, and extract.item can be called using either a number or the
  original item name

- added probtrace() function for front end users to generate probability trace functions from models

- plotting item tracelines with only two categories now omits the lowest
  category (as is more common)

- parallel option passed to calcLogLik to compute Monte Carlo log-likelihood more quickly. Can also
  be passed down the call stack from confmirt, multipleGroup, and mixedmirt

- Confidence envelopes option added to itemplot() for trace lines and information plots

- lbound and ubound parameter bounds are now available to the user for restricting the parameter
  estimation space

- mod2values() function added to convert an estimated mirt model into the appropriate data.frame
  used to determine parameter estimation characteristics (starting values, group names, etc)

- added imputeMissing() function to impute missing values given an estimated mirt model. Useful
  for checking item and person fit diagnostics and obtaining overall model fit statistics

- allow for Rasch itemtype in multidimensional confirmatory models

- oblimin the new default exploratory rotation (suggested by Dave Flora)

- more flexible calculation of M2 statistic in fitIndicies(), with user prompt option if the internal
  variables grow too large and cause time/RAM problems

## BUG FIXES

- read.mirt() fixed when objects contain standard errors (didn't properly
  line up before)

- mixedmirt() fix when COV argument supplied (reported by Aaron Kaat)

- fix for multipleGroup when independent groups don't contain all potential response options
  (reported by Scot McNary)

- prevent only using 'free_means' and 'free_varcov' in multipleGroup since this would not be
  identified without further constraints (reported by Ken Beath)


# Changes in mirt 0.5.0


- all dichotomous, graded rating scale, (generalized) partial credit, rating scale, and nominal models
  have been better optimized

- wald() will now support information matrices that contain constrained parameters

- confmirt.model() can accept a string inputs, which may be useful for knitr/sweave documents
  since the scan() function tends to hang

- multipleGroup() now has the logical options bfactor = TRUE to use the dimensional reduction algorithm
  for when the factor pattern is structured like a bifactor model

- new fitIndices() function added to compute additional model fit statistics such as M2

- testinfo() function added for test information

- lower bound parameters under more stringent control during estimation and are bounded to never
  be higher than .6

- infit and outfit stats in itemfit() now work for Rasch partial credit and rating scale models

- Rasch rating scale models can now be estimated with potential rsm.blocks (same as grsm model).
  "Generalized" rating scale models can also be estimated, though this requires manipulating the
  starting values directly

- added sample size adjusted BIC (SABIC) information statistics

- new mixedmirt() function for estimating IRT models with person and item level (e.g., LLTM) covariates.
  Currently only supports fixed effect predictors, but random effect predictors are being developed

- more structured output when using the anova() generic

# Changes in mirt 0.4.2

- item probability functions now only permit permissible values, and models may converge even when
  the log-likelihood decreases during estimation. In the EM if the model does not have a strictly
  increasing log-likelihood then a warning message will be printed

- infit and outfit statistics are now only applicable to Rasch models (as they should be),
  and in itemfit/personfit() a 'method' argument has been added to specify which factor score
  estimates should be used

- read.mirt() re-added into the package to allow for translating estimated models into a format
  usable by the plink package

- test standard error added to plot() generic using type = 'SE', and expected score plot added
  to itemplot() using type = 'score'

- weighted likelihood estimation (WLE) factor scores now available (without standard errors)

- removed the allpars option to coef() generics and only return a named list with the (possibly
  rotated) item and group coefficients

- information functions slightly positively biased due to logistic constant adjustment,
  fixed for all models. Also, information functions are now available for almost all item response models
  (mcm items missing)

- constant (D) used in estimating logistic functions can now be modified (default is still 1.702)

- partcomp models recently broken, fixed now

- more than one parameter can now be passed to parprior to make specifying identical priors more
  convenient

# Changes in mirt 0.4.1

- relative efficiency plots added to itemplot(). Works directly for multipleGroup analysis
  and for comparing different item types (e.g., 1PL vs 2PL) can be wrapped into a named list

- infit and outfit statistics added to personfit() and itemfit()

- empirical reliability printed for each dimension when fscores(..., fulldata = FALSE) called

- better system to specify fixed/free parameters and starting values using
  pars = 'values'. Should allow for much better simulation based work

- graded model type rating scale added (Muraki, 1990) with optional estimation 'blocks'. Use
  itemtype = 'grsm', and the grsm.block option

- for multipleGroup(), optional input added to change the current freely estimated parameters
  to values of a previously computed model. This will save needless iterations in the EM and MHRM
  since these parameters should be much closer to the new ML estimates

- itemplot() supports multipleGroup objects now

- analytical derivatives much more stable, although some are not yet optimized

- estimation bug fix in bfactor(), and slight bias fix in mirt() estimation (introduced in version
  0.4.0 when multipleGroup() added)

- updated documentation and beamer slide show included for some background on MIRT and some
  of the packages capabilities

- labels added to coef() when standard errors not computed. Also allpars = TRUE is now the default

- kernel estimation moved entirely to one method. Much easier to maintain and guarantees consistency
  across methods (i.e., no more quasi-Newton algorithms used)

# Changes in mirt 0.4.0

- Added itemfit() and personfit() functions for uni and multidimensional models. Within itemfit
  empirical response curves can also be plotted for unidimensional models

- Wrapped itemplot() and fscores() into S3 function for better documentation. Also response curve
  now are all contained in individual plots

- Added free.start list option for all estimation functions. Allows a quicker way to
  specify free and fixed parameters

- Added iteminfo() and extract.item() to calculate the item information and extract
  desired items

- Multiple group estimation available with the multipleGroup() function. Uses the EM and MHRM
  as the estimation engines. The MHRM seems to be faster at two factors+
  though and naturally should be more accurate, therefore it is set as the default

- wald() function added for testing linear constraints. Useful in situations
  for testing sets of parameters rather than estimating a new model for a likelihood ratio test

- Methods that use the MHRM can now estimate the nominal, gpcm, mcm, and 4PL models

- fscores computable for multiple group objects and in general play nicer with missing data
  (reported by Judith Conijn). Also, using the options full.scores = TRUE has been optimized
  with Rcpp

- Oblique rotation bug fix for fscores and coef (reported by Pedro A. Barbetta)

- Added the item probability equations in the ?mirt documentation for reference

- General bug fixes as usual that were spawned from all the added features. Overall, stay frosty.

# Changes in mirt 0.3.1

- Individual classes now correspond to the type of methods: ExploratoryClass,
  ConfirmatoryClass, and MultipleGroupClass

- plot and itemplot now works for confmirt objects

- mirt can now make use of confmirt.model specified objects and hence be confirmatory as well

- stochastic estimation of factor scores removed entirely, now only quadrature based methods
  for all objects. Also, bfactor returned objects now will estimate all the factors scores instead
  of just the general dimension

- Standard errors for mirt now automatically calculated (borrowed from running a tweaked
  MHRM run)

# Changes in mirt 0.3.0

- radically changed the underlying mechanisms for the
  estimation functions and in doing so have decided that polymirt() was
  redundant and could be replaced completely by calling confmirt(data, number_of_factors). The
  reason for the change was to facilitate a wider range or MIRT models and to allow for easier
  extensions to future multiple group analysis and multilevel modelling

- new univariate and MV models are available, including the 1-4 parameter logistic
  generalized partial credit, nominal, and multiple choice models. These are called by specifying a character
  vector called 'itemtype' of length nitems with the options '2PL','3PL','4PL','graded','gpcm',
  'nominal', or 'mcm'; use 'PC2PL' and 'PC3PL' for partially-compensatory items. If itemtype = '1PL' or 'Rasch',
  then the 1-parameter logistic/1-parameter ordinal or Rasch/partial credit models are estimated for
  all the data. The default assumes that items are either '2PL' or 'graded', as before.

- flexible user defined linear equality restrictions may be imposed on all estimation functions,
  so too can prior parameter distributions, start values, and choice of which parameters to
  estimate. These all follow these general 2 steps:

    1) Call the function as you would normally would but use, for example,
       mirt(data, 1, startvalues = 'index') to return the start values as they are indexed
    2) Edit them as you please (without changing the structure), then input them back into
       the function as mirt(data, 1, startvalues = editedstartvalues).

  This is true for the parprior (MAP priors), constrain (linear equality constraints), and
  freepars (parameters freely estimated), each with their own little quirk. All inputs are lists
  with named parameters for easy identification and manipulation. Note that this means that
  the partial credit model and Rasch models may be calculated
  as well by modifying either the start values and constraints accordingly (e.g., constrain all
  slopes to be equal to 1/1.702 and not freely estimated for the classical Rasch model, or all equal
  but estimated for the 1PL model)

- number of confmirt.model() options decreased due to the new way to specify item types, startvalues, prior
  parameter distributions, and constraints

- plink package has not kept up with item information curves, so I'll implement my own for now.
  Replaced plink item plots from 'itemplots' function with ones that I rolled

- package descriptions and documentation updated

- coef() now prints slightly different output, with the new option 'allpars = TRUE' to display all
  the item and group parameters, returned as a list

- simdata() updated to support new item types

- more accurate standard errors for MAP and ML factor scores, and specific factors in bfactorClass
  objects can now be estimated for all methods

# Changes in mirt 0.2.6-1

- dropped the ball and had lots of bug fixes this round. Future
  commits will avoid this problem by utilizing the testthat package
  to test code extensively before release

- internal change in confmirt function to move MHRM engine outside the
  function for better maintenance

- theta_angle added to mirt and polymirt plots for changing the viewing angle
  w.r.t theta_1

- null model no longer calculated when missing data present

- fixed item slope models estimated in mirt() with associated standard errors

# Changes in mirt 0.2.6

- null model computed, allowing for model statistics such as TLI

- documentation changes

- many back end technical details about estimation moved to technical lists

- support for all GPArotation methods and options, including Target rotations

- polymirt() uses confmirt() estimation engine

- 4PL support for mirt() and bfactor(), treating the upper bound as fixed

- coef() now has a rotate option for returning rotated IRT parameters

# Changes in mirt 0.2.5

- Fixed translation bug in the C++ code from bfactor() causing illegal
  vector length throw

- Fixed fscores() bug when using polychotomous items for mirt() and
  bfactor()

- pass rotate='rotation' from mirt and polymirt to override default
  'varimax' rotation at estimation time (suggested by
  Niels Waller)

- RMSEA, G^2, and p set to NaN instead of internal placeholder when
  there are missing data

- df adjusted when missing data present

- oblique rotations return invisible factor correlation matrix

# Changes in mirt 0.2.4

- degrees of freedom correctly adjusted when using noncompensatory
  items

- confmirtClass reorganized to work with S4 methods, now work
  more consistently with methods.

- fixed G^2 and log-likelihood in logLik() when product terms included

- bugfix in drawThetas when noncompensatory items used

# Changes in mirt 0.2.3

- bugfixes for fscores, itemplot, and generic functions

- read.mirt() added for creating a suitable plink object

- mirt() and bfactor() can now accommodate polychotomous items using
  an ordinal IRT scheme

- itemplot() now makes use of the handy plink package plots, giving a
  good deal of flexibility.

- Generic plot()'s now use lattice plots extensively

# Changes in mirt 0.2.2

- Ported src code into Rcpp for future tweaking.

- Added better fitted() function when missing data exist (noticed by
  Erin Horn)

# Changes in mirt 0.2.1

- ML estimation of factor scores for mirt and bfactor

- RMSEA statistic added for all fitted models

- Nonlinear polynomial estimation specification for confmirt models, now
  with more consistent returned labels

- Provide better identification criteria for confmirt() (suggested by
  Hendrik Lohse)

# Changes in mirt 0.2.0

- parameter standard errors added for mirt() (1 factor only) and bfactor() models

- bfactor() values that are ommited are recoded to NA in summary and coef
  for better viewing

- 'technical' added for confmirt function, allowing for various tweaks and
  varying beta prior weights

- product relations added for confmirt.model(). Specified by enclosing in brackets
  and using an asterisk

- documentation fixes with roxygenize

# Changes in mirt 0.1.20

- allow lower bound beta priors to vary over items (suggested by James Lee)

# Changes in mirt 0.1.6

- bias fix for mirt() function (noticed by Pedro Barbetta)
